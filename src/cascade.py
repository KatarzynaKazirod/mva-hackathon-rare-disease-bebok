#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.14"
# dependencies = ["cyvcf2>=0.34,<0.35"]
# ///
"""Step 6 - genome-wide candidate cascade over the annotated VCF.

Two arms, because the phenotype ranking admits two inheritance models:

  recessive  a gene carrying a homozygous qualifying variant, or two heterozygous ones
  dominant   a single rare heterozygous variant already classified P/LP in ClinVar

The dominant arm is not decoration. The top of the phenotype panel can be a dominant,
usually de novo syndrome; a recessive-only cascade cannot see it at all, and without
parents there is no way to establish de novo - so the evidence has to come from the variant
itself being a known pathogenic allele.

Filter ORDER is load-bearing. The genotype quality gate runs FIRST: in a recessive search a
false heterozygous call destroys a real compound-heterozygous pair by burying the true second
allele among a hundred spurious ones, and a false homozygous call invents a pair that does
not exist. With one sample and no parents there is no other check on either.

    ./src/cascade.py --vcf PATIENT_DATA/work/annotated.vcf.gz

Outputs (default <repo>/results):
    candidates.tsv             one row per (variant, gene) that survived, with the arm
    cascade_counts.json        how many records survive each stage - feeds the funnel
    set_aside_checklist.tsv    variants in the checklist genes dropped ONLY by the
                               consequence filter, kept for step 9 (hidden second allele)
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from cyvcf2 import VCF

ROOT = Path(__file__).resolve().parents[1]

#: Genotype quality gate. A recessive search is far more sensitive to genotyping error than
#: a dominant one, because it reasons about PAIRS of calls in the same gene.
MIN_GQ = 20
MIN_DP = 10
HET_AB = (0.25, 0.75)

#: A "1/1" call with 55% alt reads is a miscalled heterozygote, or a deletion on the
#: other allele. The docstring says a false homozygous call invents a pair, so the
#: homozygous side needs its own floor, not just the heterozygous window.
HOM_MIN_AB = 0.85

#: Population frequency ceilings. Homozygotes are allowed a looser ceiling than
#: compound-heterozygous candidates because a recessive allele can sit at appreciable
#: frequency; what actually excludes it is the count of healthy homozygotes (nhomalt),
#: which the af-only resource does not carry - hence step 10.
MAX_AF_COMPHET = 0.001
MAX_AF_HOM = 0.005

#: Ceiling for a variant ClinVar already calls pathogenic. Founder alleles in
#: recessive disease sit well above the de-novo ceilings above, and dropping a known
#: pathogenic allele for being "too common" is how the real second hit goes missing.
MAX_AF_PLP = 0.02

#: Consequences kept. These are bcftools csq's own terms, not VEP's.
LOF = {"stop_gained", "frameshift", "splice_acceptor", "splice_donor",
       "start_lost", "stop_lost", "transcript_ablation"}
MISSENSE = {"missense"}
INFRAME = {"inframe_deletion", "inframe_insertion", "inframe_altering"}
KEEP = LOF | MISSENSE | INFRAME

#: Severity order for reporting only - it never decides what is kept.
SEVERITY = {c: i for i, c in enumerate(
    ["transcript_ablation", "splice_acceptor", "splice_donor", "stop_gained",
     "frameshift", "start_lost", "stop_lost", "inframe_deletion", "inframe_insertion",
     "inframe_altering", "missense", "splice_region", "start_retained",
     "stop_retained", "synonymous", "5_prime_utr", "3_prime_utr", "non_coding",
     "intron", "feature_elongation", "feature_truncation", "coding_sequence",
     "intergenic"])}

#: chrY and chrM have no place in an autosomal-recessive search, and a heterozygous
#: call on either is an artefact that can only pad the candidate list.
EXCLUDE_CONTIGS = {"Y", "chrY", "MT", "chrM", "chrMT"}

#: ClinVar tokens that count as pathogenic. Matching whole tokens rather than substrings is
#: deliberate: "Conflicting_classifications_of_pathogenicity" contains "pathogenic" and must
#: NOT be read as a pathogenic assertion.
CLNSIG_PATHOGENIC = {"pathogenic", "likely_pathogenic"}
CLNSIG_BENIGN = {"benign", "likely_benign"}


def log(msg):
    print(msg, flush=True)


def die(msg):
    sys.exit(f"ERROR: {msg}")


def write_atomic(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


# --------------------------------------------------------------------------- parsing

def first(value):
    """INFO fields declared Number=A or Number=. come back as a tuple from cyvcf2."""
    if isinstance(value, (tuple, list)):
        return value[0] if value else None
    return value


def parse_bcsq(raw):
    """INFO/BCSQ -> list of dicts.

    Format written by `bcftools csq -l`:
        consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change
    comma-separated across transcripts, with compound consequences joined by '&'.
    Entries for intergenic and regulatory records are shorter, and `-l` never emits the
    '@position' back-references that haplotype-aware mode produces.
    """
    out = []
    if not raw:
        return out
    for item in str(raw).split(","):
        if not item or item.startswith("@"):
            continue
        f = item.split("|")
        if len(f) < 2 or not f[1]:
            continue                      # no gene symbol - nothing to group on
        out.append({"consequence": f[0],
                    "terms": set(f[0].split("&")),
                    "gene": f[1],
                    "transcript": f[2] if len(f) > 2 else "",
                    "biotype": f[3] if len(f) > 3 else "",
                    "aa": f[5] if len(f) > 5 else "",
                    "dna": f[6] if len(f) > 6 else ""})
    return out


def clnsig_class(raw):
    """Return 'P/LP', 'B/LB', 'other' or '' from a ClinVar CLNSIG value."""
    v = first(raw)
    if not v:
        return ""
    tokens = {t.strip().lower()
              for part in str(v).replace("|", "/").replace(",", "/").split("/")
              for t in [part] if t.strip()}
    if tokens & CLNSIG_PATHOGENIC:
        return "P/LP"
    if tokens & CLNSIG_BENIGN:
        return "B/LB"
    return "other"


def clnrev_stars(raw):
    """ClinVar star rating from CLNREVSTAT. Decides how much an assertion is worth.

    The no-assertion values are tested FIRST because `no_assertion_criteria_provided`
    contains `criteria_provided` as a substring - the same trap as
    `Conflicting_classifications_of_pathogenicity` containing `pathogenic`.
    """
    v = str(first(raw) or "").lower()
    if v.startswith("no_assertion") or v.startswith("no_classification"):
        return 0
    if "practice_guideline" in v:
        return 4
    if "reviewed_by_expert_panel" in v:
        return 3
    if "multiple_submitters" in v and "no_conflicts" in v:
        return 2
    if "criteria_provided" in v:
        return 1
    return 0


def load_checklist(mva_tsv, extra_tsv):
    """Step 9 depends on this list, so an empty one is a failure, not a warning."""
    genes = set()
    for path, skip_header in ((mva_tsv, True), (extra_tsv, True)):
        if not path.exists():
            log(f"  ! {path} missing")
            continue
        for line in path.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            if skip_header and f[0].strip() == "gene":
                continue
            genes.add(f[0].strip())
    if not genes:
        die(f"no checklist genes read from {mva_tsv} or {extra_tsv} - step 9 depends "
            f"on this list, and an empty one makes a successful-looking empty result")
    return genes


# ------------------------------------------------------------------------- the pass

def scan(vcf_path, regions, checklist, counts, min_gq=MIN_GQ, min_dp=MIN_DP):
    """One pass over the annotated VCF. Returns (per_gene, set_aside).

    per_gene maps gene symbol -> list of row dicts that passed every per-variant filter.
    set_aside holds variants in checklist genes that passed everything EXCEPT the
    consequence filter, which is where a deep-intronic or UTR second allele would be.
    """
    vcf = VCF(str(vcf_path))
    if len(vcf.samples) != 1:
        die(f"expected a single-sample VCF, found {len(vcf.samples)}: {vcf.samples}")
    header = vcf.raw_header
    for field in ("BCSQ", "gnomAD_AF", "CLNSIG"):
        if f"##INFO=<ID={field}," not in header:
            die(f"INFO/{field} is not declared in the header - annotate the VCF first")

    per_gene = defaultdict(list)
    set_aside = []
    iterator = (vcf(r) for r in regions) if regions else (vcf,)

    for stream in iterator:
        for v in stream:
            counts["all records"] += 1
            # Everything below reads ALT[0], the first value of a Number=A field and a
            # single allele balance. On an unsplit 1/2 record all three are wrong with
            # no error at all, so the invariant is asserted rather than assumed.
            if len(v.ALT) > 1:
                die(f"multiallelic record at {v.CHROM}:{v.POS} - this step needs a VCF "
                    f"split with `bcftools norm -m -any` (annotate_vcf.py does it)")
            if v.CHROM in EXCLUDE_CONTIGS:
                counts["excluded: chrY / chrM"] += 1
                continue
            if counts["all records"] % 500_000 == 0:
                log(f"    {counts['all records']:,} records")

            # cyvcf2 returns None for FILTER when the record is PASS.
            if v.FILTER is not None:
                continue
            counts["FILTER = PASS"] += 1

            gt = v.gt_types[0]                      # 0 HOM_REF, 1 HET, 2 UNKNOWN, 3 HOM_ALT
            if gt not in (1, 3):
                continue
            counts["genotype non-reference"] += 1

            gq = float(v.gt_quals[0])
            dp = int(v.gt_depths[0])
            ref_d, alt_d = int(v.gt_ref_depths[0]), int(v.gt_alt_depths[0])
            have_ad = (ref_d >= 0 and alt_d >= 0 and ref_d + alt_d > 0)
            ab = alt_d / (ref_d + alt_d) if have_ad else float("nan")
            if gq < min_gq or dp < min_dp:
                continue
            if not have_ad:
                counts["dropped: no AD"] += 1
                continue
            if gt == 1 and not (HET_AB[0] <= ab <= HET_AB[1]):
                continue
            if gt == 3 and ab < HOM_MIN_AB:
                counts["dropped: hom with low allele balance"] += 1
                continue
            counts["quality gate"] += 1

            # ClinVar is consulted BEFORE the frequency and consequence gates, because
            # both of them can drop a known pathogenic allele: a founder allele for
            # being common, and a splice-region or deep-intronic one for its term.
            clnsig_raw = first(v.INFO.get("CLNSIG")) or ""
            cls = clnsig_class(clnsig_raw)
            stars = clnrev_stars(v.INFO.get("CLNREVSTAT"))
            plp, blb = cls == "P/LP", cls == "B/LB"

            af_raw = first(v.INFO.get("gnomAD_AF"))
            af = float(af_raw) if af_raw is not None else None
            # An annotator writing "." can reach cyvcf2 as nan, and nan > ceiling is
            # False, which would keep the variant while labelling it measured.
            if af is not None and math.isnan(af):
                af = None
            ceiling = MAX_AF_HOM if gt == 3 else MAX_AF_COMPHET
            # A missing frequency is "this resource has no entry", not "absent from the
            # population": af-only is a liftover of an older gnomAD release with uneven
            # non-coding coverage. Such variants are kept and FLAGGED, never zeroed.
            rescued = False
            if af is not None and af > ceiling:
                if plp and af <= MAX_AF_PLP:
                    rescued = True
                    counts["kept: ClinVar P/LP above the AF ceiling"] += 1
                else:
                    continue
            # A rare BENIGN missense completing a "pair" is a pure waste of attention.
            # Only act on an assertion that carries at least one review star.
            if blb and stars >= 1:
                counts["dropped: ClinVar B/LB"] += 1
                continue
            counts["population frequency"] += 1

            csq = parse_bcsq(first(v.INFO.get("BCSQ")))
            if not csq:
                continue

            row_base = {
                "chrom": v.CHROM, "pos": v.POS, "ref": v.REF,
                "alt": v.ALT[0] if v.ALT else "",
                "gt": "hom" if gt == 3 else "het",
                "gq": f"{gq:.0f}", "dp": str(dp), "ab": f"{ab:.3f}",
                "gnomad_af": "" if af is None else f"{af:.3g}",
                "af_status": "no_data" if af is None else "measured",
                "clnsig": str(clnsig_raw).replace("\t", " "),
                "clnsig_class": cls,
                "clnrev_stars": str(stars),
                "clnrevstat": str(first(v.INFO.get("CLNREVSTAT")) or ""),
                "rescued_by_clinvar": "af" if rescued else "",
                "chrom_note": "chrX" if v.CHROM in ("X", "chrX") else "",
                "clndn": str(first(v.INFO.get("CLNDN")) or "").replace("\t", " ")[:120],
            }

            # One variant hits several transcripts of the same gene. Group on the gene
            # symbol and keep the most severe kept consequence, otherwise a single variant
            # counted per transcript looks like a pair all by itself.
            by_gene = defaultdict(list)
            for c in csq:
                by_gene[c["gene"]].append(c)

            kept_any = False
            for gene, entries in by_gene.items():
                keepers = [c for c in entries if c["terms"] & KEEP]
                rescue = not keepers and plp
                if rescue:
                    keepers = entries           # a known pathogenic allele outranks
                    counts["kept: ClinVar P/LP outside the kept consequences"] += 1
                # A protein-coding transcript is REQUIRED, not merely preferred. bcftools
                # csq reports a consequence per transcript, and a gene whose only
                # qualifying hit is on an NMD, retained-intron or non-stop-decay transcript
                # has no protein-level consequence to report - the call is an annotation
                # artefact dressed as a missense. Preferring protein_coding (the old
                # behaviour) silently fell back to the artefact whenever no coding
                # transcript qualified: that is how a ribosomal-protein pseudogene locus
                # produced a "compound heterozygous pair of missense variants" in which one
                # allele was synonymous and neither touched MANE Select.
                coding = [c for c in keepers if c["biotype"] == "protein_coding"]
                if keepers and not coding and plp:
                    # A ClinVar P/LP assertion is about the variant, not about our choice
                    # of transcript, so it still outranks. Flagged, never silent.
                    coding = keepers
                    counts["kept: ClinVar P/LP on a non-coding transcript only"] += 1
                elif keepers and not coding:
                    counts["dropped: no protein-coding transcript"] += 1
                    if gene in checklist:
                        worst = min(keepers, key=lambda c: min(
                            (SEVERITY.get(t, 99) for t in c["terms"]), default=99))
                        set_aside.append({**row_base, "gene": gene,
                                          "consequence": worst["consequence"],
                                          "transcript": worst["transcript"],
                                          "biotype": worst["biotype"],
                                          "aa_change": worst["aa"],
                                          "dna_change": worst["dna"]})
                if coding:
                    best = min(coding, key=lambda c: min(
                        (SEVERITY.get(t, 99) for t in c["terms"]), default=99))
                    per_gene[gene].append({**row_base,
                                           "gene": gene,
                                           "rescued_by_clinvar": "consequence" if rescue
                                           else row_base["rescued_by_clinvar"],
                                           "consequence": best["consequence"],
                                           "transcript": best["transcript"],
                                           "biotype": best["biotype"],
                                           "aa_change": best["aa"],
                                           "dna_change": best["dna"]})
                    kept_any = True
                elif not keepers and gene in checklist:
                    # No qualifying consequence at all. `not keepers` is load-bearing: the
                    # keepers-but-no-coding-transcript case is set aside above, and without
                    # this guard a checklist gene would be written to set_aside twice.
                    # Most severe of the non-qualifying terms, not simply the first
                    # transcript listed - step 9 reads this column.
                    worst = min(entries, key=lambda c: min(
                        (SEVERITY.get(t, 99) for t in c["terms"]), default=99))
                    set_aside.append({**row_base,
                                      "gene": gene,
                                      "consequence": worst["consequence"],
                                      "transcript": worst["transcript"],
                                      "biotype": worst["biotype"],
                                      "aa_change": worst["aa"],
                                      "dna_change": worst["dna"]})
            if kept_any:
                counts["protein-altering consequence"] += 1

    vcf.close()
    return per_gene, set_aside


def assign_arms(per_gene, many_alleles=6):
    """Label each gene's surviving variants with the inheritance arm they qualify under.

    The two arms are evaluated INDEPENDENTLY, so one variant can appear in both. Treating
    them as mutually exclusive per gene was wrong and defeated the dominant arm entirely:
    a large gene almost always carries a second rare protein-altering heterozygous call, so
    a ClinVar P/LP allele sitting next to one became `comp_het` and lost both its label and
    its sort priority - which is exactly how the dominant syndrome at the top of the
    phenotype panel would go missing. Downstream steps must therefore key on
    (arm, chrom, pos, ref, alt), not on position alone.
    """
    rows = []
    for gene, variants in per_gene.items():
        homs = [v for v in variants if v["gt"] == "hom"]
        hets = [v for v in variants if v["gt"] == "het"]
        crowded = "yes" if len(variants) >= many_alleles else ""

        # --- recessive arm
        if homs:
            for v in homs:
                rows.append({**v, "arm": "recessive", "pattern": "hom",
                             "gene_n_alleles": len(variants), "pair_distance_bp": "",
                             "many_alleles": crowded})
            # A homozygous hit settles the gene; its heterozygous variants are reported
            # alongside rather than as a separate claim.
            for v in hets:
                rows.append({**v, "arm": "recessive", "pattern": "hom_gene_other_het",
                             "gene_n_alleles": len(variants), "pair_distance_bp": "",
                             "many_alleles": crowded})
        elif len(hets) >= 2:
            # Smallest gap between two qualifying heterozygous calls: under a few hundred
            # bp the pair is resolvable from read pairs in step 8, over tens of kb it is
            # not, so the number decides whether pulling a CRAM slice is worth it.
            positions = sorted(v["pos"] for v in hets)
            gap = min(b - a for a, b in zip(positions, positions[1:]))
            for v in hets:
                rows.append({**v, "arm": "recessive", "pattern": "comp_het",
                             "gene_n_alleles": len(hets), "pair_distance_bp": str(gap),
                             "many_alleles": crowded})

        # --- dominant arm, independent of everything above
        for v in hets:
            if v["clnsig_class"] == "P/LP":
                rows.append({**v, "arm": "dominant", "pattern": "single_het_clinvar_plp",
                             "gene_n_alleles": len(variants), "pair_distance_bp": "",
                             "many_alleles": crowded})
    return rows


COLUMNS = ["arm", "pattern", "gene", "on_checklist", "chrom", "pos", "ref", "alt",
           "gt", "consequence", "transcript", "biotype", "aa_change", "dna_change",
           "gq", "dp", "ab", "gnomad_af", "af_status",
           "clnsig_class", "clnrev_stars", "clnsig", "clnrevstat", "clndn",
           "gene_n_alleles", "pair_distance_bp", "rescued_by_clinvar",
           "many_alleles", "chrom_note"]


def to_tsv(rows, columns):
    out = ["\t".join(columns)]
    for r in rows:
        out.append("\t".join(str(r.get(c, "")) for c in columns))
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", default=ROOT / "PATIENT_DATA" / "work" / "annotated.vcf.gz",
                    type=Path, help="annotated VCF from ./src/annotate_vcf.py")
    ap.add_argument("--outdir", default=ROOT / "results", type=Path)
    ap.add_argument("--mva-genes", default=ROOT / "refs" / "mva_genes.tsv", type=Path)
    ap.add_argument("--extra-genes", default=ROOT / "config" / "mva_genes_extra.tsv",
                    type=Path)
    ap.add_argument("--regions", default=None,
                    help="comma-separated contigs, for a smoke test on one chromosome")
    ap.add_argument("--min-gq", type=int, default=MIN_GQ)
    ap.add_argument("--min-dp", type=int, default=MIN_DP)
    args = ap.parse_args()

    if not args.vcf.exists():
        die(f"{args.vcf} not found - run ./src/annotate_vcf.py first")
    regions = [r.strip() for r in args.regions.split(",")] if args.regions else None

    # Clear previous outputs BEFORE scanning. There is no cache here, so the risk is the
    # opposite of the annotation step's: all three files are written only at the end, so a
    # run that dies midway leaves the PREVIOUS run's files in place, and a chromosome-21
    # smoke test then sits in results/ looking like a genome-wide result. Only
    # cascade_counts.json records which it was. Removing them first makes a failed run
    # produce nothing instead of something misleading.
    OUTPUTS = ("candidates.tsv", "cascade_counts.json", "set_aside_checklist.tsv")
    stale = [args.outdir / f for f in OUTPUTS if (args.outdir / f).exists()]
    for f in stale:
        f.unlink()
    if stale:
        log(f"  removed {len(stale)} file(s) from a previous run")

    log("[1/3] inputs")
    checklist = load_checklist(args.mva_genes, args.extra_genes)
    log(f"  checklist genes: {len(checklist)}")
    log(f"  vcf: {args.vcf.name}" + (f" | regions: {', '.join(regions)}" if regions else ""))

    log("[2/3] scanning")
    counts = defaultdict(int)
    per_gene, set_aside = scan(args.vcf, regions, checklist, counts,
                              args.min_gq, args.min_dp)
    stages = ["all records", "FILTER = PASS", "genotype non-reference", "quality gate",
              "population frequency", "protein-altering consequence"]
    for s in stages:
        log(f"  {s:<32} {counts[s]:>12,}")
    for extra in ("excluded: chrY / chrM", "dropped: no AD",
                  "dropped: hom with low allele balance", "dropped: ClinVar B/LB",
                  "dropped: no protein-coding transcript",
                  "kept: ClinVar P/LP above the AF ceiling",
                  "kept: ClinVar P/LP outside the kept consequences",
                  "kept: ClinVar P/LP on a non-coding transcript only"):
        if counts[extra]:
            log(f"  {extra:<44} {counts[extra]:>8,}")

    log("[3/3] arms")
    rows = assign_arms(per_gene)
    for r in rows:
        r["on_checklist"] = "yes" if r["gene"] in checklist else "no"
    genes_rec = {r["gene"] for r in rows if r["arm"] == "recessive"}
    genes_dom = {r["gene"] for r in rows if r["arm"] == "dominant"}
    log(f"  recessive: {len(genes_rec):,} genes, "
        f"{sum(1 for r in rows if r['arm'] == 'recessive'):,} variants")
    log(f"  dominant:  {len(genes_dom):,} genes, "
        f"{sum(1 for r in rows if r['arm'] == 'dominant'):,} variants")
    hit = sorted((genes_rec | genes_dom) & checklist)
    log(f"  checklist genes with candidates: {', '.join(hit) if hit else 'none'}")

    # Expected scale, stated so a wrong answer is recognisable: a single WGS sample
    # normally leaves a few hundred genes with two qualifying alleles. A handful means the
    # filters are too tight; thousands mean the quality gate did not bite.
    if len(genes_rec) < 20:
        log("  ! fewer than 20 recessive genes - check the quality gate and AF ceilings")
    elif len(genes_rec) > 2000:
        log("  ! more than 2000 recessive genes - the quality gate is probably not biting")

    rows.sort(key=lambda r: (r["arm"] != "dominant", r["on_checklist"] != "yes",
                             -int(r["gene_n_alleles"]), r["gene"], r["pos"]))
    write_atomic(args.outdir / "candidates.tsv", to_tsv(rows, COLUMNS))
    write_atomic(args.outdir / "set_aside_checklist.tsv",
                 to_tsv(set_aside, [c for c in COLUMNS
                                    if c not in ("arm", "pattern", "on_checklist",
                                                 "gene_n_alleles", "pair_distance_bp")]))
    write_atomic(args.outdir / "cascade_counts.json", json.dumps({
        "built_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "regions": args.regions or "genome-wide",
        "stages": {s: counts[s] for s in stages},
        "arms": {
            "recessive": {"genes": len(genes_rec),
                          "variants": sum(1 for r in rows if r["arm"] == "recessive")},
            "dominant": {"genes": len(genes_dom),
                         "variants": sum(1 for r in rows if r["arm"] == "dominant")},
        },
        "checklist_genes_with_candidates": hit,
        "set_aside_checklist_variants": len(set_aside),
        "adjustments": {k: counts[k] for k in (
            "excluded: chrY / chrM", "dropped: no AD",
            "dropped: hom with low allele balance", "dropped: ClinVar B/LB",
            "dropped: no protein-coding transcript",
            "kept: ClinVar P/LP above the AF ceiling",
            "kept: ClinVar P/LP outside the kept consequences",
            "kept: ClinVar P/LP on a non-coding transcript only")},
        "parameters": {"min_gq": args.min_gq, "min_dp": args.min_dp,
                       "het_allele_balance": list(HET_AB),
                       "het_allele_balance_min_hom": HOM_MIN_AB,
                       "max_af_comphet": MAX_AF_COMPHET,
                       "max_af_hom": MAX_AF_HOM, "max_af_clinvar_plp": MAX_AF_PLP,
                       "consequences_kept": sorted(KEEP),
                       "require_protein_coding_transcript": True},
        "caveats": [
            "A qualifying consequence must fall on a protein_coding transcript. A gene "
            "whose only hit is on an NMD, retained-intron or non-stop-decay transcript is "
            "dropped, because bcftools csq reports per transcript and such a call has no "
            "protein-level consequence to report. The sole exception is a ClinVar P/LP "
            "assertion, which is about the variant rather than our transcript choice; "
            "those are counted separately.",
            "A missing gnomAD_AF is recorded as af_status=no_data and kept, not treated "
            "as absent from the population; the af-only resource is a liftover of an "
            "older release with uneven non-coding coverage.",
            "Phase is unknown: comp_het means two qualifying heterozygous calls in one "
            "gene, not two alleles proven to be in trans. pair_distance_bp says whether "
            "read-backed phasing from a CRAM slice could resolve it.",
            "The dominant arm depends entirely on ClinVar having already classified the "
            "variant, because de novo status cannot be established without the parents; "
            "clnrev_stars says how much that assertion is worth, and a 0-star call "
            "carries no review criteria at all.",
            "The arms are independent, so one variant can appear in both. Key on "
            "(arm, chrom, pos, ref, alt), never on position alone.",
            "Sex is not inferred. chrX rows are flagged in chrom_note; in a male, "
            "heterozygous calls outside the pseudoautosomal regions are artefacts and "
            "can form spurious comp_het pairs. chrY and chrM are excluded outright.",
            "many_alleles flags genes carrying six or more qualifying variants - the "
            "TTN / mucin / segmental-duplication class, where the count reflects gene "
            "size and mappability rather than evidence.",
            "Consequences come from bcftools csq, which carries no NMD-escape logic, so "
            "a stop_gained here is not yet evidence for PVS1 at full strength.",
        ],
    }, indent=2))

    log(f"  {args.outdir / 'candidates.tsv'}")
    log(f"  {args.outdir / 'cascade_counts.json'}")
    log(f"  {args.outdir / 'set_aside_checklist.tsv'} ({len(set_aside):,} variants)")


if __name__ == "__main__":
    main()
