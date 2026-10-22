#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.14"
# dependencies = []
# ///
"""Step 7 - turn the cascade output into a ranked, evidence-annotated candidate list.

Three jobs the earlier steps deliberately left undone:

  1. collapse per-variant cascade rows into CLAIMS - a compound-heterozygous pair is one
     claim about one gene, not two independent guesses;
  2. attach the evidence that decides classification but needs reference files the cascade
     does not open: gnomAD v4.1 (nhomalt, grpmax, faf95) and the NMD-escape check behind
     PVS1, read from the Ensembl GFF3 for the MANE Select transcript;
  3. order the claims and write them in the challenge's column layout.

What it does NOT do: classify. ACMG classification stays manual and lives in the report; this
script neither produces nor checks one.

Contig naming stays Ensembl ('15') throughout, as everywhere else in this repository. The
'chr' prefix the submission template wants is added at the single boundary where the
submission CSV is written, not here.

    ./src/rank_candidates.py
    ./src/rank_candidates.py --top 5 --proband-id <id-from-the-dataset>

Outputs (default <repo>/results):
    final_candidates.tsv     submission column layout, DRAFT epcr and notes to be reviewed
    finalists_evidence.tsv   one row per allele: gnomAD v4.1, NMD verdict, exon, ClinVar
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: Draft epcr by evidence tier. These are ORDERING weights, not calibrated probabilities,
#: and the file they land in is a draft for review - the challenge scores the ranking, so
#: what matters is that the separation is steep and the reasoning is written down.
TIERS = [
    ("A", 0.92, "two alleles in a checklist gene or a top-10 panel gene, at least one "
                "asserted P/LP in ClinVar with review criteria"),
    ("B", 0.50, "two alleles in a checklist gene or a top-10 panel gene, no ClinVar "
                "pathogenic assertion"),
    ("C", 0.15, "two alleles in a gene ranked in the top 100 of the phenotype panel"),
    ("D", 0.04, "two alleles in a gene present in the phenotype panel"),
    ("E", 0.02, "two alleles in a gene outside the phenotype panel"),
]
#: Secondary findings are, by definition, unlikely to explain the presenting phenotype;
#: epcr asks for exactly that probability, so they sit below every primary candidate.
SECONDARY_EPCR = 0.01

#: Panel rank at or above which a single-allele dominant claim is a PRIMARY candidate
#: rather than an incidental finding.
DOMINANT_PRIMARY_RANK = 100

LOF = {"stop_gained", "frameshift", "splice_acceptor", "splice_donor", "start_lost",
       "stop_lost", "transcript_ablation"}

SUBMISSION_COLUMNS = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
                      "chrom_2", "pos_2", "ref_2", "alt_2",
                      "epcr", "finding_type", "notes"]

GNOMAD_KEYS = ("AC", "AN", "AF", "nhomalt", "grpmax", "AF_grpmax", "faf95")


def log(msg):
    print(msg, flush=True)


def die(msg):
    sys.exit(f"ERROR: {msg}")


def write_atomic(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def read_tsv(path):
    rows = []
    with path.open(encoding="utf-8", errors="replace") as fh:
        header = next(fh).rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            f += [""] * (len(header) - len(f))
            rows.append(dict(zip(header, f)))
    return rows


def to_tsv(rows, columns):
    out = ["\t".join(columns)]
    for r in rows:
        out.append("\t".join(str(r.get(c, "")).replace("\t", " ").replace("\n", " ")
                             for c in columns))
    return "\n".join(out) + "\n"


# ------------------------------------------------------------------ reference lookups

def load_mane(path):
    """gene symbol -> MANE Select transcript id."""
    if not path.exists():
        log(f"  ! {path} missing - the NMD check will be skipped")
        return {}
    out = {}
    for i, line in enumerate(path.read_text().splitlines()):
        f = line.split("\t")
        if i == 0 and f[0] == "transcript":
            continue
        if len(f) >= 2:
            out[f[1].strip()] = f[0].strip()
    return out


def transcript_exons(gff, transcript, chrom):
    """(exons in TRANSCRIPTION order, strand) for one transcript, read from the GFF3.

    Transcription order, not coordinate order. On the minus strand the first exon has the
    HIGHEST coordinates, so sorting ascending numbers the exons backwards and puts the
    "last exon-exon junction" at the 5' end - which silently inverts the NMD verdict for
    every minus-strand gene.
    """
    exons, strand = [], None
    needle = f"transcript:{transcript}"
    with gzip.open(gff, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] != chrom or f[2] != "exon":
                continue
            if needle in f[8]:
                exons.append((int(f[3]), int(f[4])))
                strand = f[6]
    exons.sort(reverse=(strand == "-"))
    return exons, strand


def nmd_verdict(pos, exons, strand):
    """Does a premature stop at `pos` escape nonsense-mediated decay?

    The rule the ClinGen PVS1 decision tree uses: a premature termination codon in the
    last exon, or within ~50 nt of the last exon-exon junction, escapes NMD - the
    transcript survives and a truncated protein is made, so PVS1 is downgraded. Worth 4-6
    points of classification, which is why it is computed rather than assumed.

    `exons` must already be in transcription order; the junction is at the 3' edge of the
    penultimate exon, which is its END coordinate on the plus strand and its START on the
    minus strand.
    """
    if len(exons) < 2 or strand not in ("+", "-"):
        return "unknown", ""
    idx = next((i + 1 for i, (a, b) in enumerate(exons) if min(a, b) <= pos <= max(a, b)),
               None)
    if idx is None:
        return "position outside the transcript", ""
    detail = f"exon {idx}/{len(exons)} ({strand})"
    if idx == len(exons):
        return "escapes NMD (last exon)", detail
    junction = exons[-2][1] if strand == "+" else exons[-2][0]
    distance = (junction - pos) if strand == "+" else (pos - junction)
    if 0 <= distance <= 50:
        return f"escapes NMD ({distance} nt from the last junction)", detail
    if distance < 0:
        return f"downstream of the last junction ({-distance:,} bp)", detail
    return f"NMD-competent ({distance:,} bp upstream of the last junction)", detail


def gnomad_v41(paths, wanted):
    """{(chrom,pos,ref,alt): {source: {key: value}}} from the local v4.1 extracts."""
    found = defaultdict(dict)
    positions = {(c, p) for c, p, _, _ in wanted}
    for source, path in paths.items():
        if not path.exists():
            log(f"  ! {path} missing - v4.1 evidence will be incomplete for {source}")
            continue
        with gzip.open(path, "rt", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if (f[0], f[1]) not in positions:
                    continue
                info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
                for alt in f[4].split(","):
                    key = (f[0], f[1], f[3], alt)
                    if key in wanted:
                        found[key][source] = {k: info.get(k, "") for k in GNOMAD_KEYS}
    return found


def load_checklist(mva_tsv, extra_tsv):
    genes = set()
    for path in (mva_tsv, extra_tsv):
        if not path.exists():
            continue
        for line in path.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            g = line.split("\t")[0].strip()
            if g and g != "gene":
                genes.add(g)
    if not genes:
        die(f"no checklist genes read from {mva_tsv} or {extra_tsv}")
    return genes


# ------------------------------------------------------------------------- claims

def build_claims(cand, panel, checklist):
    """Collapse cascade rows into one claim per (gene, arm).

    A compound-heterozygous pair is a single claim about a single gene. Emitting the two
    alleles as separate guesses would spend two of the ten submission slots on one
    hypothesis and, under a metric that gives half credit for one allele of a pair, would
    also misrepresent what is being claimed.
    """
    by_key = defaultdict(list)
    for r in cand:
        by_key[(r["gene"], r["arm"])].append(r)

    claims = []
    for (gene, arm), rows in by_key.items():
        prank = int(panel[gene]["rank"]) if gene in panel else None
        on_list = gene in checklist

        if arm == "recessive":
            homs = [r for r in rows if r["pattern"] == "hom"]
            hets = [r for r in rows if r["pattern"] == "comp_het"]
            if homs:
                alleles, kind = homs[:1], "homozygous"
            elif len(hets) >= 2:
                # With more than two qualifying hets the pair is ambiguous; take the two
                # with the strongest evidence and say so in the notes.
                hets.sort(key=lambda r: (r["clnsig_class"] != "P/LP",
                                         r["consequence"] not in LOF, int(r["pos"])))
                alleles, kind = hets[:2], "compound heterozygous"
            else:
                continue
        else:
            alleles, kind = rows[:1], "single heterozygous (ClinVar P/LP)"

        plp = any(a["clnsig_class"] == "P/LP" and int(a["clnrev_stars"] or 0) >= 1
                  for a in alleles)
        claims.append({"gene": gene, "arm": arm, "kind": kind, "alleles": alleles,
                       "panel_rank": prank, "on_checklist": on_list, "plp": plp,
                       "n_qualifying": len(rows),
                       "many_alleles": alleles[0].get("many_alleles", "") == "yes"})
    return claims


def tier_of(claim):
    """Evidence tier, and therefore the draft epcr. Pure function of recorded evidence."""
    rank = claim["panel_rank"]
    top = claim["on_checklist"] or (rank is not None and rank <= 10)
    if claim["arm"] == "recessive":
        if top and claim["plp"]:
            return "A"
        if top:
            return "B"
        if rank is not None and rank <= 100:
            return "C"
        if rank is not None:
            return "D"
        return "E"
    # Dominant arm: the assertion IS the evidence, so it only reaches the top tier when
    # the gene also fits the phenotype.
    if top and claim["plp"]:
        return "A"
    if rank is not None and rank <= DOMINANT_PRIMARY_RANK:
        return "C"
    return "E"


def build_notes(claim, evidence, covered):
    """One sentence of rationale per claim, assembled from recorded evidence only."""
    g = claim["gene"]
    bits = [f"{g} {claim['kind']}."]
    for a in claim["alleles"]:
        key = (a["chrom"], a["pos"], a["ref"], a["alt"])
        ev = evidence.get(key, {})
        part = [f"{a['chrom']}:{a['pos']}{a['ref']}>{a['alt']} {a['consequence']}"]
        if a["aa_change"]:
            part.append(f"({a['aa_change']}, numbering from bcftools csq - verify on MANE)")
        if a["clnsig_class"]:
            part.append(f"ClinVar {a['clnsig_class']}, {a['clnrev_stars']} star(s)")
        v41 = ev.get("gnomad", {})
        if not covered:
            # The v4.1 files in refs/ are REGIONAL extracts over the checklist genes only.
            # Reporting "absent from gnomAD" for a gene outside them would manufacture
            # rarity evidence out of a lookup that never happened.
            part.append("outside the local gnomAD v4.1 extract (checklist genes only), "
                        "so no v4.1 frequency was looked up")
        else:
            for source in ("exomes", "genomes"):
                if source in v41:
                    d = v41[source]
                    part.append(f"gnomAD v4.1 {source} AC={d['AC']} AF={d['AF']} "
                                f"nhomalt={d['nhomalt']}")
                else:
                    part.append(f"absent from gnomAD v4.1 {source}")
        if ev.get("nmd"):
            part.append(f"{ev['nmd']}{', ' + ev['nmd_detail'] if ev['nmd_detail'] else ''}")
        bits.append(" ".join(part) + ".")
    if claim["panel_rank"] is not None:
        bits.append(f"Phenotype panel rank {claim['panel_rank']}"
                    f"{' (MVA checklist gene)' if claim['on_checklist'] else ''}.")
    else:
        bits.append("Not present in the phenotype panel.")
    if claim["kind"].startswith("compound"):
        gap = claim["alleles"][0].get("pair_distance_bp", "")
        if gap:
            bits.append(f"Phase not established; alleles {int(gap):,} bp apart"
                        + (" - within reach of read-backed phasing."
                           if int(gap) <= 1000 else
                           " - beyond short-read pairs, so PM3 stays supporting."))
    if claim["n_qualifying"] > len(claim["alleles"]):
        bits.append(f"{claim['n_qualifying']} qualifying variants in this gene; the two "
                    f"with the strongest evidence are reported.")
    if claim["many_alleles"]:
        bits.append("Gene flagged many_alleles - count may reflect gene size or "
                    "mappability rather than evidence.")
    return " ".join(bits)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default=ROOT / "results", type=Path)
    ap.add_argument("--refs", default=ROOT / "refs", type=Path)
    ap.add_argument("--extra-genes", default=ROOT / "config" / "mva_genes_extra.tsv",
                    type=Path)
    ap.add_argument("--gff", default=None, type=Path,
                    help="Ensembl GFF3 (default: the only .chr.gff3.gz in refs/)")
    ap.add_argument("--top", type=int, default=5,
                    help="primary claims to emit (default 5); secondary findings are "
                         "added after them, within the challenge's ten-row cap")
    ap.add_argument("--proband-id", default="FILL_ME",
                    help="identifier from the dataset; left as FILL_ME on purpose, because "
                         "guessing it would produce a submission scored against nobody")
    ap.add_argument("--panel-mode", choices=["full", "no-panel", "agnostic"], default="full",
                    help="how much prior knowledge the ranking may use. 'full' uses the "
                         "phenotype panel and the MVA checklist; 'no-panel' discards the "
                         "patient's phenotype and keeps the disease-group checklist; "
                         "'agnostic' discards both and ranks on variant evidence alone. "
                         "Neither mode changes which variants are CANDIDATES - the cascade "
                         "has already run genome-wide - only the order they are reported in")
    args = ap.parse_args()

    cand_path = args.results / "candidates.tsv"
    panel_path = args.results / "panel.tsv"
    if not cand_path.exists():
        die(f"{cand_path} not found - run ./src/cascade.py first")
    # The panel is an ordering input, not a selection input, so it is required only by the
    # mode that actually reads it. Running without it is a supported configuration, not a
    # degraded one: if the answer survives --panel-mode agnostic it does not rest on a
    # phenotype prior, and that is worth being able to demonstrate rather than assert.
    if args.panel_mode == "full" and not panel_path.exists():
        die(f"{panel_path} not found - run ./src/build_panel.py first, "
            "or pass --panel-mode no-panel")
    gff = args.gff or next(iter(sorted(args.refs.glob("*.chr.gff3.gz"))), None)

    log("[1/4] inputs")
    cand = read_tsv(cand_path)
    panel = ({r["gene"]: r for r in read_tsv(panel_path)}
             if args.panel_mode == "full" and panel_path.exists() else {})
    checklist = (load_checklist(args.refs / "mva_genes.tsv", args.extra_genes)
                 if args.panel_mode != "agnostic" else set())
    mane = load_mane(args.refs / "mane.tsv")
    log(f"  panel mode: {args.panel_mode}")
    log(f"  {len(cand):,} cascade rows | {len(panel):,} panel genes | "
        f"{len(checklist)} checklist genes")

    log("[2/4] claims")
    claims = build_claims(cand, panel, checklist)
    for c in claims:
        c["tier"] = tier_of(c)
    order = {t: i for i, (t, _, _) in enumerate(TIERS)}

    def sort_key(c):
        """Tier, then phenotype rank, then the variant evidence itself.

        The evidence terms matter only where the phenotype rank cannot
        discriminate - which is every claim outside the panel (they all share the
        sentinel) and every claim in --panel-mode no-panel or agnostic. Without
        them the fallback is alphabetical by gene, which would make the agnostic
        ranking meaningless and hide exactly the robustness check it exists for.
        """
        n_lof = sum(1 for a in c["alleles"] if a["consequence"] in LOF)
        return (order[c["tier"]],
                c["panel_rank"] if c["panel_rank"] is not None else 10**9,
                not c["plp"],            # ClinVar P/LP first
                -n_lof,                  # then more loss-of-function alleles
                -len(c["alleles"]),      # then a pair before a single allele
                c["gene"])

    claims.sort(key=sort_key)
    log(f"  {len(claims):,} claims; tiers: " +
        ", ".join(f"{t}={sum(1 for c in claims if c['tier'] == t)}"
                  for t, _, _ in TIERS))

    # One claim per gene among the primaries. The arms are independent by design, so a
    # gene with a P/LP het and a second qualifying het produces BOTH a dominant
    # single-allele claim and a recessive pair - and the dominant one sorts first, putting
    # a single allele at the top of the submission where the pair belongs. Under a metric
    # that awards half credit for one allele of a true pair, that halves the score. The
    # pair is a superset of the single claim, so it wins.
    eligible = [c for c in claims if c["tier"] != "E" or c["arm"] == "recessive"]
    best_per_gene = {}
    for c in eligible:
        prev = best_per_gene.get(c["gene"])
        if prev is None or len(c["alleles"]) > len(prev["alleles"]):
            best_per_gene[c["gene"]] = c
    dropped = [c for c in eligible if best_per_gene.get(c["gene"]) is not c]
    for c in dropped:
        log(f"  collapsed into the {best_per_gene[c['gene']]['kind']} claim: "
            f"{c['gene']} {c['kind']}")
    primaries = [c for c in eligible if best_per_gene.get(c["gene"]) is c][:args.top]
    chosen_genes = {c["gene"] for c in primaries}
    secondaries = [c for c in claims
                   if c["arm"] == "dominant" and c["gene"] not in chosen_genes][:10 - len(primaries)]
    for c in primaries:
        c["finding_type"] = "primary"
    for c in secondaries:
        c["finding_type"] = "secondary"
    selected = primaries + secondaries
    log(f"  selected {len(primaries)} primary + {len(secondaries)} secondary")

    log("[3/4] evidence")
    alleles = [(a["chrom"], a["pos"], a["ref"], a["alt"])
               for c in selected for a in c["alleles"]]
    v41 = gnomad_v41({"exomes": args.refs / "gnomad_v4.1_exomes.mva.ensembl.vcf.gz",
                      "genomes": args.refs / "gnomad_v4.1_genomes.mva.ensembl.vcf.gz"},
                     set(alleles))
    evidence = {}
    exon_cache = {}
    for c in selected:
        for a in c["alleles"]:
            key = (a["chrom"], a["pos"], a["ref"], a["alt"])
            ev = {"gnomad": v41.get(key, {}), "nmd": "", "nmd_detail": ""}
            if a["consequence"].split("&")[0] in LOF and gff and c["gene"] in mane:
                tx = mane[c["gene"]]
                if tx not in exon_cache:
                    exon_cache[tx] = transcript_exons(gff, tx, a["chrom"])
                ex, strand = exon_cache[tx]
                ev["nmd"], ev["nmd_detail"] = nmd_verdict(int(a["pos"]), ex, strand)
            evidence[key] = ev
    log(f"  gnomAD v4.1: {len(v41)}/{len(set(alleles))} alleles found | "
        f"NMD checked for {sum(1 for e in evidence.values() if e['nmd'])}")

    # A secondary finding goes to a clinician for action, so the bar is not "the cascade
    # kept it" but "an unambiguous loss of function in a gene where ONE such allele causes
    # disease". Three clauses, each of which does work on this data:
    #   truncating     - removes a missense and a 5'UTR allele that are carrier states
    #   NMD-competent  - removes a last-exon truncation that deletes almost nothing, which
    #                    is how the recurrent filaggrin null alleles present
    #   ClinVar P/LP with >=1 review star - already computed as claim["plp"]
    # The principled form of this test is "the gene's inheritance mechanism matches the
    # observed zygosity", which fails for a single heterozygous allele in a RECESSIVE
    # gene - a carrier state, not a finding. Implementing that needs per-gene inheritance
    # annotation the pipeline does not yet carry, so the clauses above stand in for it and
    # every rejection is logged with its reason rather than disappearing.
    kept_secondary = []
    for c in secondaries:
        a = c["alleles"][0]
        nmd = evidence.get((a["chrom"], a["pos"], a["ref"], a["alt"]), {}).get("nmd", "")
        why = []
        if a["consequence"].split("&")[0] not in LOF:
            why.append(f"not truncating ({a['consequence']})")
        if not nmd.startswith("NMD-competent"):
            why.append(f"NMD verdict: {nmd or 'not applicable'}")
        if not c["plp"]:
            why.append("no ClinVar P/LP with review criteria")
        if why:
            log(f"  not reported as a secondary finding: {c['gene']} - {'; '.join(why)}")
        else:
            kept_secondary.append(c)
    if len(kept_secondary) != len(secondaries):
        log(f"  secondary findings: {len(kept_secondary)} of {len(secondaries)} reported")
    selected = primaries + kept_secondary

    log("[4/4] writing")
    rows, ev_rows = [], []
    for i, c in enumerate(selected):
        base = {t: v for t, v, _ in [(x[0], x[1], x[2]) for x in TIERS]}
        epcr = (SECONDARY_EPCR if c["finding_type"] == "secondary"
                else base[c["tier"]] - 0.001 * i)       # strict order, no ties
        a1 = c["alleles"][0]
        a2 = c["alleles"][1] if len(c["alleles"]) > 1 else None
        rows.append({
            "proband_id": args.proband_id,
            "chrom_1": a1["chrom"], "pos_1": a1["pos"],
            "ref_1": a1["ref"], "alt_1": a1["alt"],
            "chrom_2": a2["chrom"] if a2 else "", "pos_2": a2["pos"] if a2 else "",
            "ref_2": a2["ref"] if a2 else "", "alt_2": a2["alt"] if a2 else "",
            "epcr": f"{epcr:.3f}", "finding_type": c["finding_type"],
            "notes": build_notes(c, evidence, c["on_checklist"]),
        })
        for a in c["alleles"]:
            key = (a["chrom"], a["pos"], a["ref"], a["alt"])
            ev = evidence[key]
            row = {"gene": c["gene"], "tier": c["tier"], "arm": c["arm"],
                   "finding_type": c["finding_type"],
                   "chrom": a["chrom"], "pos": a["pos"], "ref": a["ref"], "alt": a["alt"],
                   "consequence": a["consequence"], "aa_change_csq": a["aa_change"],
                   "clnsig_class": a["clnsig_class"], "clnrev_stars": a["clnrev_stars"],
                   "nmd": ev["nmd"], "nmd_exon": ev["nmd_detail"],
                   "v41_region_covered": "yes" if c["on_checklist"] else "no",
                   "panel_rank": c["panel_rank"] if c["panel_rank"] is not None else "",
                   "on_checklist": "yes" if c["on_checklist"] else "no"}
            for source in ("exomes", "genomes"):
                d = ev["gnomad"].get(source, {})
                for k in GNOMAD_KEYS:
                    row[f"v41_{source}_{k}"] = d.get(k, "")
            ev_rows.append(row)

    write_atomic(args.results / "final_candidates.tsv", to_tsv(rows, SUBMISSION_COLUMNS))
    ev_cols = list(ev_rows[0]) if ev_rows else []
    write_atomic(args.results / "finalists_evidence.tsv", to_tsv(ev_rows, ev_cols))
    write_atomic(args.results / "ranking.json", json.dumps({
        "built_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "panel_mode": args.panel_mode,
        "panel_genes_used": len(panel),
        "checklist_genes_used": len(checklist),
        "prior_knowledge_note": (
            "The phenotype panel and the MVA checklist are ORDERING inputs. Neither "
            "selects candidates: the cascade runs genome-wide and has already produced "
            "candidates.tsv before this step opens either file. --panel-mode controls how "
            "much prior knowledge the ordering may use, so that the dependence of the "
            "final answer on the phenotype can be measured rather than argued about."),
        "tiers": {t: {"epcr": v, "rule": rule} for t, v, rule in TIERS},
        "secondary_epcr": SECONDARY_EPCR,
        "dominant_primary_rank": DOMINANT_PRIMARY_RANK,
        "secondary_finding_rule": (
            "truncating consequence AND NMD-competent AND ClinVar P/LP with at least one "
            "review star. Stands in for the principled test - the gene's inheritance "
            "mechanism matching the observed zygosity - which needs per-gene inheritance "
            "annotation the pipeline does not carry. Rejections are logged with reasons."),
        "claims_total": len(claims),
        "selected": [{"gene": c["gene"], "arm": c["arm"], "tier": c["tier"],
                      "finding_type": c["finding_type"], "panel_rank": c["panel_rank"],
                      "finding_type_basis": (
                          "dominant arm, panel rank worse than "
                          f"{DOMINANT_PRIMARY_RANK} or absent"
                          if c["finding_type"] == "secondary" else
                          f"{c['arm']} arm, tier {c['tier']}"
                          + ("" if c["panel_rank"] is not None
                             else ", gene absent from the phenotype panel"))}
                     for c in selected],
        "caveats": [
            "epcr values are ordering weights derived from the evidence tier, not "
            "calibrated probabilities. The file is a DRAFT for manual review; the methods "
            "description must say so.",
            "Amino-acid numbering in notes comes from bcftools csq, which may pick a "
            "non-MANE transcript. Verify every reported change against the MANE Select "
            "transcript before submitting.",
            "Contig naming here is Ensembl. The 'chr' prefix is added only where the "
            "submission CSV is written.",
            "proband_id is left as FILL_ME unless given on the command line.",
            "The gnomAD v4.1 files are regional extracts over the checklist genes only. "
            "For any other gene the notes say the frequency was not looked up, rather "
            "than reporting it as absent - see v41_region_covered in the evidence file.",
            "finding_type is a DEFAULT, not a judgement. A recessive-arm claim is labelled "
            "primary whatever its tier, including a gene absent from the phenotype panel; "
            "a dominant-arm claim is primary only when the gene also fits the phenotype "
            "(panel rank at or above dominant_primary_rank). Read finding_type_basis above "
            "and relabel by hand where the referral indication says otherwise.",
            "Reported secondary findings are NOT ACMG SF genes. None of the candidates "
            "here appears on the ACMG SF v3.3 list (Genet Med 2025;27(8):101454), so none "
            "would be returned as a secondary finding under that standard; they are "
            "surfaced because the challenge asks for them. Say so in the report.",
        ],
    }, indent=2))

    log(f"  {args.results / 'final_candidates.tsv'} ({len(rows)} rows)")
    log(f"  {args.results / 'finalists_evidence.tsv'} ({len(ev_rows)} alleles)")
    log(f"  {args.results / 'ranking.json'}")
    for r, c in zip(rows, selected):
        log(f"    {r['epcr']}  {c['finding_type']:<9} {c['tier']}  {c['gene']:<10} "
            f"{c['kind']}")
    if args.proband_id == "FILL_ME":
        log("  ! proband_id is FILL_ME - pass --proband-id before building the submission")


if __name__ == "__main__":
    main()
