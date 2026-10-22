#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.14"
# dependencies = ["pysam>=0.24,<0.25", "cyvcf2>=0.34,<0.35"]
# ///
"""
audit_vcf.py - what did I actually get? Run this FIRST, before downloading anything.

    ./src/audit_vcf.py --vcf PATIENT_DATA/raw/<your_sample>.vcf.gz

The file name is arbitrary: --vcf is required and has no default, so nothing in these
scripts expects a particular one. The only hard requirement is that the index sits next to
the data with the same base name (<file>.vcf.gz plus <file>.vcf.gz.tbi) - tools find the
index by appending the suffix, so renaming one without the other breaks every region query.

Nothing here needs a reference file, so this runs while the 4.3 GB is still downloading.
annotate_vcf.py detects and fixes contig naming on its own; what this script is for is
everything that naming detection does not tell you and that later steps depend on:

  * symbolic alleles <DEL>, <DUP>, <NON_REF> - a structural call cannot be annotated by a
    coordinate join, and an already-called exon deletion would change the whole hunt;
  * phasing tags PGT, PID, PS, HP - their absence settles the phasing step before it is
    written: with none of them and a single sample, read-backed phasing from the CRAM or
    "phase undetermined" are the only remaining routes;
  * whether INFO/AC, AN and AF describe this child or a cohort this child was subset from;
  * whether the file is already annotated (CSQ, ANN, BCSQ) - if it is, VEP and its plugins
    come for free and the bcftools-versus-VEP question dissolves;
  * whether the file can be used at all: plain gzip instead of bgzip is a hard error, a
    missing index is simply built. The index delivered with the data is used as-is and its
    modification time is never examined - mtime cannot tell a stale index apart from the
    order the files were copied in. --reindex rebuilds it anyway if you want that.

Output: results/vcf_audit.json plus a readable summary on stderr.
"""

import argparse
import json
import sys
from pathlib import Path

import pysam
from cyvcf2 import VCF
from pysam import bcftools

#: Repository root; the script lives in src/, so defaults are anchored here, not to cwd.
ROOT = Path(__file__).resolve().parents[1]

PRIMARY = {str(i) for i in range(1, 23)} | {"X", "Y", "MT", "M"}
PHASE_TAGS = ("PGT", "PID", "PS", "HP")
ANNOTATION_TAGS = ("CSQ", "ANN", "BCSQ")


def log(msg):
    print(f"[audit] {msg}", file=sys.stderr, flush=True)


def warn(msg):
    print(f"[audit] WARNING: {msg}", file=sys.stderr, flush=True)


def die(msg):
    print(f"[audit] ERROR: {msg}", file=sys.stderr, flush=True)
    sys.exit(1)


def _json_default(obj):
    """cyvcf2 hands back numpy scalars and json does not know them.

    .item() converts a numpy scalar to the matching Python type, so the number survives as
    a number. Anything else is a genuine bug and should still raise.
    """
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {obj.__class__.__name__} is not JSON serializable")


def is_primary(contig):
    return contig.removeprefix("chr") in PRIMARY


def check_file(path, reindex=False):
    """The three index failures, all of which are silent in different ways."""
    if path.name.endswith((".tbi", ".csi")):
        die(f"{path.name} is an INDEX, not data. Pass the .vcf.gz itself - tools find the "
            "index by appending the suffix. The index is 2-3 MB where the data is 1-5 GB, "
            "so compare the file sizes if you are unsure which is which.")
    if not path.exists():
        die(f"{path} does not exist")
    with path.open("rb") as fh:
        if fh.read(4) != b"\x1f\x8b\x08\x04":
            die(f"{path.name} is not bgzip-compressed. Plain gzip looks identical from "
                "the outside (same .gz, same first two magic bytes) but tabix refuses it. "
                f"Fix: gunzip -c {path.name} | bgzip -c > fixed.vcf.gz")
    indexes = [Path(f"{path}{ext}") for ext in (".tbi", ".csi")]
    present = [i for i in indexes if i.exists()]
    if not present:
        log(f"{path.name} has no index - building one")
        pysam.tabix_index(str(path), preset="vcf", force=True)
        present = [Path(f"{path}.tbi")]
    elif reindex:
        # Index policy, matching annotate_vcf.py: the index that came with the data is
        # used as-is and its modification time is never examined, because mtime cannot tell
        # a genuinely stale index apart from the order the files happened to be copied in
        # (`cp` stamps each destination with the time of the copy). --reindex rebuilds it
        # on request, and that is the only way it gets rebuilt.
        log(f"{path.name}: rebuilding the index on request (minutes for a WGS)")
        pysam.tabix_index(str(path), preset="vcf", force=True)
    try:
        with pysam.VariantFile(str(path)) as handle:
            next(iter(handle), None)
    except (OSError, ValueError) as err:
        die(f"{path.name} cannot be read as a VCF ({err})")


def index_stats(path):
    """{contig: n_records} read from the index, or {} when the index cannot say.

    Indexes written by older tabix/bcftools carry no per-contig count metadata, and
    `bcftools index -s` refuses with "does not contain any count metadata" rather than
    returning zeros. That index is still perfectly usable for region queries - it only
    cannot answer "how many records" without reading the data - so this is reported and
    stepped over, not treated as a broken file.
    """
    stats = {}
    try:
        out = bcftools.index("-s", str(path)) or ""
    except Exception as err:                       # pysam raises SamtoolsError
        if "count metadata" not in str(err):
            raise
        warn(f"{Path(path).name}: the index predates per-contig count metadata, so record "
             "counts from the index are unavailable (region queries are unaffected). "
             "Re-run with --reindex to rebuild it with the bundled bcftools.")
        return stats
    for line in out.splitlines():
        fields = line.split("\t")
        if len(fields) >= 3:
            stats[fields[0]] = int(fields[2])
    return stats


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--outdir", type=Path, default=None,
                        help="where vcf_audit.json goes (default <repo>/results)")
    parser.add_argument("--max-examples", type=int, default=10)
    parser.add_argument("--regions", default=None,
                        help="limit the streaming pass to these contigs (e.g. 21) for a "
                             "quick end-to-end test; the header and index checks are "
                             "unaffected, but the counts then describe that region only")
    parser.add_argument("--reindex", action="store_true",
                        help="rebuild the tabix index instead of stopping when it is "
                             "missing or older than the data")
    args = parser.parse_args()

    vcf_path = args.vcf
    args.outdir = args.outdir or ROOT / "results"
    args.outdir.mkdir(parents=True, exist_ok=True)
    log(f"VCF: {vcf_path}  ({vcf_path.stat().st_size / 1024**3:.2f} GB)")
    check_file(vcf_path, reindex=args.reindex)

    reader = VCF(str(vcf_path))
    samples = list(reader.samples)
    info_keys = sorted({line.info().get("ID") for line in reader.header_iter()
                        if line.type == "INFO"} - {None})
    format_keys = sorted({line.info().get("ID") for line in reader.header_iter()
                          if line.type == "FORMAT"} - {None})
    contigs = list(reader.seqnames)
    # Which of them actually carry records. The header of a GATK VCF lists thousands of
    # contigs with nothing on them, and contig names sit in the core of the index rather
    # than in the optional count metadata that older tabix versions omit.
    with_data = []
    try:
        tbx = pysam.TabixFile(str(vcf_path))
        with_data = list(tbx.contigs)
        tbx.close()
    except (OSError, ValueError):
        pass
    style = "chr" if any(c.startswith("chr") for c in contigs) else "no_chr"
    non_primary = [c for c in contigs if not is_primary(c)]
    stats = index_stats(vcf_path)
    non_primary_records = sum(stats.get(c, 0) for c in non_primary) if stats else None
    total_from_index = sum(stats.values()) if stats else None

    log(f"samples ({len(samples)}): {', '.join(samples) or '-'}")
    log(f"contig naming: {style}")
    log(f"contigs in header: {len(contigs)}, non-primary: {len(non_primary)}")
    if with_data:
        np_with_data = [c for c in with_data if not is_primary(c)]
        log(f"contigs carrying records (from the index): {len(with_data)}, "
            f"of which non-primary: {len(np_with_data)}")
        if np_with_data:
            log(f"  {', '.join(np_with_data[:args.max_examples])}")
    if non_primary and non_primary_records is not None:
        log(f"  records on non-primary contigs: {non_primary_records:,} "
            f"({non_primary_records / max(total_from_index, 1):.2%})")
    elif non_primary:
        log("  records on non-primary contigs: unknown (no count metadata in the index)")
    if non_primary:
        log(f"  examples: {', '.join(non_primary[:args.max_examples])}")
        log("  -> annotate_vcf.py drops these before norm; `norm -f` aborts on a contig "
            "the FASTA does not have")
    log(f"INFO declared:   {', '.join(info_keys)}")
    log(f"FORMAT declared: {', '.join(format_keys)}")

    annotated = [k for k in ANNOTATION_TAGS if k in info_keys]
    phasing = [t for t in PHASE_TAGS if t in format_keys]
    log(f"existing annotation fields: {', '.join(annotated) or 'none (raw caller output)'}")
    if annotated:
        warn("this file is already annotated - VEP's output may be available for free, "
             "which changes the tool choice for the whole stage")
    log(f"phasing tags: {', '.join(phasing) or 'none'}")
    if not phasing and len(samples) == 1:
        log("  -> no phasing field and a single sample: phase cannot come from the VCF, "
            "and trio segregation is out. Read-backed phasing from the CRAM, or PM3 at "
            "supporting level instead of full.")

    # One streaming pass for what neither the header nor the index can answer.
    total = symbolic = multiallelic = indels = hets = hom_alt = 0
    an_values = {}
    symbolic_examples = []
    stream = reader(args.regions) if args.regions else reader
    for variant in stream:
        total += 1
        alts = variant.ALT
        if len(alts) > 1:
            multiallelic += 1
        if any(a.startswith("<") or a == "*" for a in alts):
            symbolic += 1
            if len(symbolic_examples) < args.max_examples:
                symbolic_examples.append(f"{variant.CHROM}:{variant.POS} {','.join(alts)}")
        if variant.is_indel:
            indels += 1
        if samples:
            gt = int(variant.gt_types[0])     # cyvcf2: 0 hom-ref, 1 het, 2 unknown, 3 hom-alt
            hets += gt == 1                   # int(), because gt_types is a numpy array and
            hom_alt += gt == 3                # int + numpy.bool_ silently becomes numpy.int64
        an = variant.INFO.get("AN")
        if an is not None:
            an_values[int(an)] = an_values.get(int(an), 0) + 1
        if total % 1_000_000 == 0:
            log(f"  ... {total:,} records")

    top_an = sorted(an_values.items(), key=lambda kv: -kv[1])[:5]
    if total_from_index is not None:
        log(f"records: {total:,} (index says {total_from_index:,})")
        if total != total_from_index:
            warn("the streaming count and the index disagree - the index does not describe "
                 "this file; rebuild it with --reindex before any region query")
    else:
        log(f"records: {total:,} (index carries no count metadata, nothing to compare)")
    log(f"  multiallelic: {multiallelic:,}   indels: {indels:,}   symbolic/star: {symbolic:,}")
    for example in symbolic_examples:
        log(f"    {example}")
    if symbolic:
        warn("symbolic or star alleles are present - these cannot be annotated by a "
             "coordinate join and need to be looked at separately")
    if samples:
        log(f"  genotypes: {hets:,} het, {hom_alt:,} hom-alt")
    log(f"  most common INFO/AN: {top_an}")
    if top_an and top_an[0][0] > 2:
        warn("AN > 2 dominates: this file was joint-genotyped across a cohort and then "
             "subset to one sample. INFO/AC, AN and AF then describe THAT COHORT, not "
             "this child, and must not be read as sample frequencies. The pipeline is "
             "unaffected because frequencies come from gnomAD, but the methods section "
             "has to say so.")

    audit = {
        "vcf": str(vcf_path),
        "regions": args.regions or "whole file",
        "contigs_with_data": with_data,
        "size_bytes": vcf_path.stat().st_size,
        "samples": samples,
        "contig_style": style,
        "n_contigs": len(contigs),
        "non_primary_contigs": non_primary,
        "non_primary_records": non_primary_records,
        "info_declared": info_keys,
        "format_declared": format_keys,
        "existing_annotation_fields": annotated,
        "phasing_tags": phasing,
        "n_records": total,
        "n_multiallelic": multiallelic,
        "n_indels": indels,
        "n_symbolic_or_star": symbolic,
        "symbolic_examples": symbolic_examples,
        "n_het": hets,
        "n_hom_alt": hom_alt,
        "an_top": top_an,
    }
    out = args.outdir / "vcf_audit.json"
    out.write_text(json.dumps(audit, indent=2, default=_json_default))
    log(f"written: {out}")
    log("next: annotate_vcf.py (needs the reference files from download_refs.py)")


if __name__ == "__main__":
    main()