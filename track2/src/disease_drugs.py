#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Ask the same question of the whole disease, not of one child's genotype.

Everything else in this repository starts from the two BUB1B alleles found in
one patient. That is the right place to start and the wrong place to stop: a
repurposing hypothesis that only applies to one compound heterozygote is not a
result anybody can build on. This script therefore asks the disease-level
question in two independent ways.

  1. Is any medicine recorded against mosaic variegated aneuploidy itself?
     Every MONDO entry in the group is queried - the grouped entry
     (Orphanet:1052) and each numbered subtype - and a disease with no drugs
     produces a row saying zero rather than being left out.

  2. Walking from the disease to its genetically associated targets, and from
     those targets to drugs, what is available across the group as a whole?
     This is the part that generalises past this patient: MVA1 is BUB1B, but
     MVA2 is CEP57, MVA3 is TRIP13, and Atelis syndrome is SLF2/SMC5.

A trap this walk contains, and why the script checks for it: the targets
associated with MVA include PAK6, which is not an MVA gene. PAK6 sits at
chr15:40,217,428-40,277,489, overlapping BUB1B, and a readthrough locus
BUB1B-PAK6 is annotated across both. The association is positional, not
biological. PAK6 happens to carry no drugs, so nothing propagates here - but a
pipeline that walked disease to target to drug without looking would have been
willing to propose a PAK inhibitor for MVA on the strength of a transcript
annotation. Overlap with a known gene is therefore reported per row.

Usage:
    python3 src/disease_drugs.py
    python3 src/disease_drugs.py --min-score 0.3 --out-dir /tmp/mva
"""
from __future__ import annotations

import argparse
import csv
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import ot_client  # noqa: E402

GROUP_ROOT = "MONDO_0000141"          # mosaic variegated aneuploidy syndrome (grouped)
EXTRA_DISEASES = ["EFO_0009077"]      # premature chromatid separation trait

DISEASE_Q = """query($id:String!){
  disease(efoId:$id){
    id name dbXRefs
    children{ id name }
    descendants
    drugAndClinicalCandidates{ count rows{ maxClinicalStage
      drug{ id name maximumClinicalStage } } }
    associatedTargets(page:{index:0,size:50}){
      count rows{ score target{ id approvedSymbol biotype
        genomicLocation{ chromosome start end } } } }
  }
}"""


def disease(did: str) -> dict | None:
    return ot_client.query(DISEASE_Q, {"id": did})["disease"]


def overlap_note(symbol: str, loc: dict | None, known: dict[str, dict]) -> str:
    """Flag a target whose locus overlaps a checklist gene, or that is a
    readthrough of one. Positional overlap is the commonest way a gene with no
    role in the disease acquires an association score."""
    if "-" in symbol:
        parts = [p for p in symbol.split("-") if p in known]
        if parts:
            return f"readthrough locus involving {', '.join(parts)}"
    if not loc:
        return ""
    for gene, ref in known.items():
        if gene == symbol or not ref.get("chromosome"):
            continue
        if (loc["chromosome"] == ref["chromosome"]
                and loc["start"] <= ref["end"] and ref["start"] <= loc["end"]):
            return f"locus overlaps {gene} (chr{loc['chromosome']})"
    return ""


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-score", type=float, default=0.3,
                    help="association score floor for the target walk")
    ap.add_argument("--checklist", type=pathlib.Path,
                    default=pathlib.Path("config/mva_checklist_genes.tsv"))
    ap.add_argument("--out-dir", type=pathlib.Path, default=pathlib.Path("."))
    args = ap.parse_args()

    info = ot_client.meta()
    print(f"Open Targets data release {info['data_version']}", file=sys.stderr)

    # --- the disease group -------------------------------------------------
    root = disease(GROUP_ROOT)
    if root is None:
        print(f"{GROUP_ROOT} is not in this release", file=sys.stderr)
        return 1
    ids = [GROUP_ROOT] + [c["id"] for c in root["children"]] + EXTRA_DISEASES
    seen, diseases = set(), []
    for did in ids:
        if did in seen:
            continue
        seen.add(did)
        rec = disease(did)
        if rec is not None:
            diseases.append(rec)

    disease_rows, all_targets = [], {}
    for rec in diseases:
        drugs = {r["drug"]["id"]: r["drug"]
                 for r in (rec["drugAndClinicalCandidates"]["rows"] or []) if r["drug"]}
        approved = sum(1 for d in drugs.values() if d["maximumClinicalStage"] == "APPROVAL")
        xrefs = [x for x in (rec["dbXRefs"] or [])
                 if x.startswith(("OMIM", "Orphanet", "OMIMPS"))]
        disease_rows.append({
            "disease_id": rec["id"], "disease": rec["name"],
            "xrefs": "; ".join(xrefs[:4]),
            "n_drugs": len(drugs), "n_approved": approved,
            "drugs": "; ".join(sorted(d["name"] for d in drugs.values())) or "none",
        })
        for row in rec["associatedTargets"]["rows"]:
            tgt = row["target"]
            prev = all_targets.get(tgt["approvedSymbol"])
            if prev is None or row["score"] > prev["score"]:
                all_targets[tgt["approvedSymbol"]] = {"score": row["score"], "target": tgt}
        print(f"  {rec['id']:16s} {rec['name'][:52]:52s} drugs={len(drugs)} "
              f"approved={approved}", file=sys.stderr)

    # --- the target walk ---------------------------------------------------
    with args.checklist.open(newline="") as fh:
        checklist = {r["gene"]: r for r in csv.DictReader(fh, delimiter="\t")}
    known_loci: dict[str, dict] = {}
    for gene in checklist:
        eid = ot_client.resolve_symbol(gene)
        if not eid:
            continue
        rec = ot_client.query(
            "query($id:String!){target(ensemblId:$id){genomicLocation"
            "{chromosome start end}}}", {"id": eid})["target"]
        if rec and rec["genomicLocation"]:
            known_loci[gene] = rec["genomicLocation"]

    gene_rows = []
    for symbol, entry in sorted(all_targets.items(), key=lambda kv: -kv[1]["score"]):
        if entry["score"] < args.min_score:
            continue
        tgt = entry["target"]
        record = ot_client.target_drugs(tgt["id"])
        drugs = ot_client.unique_drugs(record, record["approvedSymbol"]) if record else {}
        approved = [d for d in drugs.values() if d["max_stage_drug"] == "APPROVAL"]
        gene_rows.append({
            "gene": symbol, "ensembl": tgt["id"],
            "association_score": round(entry["score"], 3),
            "on_mva_checklist": "yes" if symbol in checklist else "no",
            "mva_subtype": checklist.get(symbol, {}).get("mva_subtype", ""),
            "n_drugs": len(drugs), "n_approved": len(approved),
            "approved_drugs": "; ".join(sorted(d["drug"] for d in approved)) or "none",
            "annotation_caveat": overlap_note(symbol, tgt.get("genomicLocation"), known_loci),
        })
        print(f"  {symbol:14s} score={entry['score']:.2f} drugs={len(drugs):3d} "
              f"approved={len(approved)}", file=sys.stderr)

    out_dir = args.out_dir
    for rows, name in ((disease_rows, "data/mva_disease_drugs.tsv"),
                       (gene_rows, "data/mva_group_genes.tsv")):
        path = out_dir / name
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                                    lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        print(f"{len(rows):4d} rows -> {path}", file=sys.stderr)

    summary = {
        "open_targets_release": info["data_version"],
        "disease_entries": len(disease_rows),
        "disease_entries_with_any_drug": sum(1 for r in disease_rows if r["n_drugs"]),
        "approved_for_the_disease": sum(r["n_approved"] for r in disease_rows),
        "genes_above_score": len(gene_rows),
        "genes_on_checklist": sum(1 for r in gene_rows if r["on_mva_checklist"] == "yes"),
        "genes_with_any_drug": sum(1 for r in gene_rows if r["n_drugs"]),
        "approved_across_group_genes": sum(r["n_approved"] for r in gene_rows),
        "annotation_caveats": [r["gene"] for r in gene_rows if r["annotation_caveat"]],
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
