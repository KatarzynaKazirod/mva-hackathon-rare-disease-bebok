#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Pull every compound recorded against the nodes that have any.

Reads the census, keeps the nodes with a non-zero count and fetches one row per
compound: ChEMBL id, action type, mechanism of action text and highest clinical
stage reached. This is the raw material the direction gate then filters, and it
is kept as its own file so that the gate can be re-run and argued with without
re-querying the API.

Two columns of `data/candidates_raw.tsv` cannot be produced from the public
endpoint - `withdrawn` and the black-box warning - because the public `Drug`
type does not expose them. They are written as `curated` and must be filled
from the drug label. The pipeline states this rather than emitting a confident
blank.

Usage:
    python3 src/fetch_candidates.py
    python3 src/fetch_candidates.py --census data/node_drug_census.tsv --out /tmp/raw.tsv
"""
from __future__ import annotations

import argparse
import csv
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import ot_client  # noqa: E402

COLUMNS = ["target", "drug_id", "drug", "drug_type", "action_type",
           "max_stage_drug", "moa", "withdrawn"]


def fetch(census_path: pathlib.Path, *, verbose: bool = True) -> list[dict]:
    with census_path.open(newline="") as fh:
        nodes = [r for r in csv.DictReader(fh, delimiter="\t")
                 if int(r["n_drug_rows"]) > 0]

    rows: list[dict] = []
    for node in nodes:
        symbol = node["node"]
        eid = node.get("ensembl") or ot_client.resolve_symbol(symbol)
        record = ot_client.target_drugs(eid)
        drugs = ot_client.unique_drugs(record, record["approvedSymbol"])
        for drug in sorted(drugs.values(), key=lambda d: (d["max_stage_drug"] or "", d["drug"] or "")):
            # A drug with two action types on the same target is reported once per
            # type: the gate decides on direction, and two directions are two
            # different decisions that must not be averaged into one row.
            for action in sorted(drug["action_types"]) or [""]:
                rows.append({
                    "target": symbol,
                    "drug_id": drug["drug_id"],
                    "drug": drug["drug"],
                    "drug_type": "",          # not exposed by the public endpoint
                    "action_type": action,
                    "max_stage_drug": drug["max_stage_drug"],
                    "moa": "; ".join(sorted(drug["moa"])),
                    "withdrawn": "curated",   # see module docstring
                })
        if verbose:
            print(f"  {symbol:10s} {len(drugs):3d} compounds", file=sys.stderr)
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--census", type=pathlib.Path, default=pathlib.Path("data/node_drug_census.tsv"))
    ap.add_argument("--out", type=pathlib.Path, default=pathlib.Path("data/candidates_raw.tsv"))
    args = ap.parse_args()

    rows = fetch(args.census)
    # Deterministic order. The API returns compounds per target in no guaranteed
    # order, which made this table differ between runs that agreed on every
    # value - and a reproducibility check that trips on row order is one a
    # reader soon learns to ignore.
    rows.sort(key=lambda r: (r["target"], r["drug_id"]))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    approved = sum(1 for r in rows if r["max_stage_drug"] == "APPROVAL")
    print(f"{len(rows)} compound rows -> {args.out}", file=sys.stderr)
    print(f"approved among them: {approved}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
