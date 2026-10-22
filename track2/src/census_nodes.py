#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Count the drugs available against every node of the broken machinery.

This is the central result of Track 2 and it is a negative one, so it is
written as a census rather than a search: every node on the list is queried,
and a node with no drugs produces a row saying zero rather than no row at all.
A search that returns nothing and a census that returns zeros look the same in
a summary and are not the same evidence.

Two node lists are counted separately and deliberately:

  config/sac_nodes.tsv       spindle assembly checkpoint and kinetochore
  config/cohesion_nodes.tsv  sister-chromatid cohesion

The cohesion layer is separate because premature chromatid separation is the
cytogenetic hallmark of this disease, so "nothing available" has to be
established for that layer on its own rather than inherited from the first.

Usage:
    python3 src/census_nodes.py                        # both lists
    python3 src/census_nodes.py --list sac             # one list
    python3 src/census_nodes.py --out-dir /tmp/check   # write elsewhere
"""
from __future__ import annotations

import argparse
import csv
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import ot_client  # noqa: E402

LISTS = {
    "sac": ("config/sac_nodes.tsv", "data/node_drug_census.tsv", "sac_kinetochore"),
    "cohesion": ("config/cohesion_nodes.tsv", "data/cohesion_layer_census.tsv", "cohesion"),
}


def read_nodes(path: pathlib.Path) -> list[dict]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def census(nodes: list[dict], *, verbose: bool = True) -> tuple[list[dict], list[str]]:
    """One row per node: compounds recorded, and how many are approved.

    `n_drug_rows` counts unique ChEMBL ids, not connection rows - a drug listed
    against several indications is one compound. Unresolved symbols are
    returned separately instead of being written out as zeros, because a
    look-up failure and a genuine absence of drugs are different facts.
    """
    rows, unresolved = [], []
    for i, node in enumerate(nodes, 1):
        symbol = node["node"]
        eid = node.get("ensembl") or ot_client.resolve_symbol(symbol)
        if not eid:
            unresolved.append(symbol)
            if verbose:
                print(f"  [{i:3d}/{len(nodes)}] {symbol:12s} UNRESOLVED", file=sys.stderr)
            continue
        record = ot_client.target_drugs(eid)
        if record is None:
            unresolved.append(symbol)
            if verbose:
                print(f"  [{i:3d}/{len(nodes)}] {symbol:12s} NOT IN RELEASE", file=sys.stderr)
            continue
        drugs = ot_client.unique_drugs(record, record["approvedSymbol"])
        approved = sum(1 for d in drugs.values() if d["max_stage_drug"] == "APPROVAL")
        row = {"node": symbol, "ensembl": eid,
               "n_drug_rows": len(drugs), "n_approved": approved}
        if "role" in node:
            # The cohesion census is keyed `target` in the committed table; keep the
            # header stable so a regenerated file diffs cleanly against it.
            row = {"target": symbol, "role": node["role"],
                   "n_drug_rows": len(drugs), "n_approved": approved}
        rows.append(row)
        if verbose:
            print(f"  [{i:3d}/{len(nodes)}] {symbol:12s} drugs={len(drugs):3d} approved={approved}",
                  file=sys.stderr)
    return rows, unresolved


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", choices=[*LISTS, "both"], default="both")
    ap.add_argument("--out-dir", type=pathlib.Path, default=pathlib.Path("."))
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    info = ot_client.meta()
    print(f"Open Targets data release {info['data_version']} (API {info['api_version']})",
          file=sys.stderr)

    chosen = list(LISTS) if args.list == "both" else [args.list]
    summary = {"open_targets_release": info["data_version"], "lists": {}}

    for name in chosen:
        src, dst, _ = LISTS[name]
        nodes = read_nodes(pathlib.Path(src))
        print(f"\n{name}: {len(nodes)} nodes from {src}", file=sys.stderr)
        rows, unresolved = census(nodes, verbose=not args.quiet)

        out = args.out_dir / dst
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                                    lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

        no_drug = sum(1 for r in rows if r["n_drug_rows"] == 0)
        compounds = sum(r["n_drug_rows"] for r in rows)
        approved = sum(r["n_approved"] for r in rows)
        summary["lists"][name] = {
            "nodes": len(rows), "nodes_without_any_drug": no_drug,
            "compounds_total": compounds, "approved_total": approved,
            "unresolved": unresolved,
        }
        print(f"  -> {out}", file=sys.stderr)
        print(f"  {len(rows)} nodes; {no_drug} with no drug of any kind; "
              f"{compounds} compounds; {approved} approved", file=sys.stderr)

    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
