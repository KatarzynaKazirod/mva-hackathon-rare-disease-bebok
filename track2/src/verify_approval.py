#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Verify approval at PRODUCT level, and mechanism at SOURCE level.

Two checks the target database cannot settle on its own.

1. APPROVAL. A target database reports approval for a *substance*,
   jurisdiction-agnostically. That is not evidence an approved medicinal
   product exists: one NAD+ precursor returns approved status while being
   marketed as a supplement. openFDA /drug/ndc returns the actual US
   registration category per product (NDA / ANDA / BLA / bulk ingredient),
   which is what "already in human use" has to mean here.

2. MECHANISM. Open Targets does not fully mirror ChEMBL's mechanism table.
   Demonstrated: sodium phenylbutyrate returns no mechanism row in Open
   Targets but one in ChEMBL, and acetylcysteine returns none and two. Any
   count of "the proposed mechanism is not the registered one" built on the
   Open Targets field alone therefore risks reading an ABSENT row as a
   DIFFERENT mechanism. This script records both sides so the distinction is
   auditable rather than asserted.

Neither check decides whether a candidate is good. They decide what is known
about it, which is a different question and the only one a database can answer.

Usage:
    python3 src/verify_approval.py --ranked results/ranked_candidates.tsv \\
                                   --out data/approval_status_openfda.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import pathlib
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

FDA_NDC = "https://api.fda.gov/drug/ndc.json"
CHEMBL_MECH = "https://www.ebi.ac.uk/chembl/api/data/mechanism"


def _get(url: str, tries: int = 3, timeout: int = 25, wait: float = 3.0):
    """Fetch JSON. Returns {"_err": ...} rather than raising: an endpoint being
    down is a result to record, not a reason to lose the rest of the run."""
    last = ""
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=timeout) as fh:
                return json.load(fh)
        except urllib.error.HTTPError as exc:
            # openFDA answers 404 for "no product matches", which is a real
            # finding (no approved product) and not a transport failure.
            if exc.code == 404:
                return {"_notfound": True}
            last = f"HTTP {exc.code}"
        except Exception as exc:  # noqa: BLE001 - network shapes vary
            last = str(exc)[:60]
        if attempt < tries - 1:
            time.sleep(wait)
    return {"_err": last}


def openfda_products(drug: str) -> dict:
    """Registration categories for every marketed product of this substance."""
    q = urllib.parse.quote(f'active_ingredients.name:"{drug}"')
    d = _get(f"{FDA_NDC}?search={q}&limit=100")
    if "_notfound" in d:
        return {"n": 0, "categories": {}, "note": "no NDC record"}
    if "_err" in d:
        return {"n": None, "categories": {}, "note": f"openFDA unreachable ({d['_err']})"}
    cats: dict[str, int] = {}
    earliest = None
    for r in d.get("results", []):
        cat = r.get("marketing_category", "UNKNOWN")
        cats[cat] = cats.get(cat, 0) + 1
        start = r.get("marketing_start_date")
        if start and (earliest is None or start < earliest):
            earliest = start
    return {"n": len(d.get("results", [])), "categories": cats, "earliest": earliest, "note": ""}


def chembl_mechanisms(chembl_id: str) -> dict:
    d = _get(f"{CHEMBL_MECH}?format=json&molecule_chembl_id={chembl_id}&limit=25")
    if "_err" in d:
        return {"n": None, "rows": [], "note": f"ChEMBL unreachable ({d['_err']})"}
    rows = [f"{m.get('action_type')}: {m.get('mechanism_of_action')}"
            for m in d.get("mechanisms", [])]
    return {"n": len(rows), "rows": rows, "note": ""}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ranked", type=pathlib.Path,
                    default=pathlib.Path("results/ranked_candidates.tsv"))
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("data/approval_status_openfda.tsv"))
    ap.add_argument("--skip-chembl", action="store_true",
                    help="openFDA only; the ChEMBL mechanism endpoint is "
                         "intermittently unavailable and blocks long runs")
    args = ap.parse_args()

    with args.ranked.open(newline="") as fh:
        cand = [{k: (v or "").strip() for k, v in r.items()}
                for r in csv.DictReader(fh, delimiter="\t")]

    out = []
    unresolved = 0
    for row in cand:
        drug, cid = row["drug"], row["chembl_id"]
        fda = openfda_products(drug)
        mech = {"n": None, "rows": [], "note": "not queried"} if args.skip_chembl \
            else chembl_mechanisms(cid)

        cats = fda["categories"]
        real = sum(v for k, v in cats.items()
                   if "NDA" in k.upper() or "ANDA" in k.upper() or "BLA" in k.upper())
        if fda["n"] is None:
            verdict = "UNRESOLVED - openFDA unreachable"
            unresolved += 1
        elif real:
            verdict = "YES - approved US product"
        elif fda["n"]:
            verdict = "NO - records exist but none is an approved product"
        else:
            verdict = "NO - no US product"

        out.append({
            "drug": drug,
            "chembl_id": cid,
            "openfda_records": "" if fda["n"] is None else str(fda["n"]),
            "marketing_categories": "; ".join(f"{k} x{v}" for k, v in sorted(cats.items())),
            "earliest_marketing_start": fda.get("earliest") or "",
            "product_level_verdict": verdict,
            "open_targets_stage": row.get("approval_source", "")[:60],
            "chembl_mechanism_count": "" if mech["n"] is None else str(mech["n"]),
            "chembl_mechanisms": "; ".join(mech["rows"]),
            "note": "; ".join(x for x in (fda["note"], mech["note"]) if x),
        })
        print(f"  {drug[:26]:26s} openFDA={fda['n']} approved-product={bool(real)} "
              f"chembl_moa={mech['n']}", file=sys.stderr)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)

    yes = sum(1 for r in out if r["product_level_verdict"].startswith("YES"))
    print(f"\n{len(out)} candidates: {yes} with an approved US product, "
          f"{unresolved} unresolved -> {args.out}", file=sys.stderr)
    # An unreachable endpoint must not be silently recorded as an absent product.
    return 1 if unresolved else 0


if __name__ == "__main__":
    raise SystemExit(main())
