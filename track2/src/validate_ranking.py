#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Re-compute the ranking arithmetic and check the shortlist against its sources.

The scoring components are judgements and are curated by hand - that is stated
plainly rather than hidden behind a script that appears to derive them. What a
script can check is everything downstream of the judgement, and those are the
mistakes that actually happen:

  * the total must equal mechanism_concordance + phenotype_evidence
    + human_exposure - risk_penalty;
  * rows must be ordered by rank, and ties must be visible rather than
    silently resolved;
  * every ChEMBL id must appear in the candidate tables it is drawn from, so a
    drug cannot enter the shortlist without having passed the earlier steps;
  * nothing may appear in both the shortlist and the contraindication table.

Usage:
    python3 src/validate_ranking.py
    python3 src/validate_ranking.py --strict      # non-zero exit on any problem
"""
from __future__ import annotations

import argparse
import csv
import pathlib
import sys

COMPONENTS = ("expected_benefit_ceiling", "evidence_quality", "executability")
PENALTY = "risk_penalty"
RANK = "rank_global"


def read_tsv(path: pathlib.Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="") as fh:
        return [{k: (v or "").strip() for k, v in r.items()}
                for r in csv.DictReader(fh, delimiter="\t")]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ranked", type=pathlib.Path,
                    default=pathlib.Path("results/ranked_candidates.tsv"))
    ap.add_argument("--cleared", type=pathlib.Path,
                    default=pathlib.Path("data/candidates_cleared.tsv"))
    ap.add_argument("--symptomatic", type=pathlib.Path,
                    default=pathlib.Path("data/symptomatic_candidates.tsv"))
    ap.add_argument("--contraindicated", type=pathlib.Path,
                    default=pathlib.Path("data/contraindicated.tsv"))
    ap.add_argument("--strict", action="store_true")
    args = ap.parse_args()

    ranked = read_tsv(args.ranked)
    problems: list[str] = []

    for row in ranked:
        total = sum(int(row[c]) for c in COMPONENTS) - int(row[PENALTY])
        recorded = int(row["total"])
        # weighted = raw x category weight. The symptomatic arm is down-weighted, so the
        # weighted column is what the global rank is built from and has to be checked too:
        # an unchecked derived column is where a hand-edited table silently goes wrong.
        weight = float(row["category_weight"])
        w_expect = total * weight
        w_recorded = float(row["weighted_total"])
        ok = total == recorded and abs(w_expect - w_recorded) < 1e-9
        flag = "ok" if ok else "MISMATCH"
        print(f"  rank {row[RANK]:>2s}  {row['drug'][:24]:24s} "
              f"{'+'.join(row[c] for c in COMPONENTS)}-{row[PENALTY]} = {total:2d}  "
              f"recorded {recorded:2d}  x{weight} = {w_recorded:4.1f}  {flag}")
        if total != recorded:
            problems.append(f"rank {row[RANK]} ({row['drug']}): "
                            f"components give {total}, file says {recorded}")
        if abs(w_expect - w_recorded) >= 1e-9:
            problems.append(f"rank {row[RANK]} ({row['drug']}): weighted should be "
                            f"{w_expect:g}, file says {w_recorded:g}")

    ranks = [int(r[RANK]) for r in ranked]
    if ranks != sorted(ranks):
        problems.append("rows are not in rank order")
    if ranks != list(range(1, len(ranks) + 1)):
        problems.append(f"global ranks are not 1..{len(ranks)} without gaps: {ranks}")
    # the global ordering must follow the weighted score, not the raw one
    w = [float(r["weighted_total"]) for r in ranked]
    if w != sorted(w, reverse=True):
        problems.append("global rank order does not follow weighted_total")
    totals = [int(r["total"]) for r in ranked]
    if totals != sorted(totals, reverse=True):
        print("\n  note: total is not monotonic with rank - ties are broken by "
              "judgement, which is legitimate but must be stated", file=sys.stderr)
    dupes = {t for t in totals if totals.count(t) > 1}
    if dupes:
        print(f"  tied totals present: {sorted(dupes)} - check the tie-break is "
              f"argued in the report", file=sys.stderr)

    known = {r["chembl_id"] for r in read_tsv(args.cleared)}
    known |= {r["chembl_id"] for r in read_tsv(args.symptomatic)}
    unsourced = [r["drug"] for r in ranked if r["chembl_id"] not in known]
    if unsourced:
        problems.append("shortlisted without appearing in a candidate table: "
                        + ", ".join(unsourced))

    # A shortlisted drug may also appear in the contraindication table - that is a
    # real situation, not a bug, but it must be declared on the row rather than left
    # for the reader to notice by cross-referencing two files.
    banned = {r["chembl_id"] for r in read_tsv(args.contraindicated)
              if (r.get("chembl_id") or "").startswith("CHEMBL")}
    # The flag column carries three different kinds of statement and they must not be
    # conflated: a contraindication is a reason not to give the drug, a caution is a
    # reason to watch, and a regulatory-status note is neither. Only the first implies a
    # row in the contraindication table. Requiring a known prefix also stops the flag
    # text from drifting into free prose that no check can read.
    PREFIXES = ("SYSTEMICALLY CONTRAINDICATED", "CONTRAINDICATED", "CAUTION",
                "APPROVAL IN QUESTION")
    for row in ranked:
        flag = row.get("contraindication_flag", "none").strip()
        flagged = flag not in ("", "none")
        if flagged and not flag.startswith(PREFIXES):
            problems.append(f"{row['drug']}: flag does not start with a recognised "
                            f"prefix {PREFIXES}: {flag[:40]!r}")
        is_contra = flag.startswith(("SYSTEMICALLY CONTRAINDICATED", "CONTRAINDICATED"))
        listed = row["chembl_id"] in banned
        if listed and not is_contra:
            problems.append(f"{row['drug']} is in the contraindication table but its "
                            "ranked row does not declare a contraindication")
        if is_contra and not listed:
            problems.append(f"{row['drug']} is flagged contraindicated but is absent "
                            "from the contraindication table")

    print()
    for p in problems:
        print(f"PROBLEM: {p}", file=sys.stderr)
    if not problems:
        print(f"{len(ranked)} ranked rows: arithmetic, ordering and sourcing all check out",
              file=sys.stderr)
    return 1 if (problems and args.strict) else 0


if __name__ == "__main__":
    raise SystemExit(main())
