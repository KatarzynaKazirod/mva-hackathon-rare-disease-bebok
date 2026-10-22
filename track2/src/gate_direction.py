#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Gate candidates on whether the drug pushes the node the way the defect needs.

Direction of action is a gate, not a scoring dimension. A compound that acts on
the right protein in the wrong direction is not a weak candidate to be ranked
below the others; it is the wrong intervention and is removed, with the reason
recorded so the removal can be argued with.

The gate runs each compound against *both* tissue contexts, because the
required direction differs between them: somatic tissue needs segregation
fidelity restored, while the tumour needs pushing past the threshold it can
tolerate. The same drug can therefore be contraindicated systemically and
rational against the sarcoma, and the two verdicts are kept apart.

Inputs
    data/candidates_raw.tsv       compounds, one row per (target, action type)
    data/signed_network.tsv       per node: required direction in each context
    config/action_type_effects.tsv  action type -> direction on target function

Outputs
    data/candidates_gated.tsv     passed, with the context it passed in
    data/rejected_direction.tsv   rejected, with the reason

The acceptance test is built in: MPS1/TTK inhibitors must fail in somatic
tissue and pass for the tumour. A gate that passes everything, or rejects
everything, is not a gate, and `--self-test` fails loudly if that happens.

Usage:
    python3 src/gate_direction.py
    python3 src/gate_direction.py --self-test
"""
from __future__ import annotations

import argparse
import csv
import pathlib
import sys

CONTEXTS = {"somatic": "required_normal", "tumour": "required_tumour"}


def read_tsv(path: pathlib.Path) -> list[dict]:
    with path.open(newline="") as fh:
        return [{k: (v or "").strip() for k, v in r.items()}
                for r in csv.DictReader(fh, delimiter="\t")]


def satisfies(effect: str, required: str) -> bool:
    """Does `effect` meet `required`?

    `increase_or_neutral` admits an increase or no change but not a decrease;
    `neutral` admits only no change. An unknown effect never satisfies a
    requirement - absence of a direction is not evidence of the right one.
    """
    if required in ("", "not_applicable"):
        return False
    if effect == "unknown":
        return False
    if required == "increase_or_neutral":
        return effect in ("increase", "neutral")
    return effect == required


def gate(candidates: list[dict], network: list[dict],
         effects: dict[str, str]) -> tuple[list[dict], list[dict]]:
    required = {n["node"]: n for n in network}
    passed: list[dict] = []
    rejected: list[dict] = []

    for cand in candidates:
        target = cand["target"]
        node = required.get(target)
        effect = effects.get(cand["action_type"], "unknown")

        if node is None:
            rejected.append({**cand, "drug_effect": effect, "context": "n/a",
                             "required": "n/a",
                             "reject_reason": f"{target} is not signed in the network; "
                                              "direction cannot be decided"})
            continue

        for context, column in CONTEXTS.items():
            req = node.get(column, "")
            base = {"target": target, "drug_id": cand["drug_id"], "drug": cand["drug"],
                    "action_type": cand["action_type"], "drug_effect": effect,
                    "context": context, "required": req,
                    "max_stage_drug": cand["max_stage_drug"]}
            if req == "not_applicable":
                continue
            if satisfies(effect, req):
                passed.append(base)
            else:
                rejected.append({**base, "reject_reason":
                                 f"{target} requires {req} in {context} tissue; "
                                 f"{cand['action_type']} produces {effect}"})
    return passed, rejected


def self_test(passed: list[dict], rejected: list[dict]) -> list[str]:
    """Cases whose answer is known in advance. Failing any of these voids the gate."""
    problems = []
    if not passed:
        problems.append("nothing passed - a gate that rejects everything is not a gate")
    if not rejected:
        problems.append("nothing was rejected - the gate is not filtering")

    ttk_som = [r for r in rejected if r["target"] == "TTK" and r.get("context") == "somatic"]
    ttk_tum = [r for r in passed if r["target"] == "TTK" and r["context"] == "tumour"]
    if not ttk_som:
        problems.append("TTK inhibitors were not rejected for somatic tissue")
    if not ttk_tum:
        problems.append("TTK inhibitors did not pass for the tumour")

    approved = [r for r in passed if r["max_stage_drug"] == "APPROVAL"]
    if approved:
        problems.append(f"{len(approved)} approved drugs passed - the reported result "
                        "is that none does; re-check before publishing")
    return problems


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidates", type=pathlib.Path, default=pathlib.Path("data/candidates_raw.tsv"))
    ap.add_argument("--network", type=pathlib.Path, default=pathlib.Path("data/signed_network.tsv"))
    ap.add_argument("--effects", type=pathlib.Path,
                    default=pathlib.Path("config/action_type_effects.tsv"))
    ap.add_argument("--out-dir", type=pathlib.Path, default=pathlib.Path("."))
    ap.add_argument("--self-test", action="store_true",
                    help="exit non-zero if a known-answer case comes out wrong")
    args = ap.parse_args()

    effects = {r["action_type"]: r["effect_on_target_function"]
               for r in read_tsv(args.effects)}
    passed, rejected = gate(read_tsv(args.candidates), read_tsv(args.network), effects)

    pass_cols = ["target", "drug_id", "drug", "action_type", "drug_effect",
                 "context", "required", "max_stage_drug"]
    for rows, cols, name in ((passed, pass_cols, "data/candidates_gated.tsv"),
                             (rejected, pass_cols + ["reject_reason"],
                              "data/rejected_direction.tsv")):
        out = args.out_dir / name
        out.parent.mkdir(parents=True, exist_ok=True)
        # Sort before writing. Row order here is otherwise inherited from the
        # input and from dict iteration, which made these two tables differ
        # between runs that agreed on every value - a reproducibility check
        # that fails on row order is a check nobody will trust for long.
        rows = sorted(rows, key=lambda r: (r.get("context", ""), r.get("target", ""),
                                           r.get("drug_id", "")))
        with out.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                                    lineterminator="\n", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        print(f"{len(rows):4d} rows -> {out}", file=sys.stderr)

    by_context = {}
    for row in passed:
        by_context[row["context"]] = by_context.get(row["context"], 0) + 1
    print(f"passed by context: {by_context}", file=sys.stderr)
    print(f"approved among passed: "
          f"{sum(1 for r in passed if r['max_stage_drug'] == 'APPROVAL')}", file=sys.stderr)

    problems = self_test(passed, rejected)
    for p in problems:
        print(f"SELF-TEST: {p}", file=sys.stderr)
    if args.self_test and problems:
        return 1
    if args.self_test:
        print("SELF-TEST: all known-answer cases behaved as expected", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
