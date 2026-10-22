#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Run the whole Track 2 search, or check that the committed tables still hold.

The pipeline has two modes and the second one is the point of it.

`--mode run` regenerates the tables from the live Open Targets API.

`--mode check` regenerates them into a temporary directory and diffs them
against what is committed, without touching the repository. This turns every
number in the report into something a reader can falsify in one command, and it
is also how a change in an upstream data release gets noticed: Open Targets
versions its releases, so a count taken in one release is not guaranteed in the
next. A silent drift is the failure mode that matters here, because the central
claim of this work is a *negative* one - that no approved medicine acts on any
node of the broken machinery - and a negative claim rots quietly when the
database underneath it moves.

Stages:
    1  census      every node counted, zeros included        (needs network)
    2  candidates  compounds for the nodes that have any     (needs network)
    3  gate        direction of action, per tissue context   (offline)
    4  ranking     arithmetic, ordering and sourcing checks  (offline)
    5  citations   every PMID resolved against PubMed        (needs network)

Usage:
    python3 src/run_pipeline.py --mode check
    python3 src/run_pipeline.py --mode check --offline       # stages 3 and 4 only
    python3 src/run_pipeline.py --mode run
"""
from __future__ import annotations

import argparse
import pathlib
import shutil
import subprocess
import sys
import tempfile

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent

REGENERATED = [
    "data/node_drug_census.tsv",
    "data/cohesion_layer_census.tsv",
    "data/candidates_raw.tsv",
    "data/candidates_gated.tsv",
    "data/rejected_direction.tsv",
    "data/approval_status_openfda.tsv",
]


def run(cmd: list[str], *, cwd: pathlib.Path = ROOT) -> int:
    print(f"\n$ {' '.join(cmd)}", flush=True)
    return subprocess.run([sys.executable, *cmd], cwd=cwd).returncode


# Columns that CANNOT be regenerated from the public endpoint and are carried as
# curated values instead. The public v4 `Drug` type exposes only id, name,
# tradeNames, maximumClinicalStage and mechanismsOfAction - there is no
# hasBeenWithdrawn and no drug_type - so a fresh run leaves these blank or marks
# them "curated". Comparing them would report a difference on every row forever
# and train the reader to ignore the check, which is worse than not having one.
# They are excluded from the comparison and named here so the exclusion is visible.
CURATED_COLUMNS = {
    "data/candidates_raw.tsv": ("drug_type", "withdrawn"),
}


def normalise(path: pathlib.Path, ignore: tuple[str, ...] = ()) -> list[str]:
    """Compare on content, not line endings: the committed tables were written
    with the csv module's default CRLF and the scripts here emit LF.

    `ignore` drops named columns before comparing, for fields that are curated
    rather than regenerated (see CURATED_COLUMNS)."""
    if not path.exists():
        return []
    lines = [line.rstrip("\r\n") for line in path.read_text().splitlines()]
    if not ignore or not lines:
        return lines
    header = lines[0].split("\t")
    keep = [i for i, c in enumerate(header) if c not in ignore]
    return ["\t".join(row.split("\t")[i] for i in keep if i < len(row.split("\t")))
            for row in lines]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["run", "check"], default="check")
    ap.add_argument("--offline", action="store_true",
                    help="skip the stages that need the network")
    ap.add_argument("--keep", action="store_true", help="keep the scratch directory")
    args = ap.parse_args()

    failures: list[str] = []
    scratch = pathlib.Path(tempfile.mkdtemp(prefix="track2-check-"))
    out = ROOT if args.mode == "run" else scratch
    print(f"mode={args.mode}  writing to {out}")

    if not args.offline:
        if run(["src/census_nodes.py", "--out-dir", str(out), "--quiet"]):
            failures.append("census")
        census = out / "data/node_drug_census.tsv"
        if census.exists() and run(["src/fetch_candidates.py", "--census", str(census),
                                    "--out", str(out / "data/candidates_raw.tsv")]):
            failures.append("candidates")

    raw = out / "data/candidates_raw.tsv"
    if not raw.exists():
        raw = ROOT / "data/candidates_raw.tsv"
        print(f"\n(using the committed {raw.relative_to(ROOT)} for the gate)")
    if run(["src/gate_direction.py", "--candidates", str(raw),
            "--out-dir", str(out), "--self-test"]):
        failures.append("direction gate self-test")

    if run(["src/validate_ranking.py", "--strict"]):
        failures.append("ranking validation")

    # Product-level approval. Kept behind --offline because it needs openFDA, and
    # separated from the citation check because its failure mode is different: an
    # unreachable endpoint must surface as UNRESOLVED rather than be mistaken for
    # an absent product. --skip-chembl because the ChEMBL mechanism endpoint is
    # intermittently unavailable and would otherwise stall the whole run.
    if not args.offline:
        (out / "data").mkdir(parents=True, exist_ok=True)
        if run(["src/verify_approval.py", "--skip-chembl",
                "--ranked", str(ROOT / "results/ranked_candidates.tsv"),
                "--out", str(out / "data/approval_status_openfda.tsv")]):
            failures.append("product-level approval check")

    if not args.offline and run(["src/verify_citations.py", "--fail-on-missing"]):
        failures.append("citation check")

    if args.mode == "check":
        print("\n--- regenerated vs committed ---")
        for rel in REGENERATED:
            ignore = CURATED_COLUMNS.get(rel, ())
            fresh = normalise(out / rel, ignore)
            committed = normalise(ROOT / rel, ignore)
            if not fresh:
                print(f"  {rel:38s} not regenerated in this run")
                continue
            if fresh == committed:
                note = f", {len(ignore)} curated cols excluded" if ignore else ""
                print(f"  {rel:38s} identical ({len(fresh) - 1} rows{note})")
            else:
                print(f"  {rel:38s} DIFFERS "
                      f"({len(committed) - 1} committed vs {len(fresh) - 1} regenerated)")
                failures.append(f"{rel} differs from the committed table")

    if not args.keep and args.mode == "check":
        shutil.rmtree(scratch, ignore_errors=True)
    elif args.keep:
        print(f"\nscratch kept at {scratch}")

    print()
    if failures:
        for f in failures:
            print(f"FAILED: {f}", file=sys.stderr)
        return 1
    print("all stages passed", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
