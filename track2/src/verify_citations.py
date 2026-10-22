#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Check that every PMID cited in this repository exists and says what is claimed.

This enforces, in code, the rule the whole track is built on: no drug name and
no mechanism enters an output from memory. A language model will produce a
plausible PMID for a claim it has invented, and a plausible PMID is
indistinguishable from a real one until it is looked up. So every identifier in
the repository is extracted and resolved against PubMed, and the title and year
are printed next to the file that cites it, for a human to read.

What this does and does not establish: resolving a PMID proves the article
exists and tells you its title. It does not prove the article supports the
sentence citing it. The output is a worklist for that check, not a substitute.

Usage:
    python3 src/verify_citations.py
    python3 src/verify_citations.py --fail-on-missing     # for CI
"""
from __future__ import annotations

import argparse
import collections
import csv
import io
import json
import pathlib
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
PMID_RE = re.compile(r"PMID[:\s]*(\d{7,8})")
BARE_PMID_RE = re.compile(r"\b(\d{7,8})\b")
CHEMBL_RE = re.compile(r"\b(CHEMBL\d+)\b")
SEARCH_EXT = ("*.md", "data/*.tsv", "results/*.tsv", "config/*.tsv")


def collect(root: pathlib.Path) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    pmids: dict[str, set[str]] = collections.defaultdict(set)
    chembl: dict[str, set[str]] = collections.defaultdict(set)
    for pattern in SEARCH_EXT:
        for path in sorted(root.glob(pattern)):
            text = path.read_text(errors="replace")
            for pmid in PMID_RE.findall(text):
                pmids[pmid].add(str(path.relative_to(root)))
            for cid in CHEMBL_RE.findall(text):
                chembl[cid].add(str(path.relative_to(root)))
            # Tables cite bare identifiers in a pmid column, with no "PMID" prefix
            # for the pattern above to find. Without this branch the claim that
            # every identifier in the repository is re-resolved is not true.
            if path.suffix == ".tsv":
                reader = csv.DictReader(io.StringIO(text), delimiter="\t")
                cols = [c for c in (reader.fieldnames or []) if "pmid" in c.lower()]
                if cols:
                    for row in reader:
                        for col in cols:
                            for pmid in BARE_PMID_RE.findall(row.get(col) or ""):
                                pmids[pmid].add(str(path.relative_to(root)))
    return pmids, chembl


def fetch_summaries(pmids: list[str], email: str | None = None) -> dict[str, dict]:
    """Resolve PMIDs in batches. Unresolved ids are returned as None entries."""
    out: dict[str, dict] = {}
    for i in range(0, len(pmids), 150):
        batch = pmids[i:i + 150]
        params = {"db": "pubmed", "id": ",".join(batch), "retmode": "xml"}
        if email:
            params["email"] = email
        url = f"{ESUMMARY}?{urllib.parse.urlencode(params)}"
        with urllib.request.urlopen(url, timeout=90) as fh:
            root = ET.parse(fh).getroot()
        for doc in root.findall(".//DocSum"):
            pmid = doc.findtext("Id")
            fields = {item.get("Name"): (item.text or "")
                      for item in doc.findall("Item")}
            out[pmid] = {"title": fields.get("Title", ""),
                         "source": fields.get("Source", ""),
                         "pubdate": fields.get("PubDate", "")}
        for pmid in batch:
            out.setdefault(pmid, None)
        if i + 150 < len(pmids):
            time.sleep(0.4)          # NCBI asks for <=3 requests/second
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=pathlib.Path, default=pathlib.Path("."))
    ap.add_argument("--email", default=None,
                    help="contact address for NCBI; optional, omit if you have none")
    ap.add_argument("--fail-on-missing", action="store_true")
    ap.add_argument("--json", type=pathlib.Path, default=None)
    args = ap.parse_args()

    pmids, chembl = collect(args.root)
    print(f"{len(pmids)} distinct PMIDs and {len(chembl)} ChEMBL ids cited\n", file=sys.stderr)
    if not pmids:
        return 0

    summaries = fetch_summaries(sorted(pmids), args.email)
    missing = [p for p, s in summaries.items() if s is None]

    for pmid in sorted(pmids, key=int):
        info = summaries.get(pmid)
        where = ", ".join(sorted(pmids[pmid]))
        if info is None:
            print(f"UNRESOLVED  PMID:{pmid}  cited in {where}")
            continue
        title = info["title"][:96]
        print(f"ok  PMID:{pmid}  {info['source']} {info['pubdate']}\n"
              f"    {title}\n    cited in {where}")

    print(f"\nresolved {len(pmids) - len(missing)}/{len(pmids)}", file=sys.stderr)
    if missing:
        print(f"UNRESOLVED: {', '.join(missing)}", file=sys.stderr)

    if args.json:
        args.json.write_text(json.dumps(
            {"pmids": {p: {"summary": summaries.get(p), "cited_in": sorted(v)}
                       for p, v in pmids.items()},
             "chembl_ids": {c: sorted(v) for c, v in chembl.items()}}, indent=2))

    return 1 if (missing and args.fail_on_missing) else 0


if __name__ == "__main__":
    raise SystemExit(main())
