#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Test whether public expression data can speak to the mechanism at all.

Two datasets were examined and the interesting part of the answer is that one
of them cannot answer the question, which is only visible once the data are in
front of you:

  GSE134781  the Sieben allelic series - mouse skeletal muscle, three
             genotypes x three replicates. Skeletal muscle is post-mitotic, so
             a spindle assembly checkpoint is barely expressed there. This
             script computes CPM for a proliferation marker and for the node
             set and shows how many nodes clear a 1-CPM floor. A fold change
             computed on counts in the single digits is noise with a p-value
             attached, and this is how that was caught: a roughly two-fold
             Bub1b reduction was withdrawn after reading the raw counts, and
             the equal drop in the missense genotype excludes the degradation
             mechanism it was being used to support.

  GSE86928   A549 treated with sodium phenylbutyrate. The checkpoint genes do
             fall, but so do the proliferation markers, so the script reports
             both side by side: a drop that comes with MKI67, CCNA2, FOXM1 and
             TOP2A falling is a cell-cycle arrest signature and is not specific
             to the checkpoint. Treating it as target engagement would be
             reading the wrong thing.

Standard library only. Set --no-download to work from files already present.

Usage:
    python3 src/expression_evidence.py --dataset GSE134781
    python3 src/expression_evidence.py --dataset GSE134781 --cache /tmp/geo
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import pathlib
import statistics
import sys
import urllib.request

GEO_SUPP = ("https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}nnn/{acc}/suppl/")

# Proliferation markers. If these are flat or absent the tissue is not cycling
# and nothing about a mitotic checkpoint can be concluded from it.
PROLIFERATION = ["Mki67", "Ccna2", "Ccnb1", "Top2a", "Foxm1", "Plk1", "Bub1b",
                 "MKI67", "CCNA2", "CCNB1", "TOP2A", "FOXM1", "PLK1", "BUB1B"]


def geo_supp_url(accession: str) -> str:
    stub = accession[:-3]
    return GEO_SUPP.format(stub=stub, acc=accession)


def download(url: str, dest: pathlib.Path) -> pathlib.Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"  fetching {url}", file=sys.stderr)
    with urllib.request.urlopen(url, timeout=300) as fh, dest.open("wb") as out:
        while chunk := fh.read(1 << 20):
            out.write(chunk)
    return dest


def read_counts(path: pathlib.Path) -> tuple[list[str], dict[str, list[float]]]:
    """Read a gene x sample counts matrix, gzipped or not."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        delim = "\t" if "\t" in sample.splitlines()[0] else ","
        reader = csv.reader(fh, delimiter=delim)
        header = next(reader)
        samples = header[1:]
        counts: dict[str, list[float]] = {}
        for row in reader:
            if not row or not row[0]:
                continue
            try:
                counts[row[0]] = [float(x or 0) for x in row[1:]]
            except ValueError:
                continue
    return samples, counts


def cpm(counts: dict[str, list[float]]) -> tuple[dict[str, list[float]], list[float]]:
    """Counts per million. Library size is the only normalisation applied here;
    this is a diagnostic of whether the tissue is cycling, not a DE analysis."""
    n = len(next(iter(counts.values())))
    totals = [sum(v[i] for v in counts.values()) for i in range(n)]
    out = {g: [1e6 * v[i] / totals[i] if totals[i] else 0.0 for i in range(n)]
           for g, v in counts.items()}
    return out, totals


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", default="GSE134781")
    ap.add_argument("--counts", type=pathlib.Path, default=None,
                    help="local counts matrix; skips the download")
    ap.add_argument("--cache", type=pathlib.Path, default=pathlib.Path("geo_cache"))
    ap.add_argument("--nodes", type=pathlib.Path, default=pathlib.Path("config/sac_nodes.tsv"))
    ap.add_argument("--cpm-floor", type=float, default=1.0)
    args = ap.parse_args()

    path = args.counts
    if path is None:
        print(f"Supplementary files for {args.dataset} are listed at:\n"
              f"  {geo_supp_url(args.dataset)}\n"
              f"Download the counts matrix and pass it with --counts.\n"
              f"Filenames differ per series, so this is left explicit rather than guessed.",
              file=sys.stderr)
        return 2

    samples, counts = read_counts(path)
    norm, totals = cpm(counts)
    print(f"{len(counts)} genes x {len(samples)} samples", file=sys.stderr)
    print(f"library sizes (millions): "
          f"{[round(t/1e6, 1) for t in totals]}", file=sys.stderr)

    with args.nodes.open(newline="") as fh:
        wanted = {r["node"] for r in csv.DictReader(fh, delimiter="\t")}
    index = {g.upper(): g for g in counts}
    mapped = {w: index[w.upper()] for w in wanted if w.upper() in index}

    print(f"\nnode set: {len(mapped)}/{len(wanted)} symbols found in the matrix")
    above = [g for g in mapped.values()
             if statistics.median(norm[g]) >= args.cpm_floor]
    print(f"nodes with median CPM >= {args.cpm_floor}: {len(above)}/{len(mapped)}")
    med = statistics.median([statistics.median(norm[g]) for g in mapped.values()]) \
        if mapped else float("nan")
    print(f"median of node median CPM: {med:.2f}")

    print("\nproliferation markers (median CPM, then raw counts per sample):")
    for marker in PROLIFERATION:
        gene = index.get(marker.upper())
        if gene is None:
            continue
        raw = "  ".join(f"{int(v):>5d}" for v in counts[gene])
        print(f"  {gene:8s} CPM={statistics.median(norm[gene]):8.2f}   counts: {raw}")

    print("\nsamples:", ", ".join(samples))
    print("\nRead the raw counts before believing any fold change computed from them:\n"
          "single-digit and low-double-digit counts on libraries of tens of millions\n"
          "of reads do not support a direction, whatever p-value accompanies them.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
