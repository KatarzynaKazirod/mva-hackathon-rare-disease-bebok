"""Download and prepare all reference files for annotating a GRCh38 VCF.

Everything is public data; no patient data is sent anywhere. Run once:

    ./src/download_refs.py                        # everything into refs/
    ./src/download_refs.py --ensembl-release 116  # pin the Ensembl release
    ./src/download_refs.py --update-hpo           # replace HPO files with the latest release
    ./src/download_refs.py --skip-gnomad-af       # without the genome-wide gnomAD AF file
    ./src/download_refs.py --skip-mva-gnomad      # without gnomAD v4.1 for the MVA genes
    ./src/download_refs.py --extra-genes my.tsv   # another file of manually added genes
    ./src/download_refs.py --keep-chr-naming      # leave sources in their own contig naming
    ./src/download_refs.py --keep-original-naming # keep the chr-named 3 GB gnomAD AF file too

What it downloads (all GRCh38):
  1. Ensembl primary-assembly FASTA     -> GRCh38.fa (+ .fai), decompressed
  2. Ensembl gene models, chromosomes    -> Homo_sapiens.GRCh38.<release>.chr.gff3.gz
     (same release as the FASTA)            + mane.tsv (MANE Select transcripts)
                                            + gene_coords.tsv (every gene, with coordinates,
                                              so a phenotype-derived gene list never has to
                                              be sent to a remote API to be resolved)
  3. HPO annotations + ontology          -> hpo/phenotype.hpoa, genes_to_phenotype.txt,
     (pinned GitHub release)                phenotype_to_genes.txt, hp.obo
                                            + mva_genes.tsv (genes linked to MVA in HPO,
                                              plus manual additions from config/mva_genes_extra.tsv)
                                            + mva_checklist.bed (the same padded regions the
                                              gnomAD v4.1 extracts below are cut from)
  4. ClinVar                             -> clinvar.vcf.gz (+ .tbi), MD5-verified
  5. gnomAD allele frequencies (GATK)    -> af-only-gnomad.hg38.vcf.gz (+ .tbi)
  6. gnomAD v4.1 for the MVA genes       -> gnomad_v4.1_{exomes,genomes}.mva.vcf.gz (+ .tbi)
     (only gene coordinates are requested from the gnomAD server)
  7. Ensembl contig naming for every     -> <name>.ensembl.vcf.gz (+ .tbi) for every source
     annotation source                      whose contigs are chr-prefixed

Chromosome naming differs between sources (Ensembl and ClinVar use "1", gnomAD uses "chr1").
`bcftools annotate` joins on the CHROM string, so a mismatch raises no error at all - it
silently leaves the annotation columns empty. Recording the naming is therefore not enough,
and this script does two things with it: manifest.json documents the naming of every file,
and step 7 REWRITES every annotation-source VCF whose contigs are chr-prefixed into an
Ensembl-named copy (<name>.ensembl.vcf.gz, indexed). After a successful run the invariant
for the rest of the pipeline is "everything in refs/ uses Ensembl naming", so the only file
still to convert is the patient VCF itself. Pass --keep-chr-naming to skip the conversion.

Downloads resume after an interruption and finished files are skipped on re-runs.
manifest.json records URLs, releases, dates, sizes, MD5 sums and chromosome naming.
"""

from __future__ import annotations

import argparse
import contextlib
import gzip
import hashlib
import json
import re
import shutil
import sys
import tempfile
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

import pysam
from pysam import bcftools

# ---------------------------------------------------------------- sources

ALLOWED_HOSTS = {
    "ftp.ensembl.org", "ftp.ncbi.nlm.nih.gov", "storage.googleapis.com",
    "github.com", "release-assets.githubusercontent.com", "objects.githubusercontent.com",
}

ENSEMBL_LISTING = "https://ftp.ensembl.org/pub/current_gff3/homo_sapiens/"
ENSEMBL_FASTA = ("https://ftp.ensembl.org/pub/release-{rel}/fasta/homo_sapiens/dna/"
                 "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")
ENSEMBL_GFF = ("https://ftp.ensembl.org/pub/release-{rel}/gff3/homo_sapiens/"
               "Homo_sapiens.GRCh38.{rel}.chr.gff3.gz")

HPO_REPO = "https://github.com/obophenotype/human-phenotype-ontology/releases"
HPO_FILES = ("phenotype.hpoa", "genes_to_phenotype.txt", "phenotype_to_genes.txt", "hp.obo")

CLINVAR = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"

GNOMAD_AF = "https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"

GNOMAD_V4 = {
    "exomes": "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/"
              "gnomad.exomes.v4.1.sites.chr{chrom}.vcf.bgz",
    "genomes": "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/"
               "gnomad.genomes.v4.1.sites.chr{chrom}.vcf.bgz",
}
# Disease-ID allow-list. Name matching alone is NOT sufficient, and this is
# measured rather than assumed: on HPO v2026-09-01 the pattern below matches 6
# diseases yielding 7 genes (BUB1B, CEP57, TRIP13, CENATAC, MAD1L1, and BUB1 +
# BUB3 via ORPHA:1052) and SILENTLY misses SLF2 and SMC5, because OMIM 620184
# and 620185 have the preferred titles "Atelis syndrome 1/2" and MVA5/MVA6
# survive only as alternative titles. CEP192 has no HPO annotation at all and
# is supplied through config/mva_genes_extra.tsv instead.
MVA_DISEASE_IDS = {
    "OMIM:257300",  # MVA1, BUB1B
    "OMIM:614114",  # MVA2, CEP57
    "OMIM:617598",  # MVA3, TRIP13
    "OMIM:620153",  # MVA4 (provisional "?" in OMIM), CENATAC
    "OMIM:620184",  # Atelis syndrome 1 = MVA5 (alternative title), SLF2
    "OMIM:620185",  # Atelis syndrome 2 = MVA6 (alternative title), SMC5
    "OMIM:620189",  # MVA7, MAD1L1
    "ORPHA:1052",   # grouped MVA entry: BUB1B, CEP57, TRIP13, BUB1, BUB3
}
MVA_DISEASE_PATTERN = re.compile(
    r"mosaic variegated (aneuploidy|hyperploidy)|atelis syndrome", re.IGNORECASE)
# The full checklist, not three genes: families have been reported with a
# negative BUB1B/CEP57/TRIP13 screen, so a 3-gene fallback is a silent panel
# collapse in exactly the situation where the fallback fires.
MVA_FALLBACK_GENES = ("BUB1B", "CEP57", "TRIP13", "CENATAC", "SLF2", "SMC5",
                      "MAD1L1", "BUB1", "BUB3", "CEP192")
MVA_PADDING = 5000

#: Repository root: this script lives in src/, and every default path below is anchored
#: here rather than to the current working directory - otherwise `cd src && ./download_refs.py`
#: would quietly create a second refs/ tree inside src/.
ROOT = Path(__file__).resolve().parents[1]

GRCH38_CHR1 = 248_956_422
PRIMARY_CONTIGS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
MIN_FREE_GB = 20          # refs/ settles at ~6.7 GB; the peak during step 7 is ~10.5 GB
BCFTOOLS_THREADS = 4


# ---------------------------------------------------------------- helpers

def check_url(url: str) -> None:
    parts = urllib.parse.urlparse(url)
    if parts.scheme != "https" or parts.hostname not in ALLOWED_HOSTS:
        sys.exit(f"Refusing to use a non-allow-listed source: {url}")


def human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:,.1f} {unit}"
        n /= 1024
    return f"{n:,.1f} TB"


def download(url: str, dest: Path, overwrite: bool = False) -> Path:
    """Download url to dest, resuming a partial .part file; skip if dest exists (unless overwrite)."""
    if dest.exists() and not overwrite:
        print(f"  ✔ already present: {dest.name}")
        return dest
    check_url(url)
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    start = part.stat().st_size if part.exists() else 0
    req = urllib.request.Request(url, headers={"Range": f"bytes={start}-"} if start else {})
    print(f"  ↓ {url}" + (f" (resuming at {human(start)})" if start else ""))
    with urllib.request.urlopen(req, timeout=120) as resp:
        check_url(resp.geturl())                  # the redirect target must be allow-listed too
        if start and resp.status != 206:          # server ignored the range: start over
            start = 0
        total = int(resp.headers.get("Content-Length", 0)) + start
        done, next_report = start, 0.0
        with part.open("ab" if start else "wb") as out:
            while chunk := resp.read(1 << 20):
                out.write(chunk)
                done += len(chunk)
                if total and done / total >= next_report:
                    print(f"    {done / total:6.1%}  {human(done)} / {human(total)}", flush=True)
                    next_report += 0.1
    part.replace(dest)
    return dest


def md5sum(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        while chunk := fh.read(1 << 22):
            h.update(chunk)
    return h.hexdigest()


def file_info(path: Path, previous: dict | None) -> dict:
    """Path, size and MD5; reuse the previous MD5 when the file size has not changed."""
    size = path.stat().st_size
    if previous and previous.get("path") == str(path) and previous.get("size") == size:
        return {k: previous[k] for k in ("path", "size", "md5")}
    return {"path": str(path), "size": size, "md5": md5sum(path)}


def latest_ensembl_release() -> int:
    check_url(ENSEMBL_LISTING)
    with urllib.request.urlopen(ENSEMBL_LISTING, timeout=60) as resp:
        html = resp.read().decode()
    releases = {int(r) for r in re.findall(r"Homo_sapiens\.GRCh38\.(\d+)\.chr\.gff3\.gz", html)}
    if not releases:
        sys.exit("Could not determine the current Ensembl release; pass --ensembl-release.")
    return max(releases)


def latest_hpo_release() -> str:
    """Tag of the latest HPO GitHub release, e.g. 'v2026-09-01' (from the /latest redirect)."""
    url = f"{HPO_REPO}/latest"
    check_url(url)
    with urllib.request.urlopen(url, timeout=60) as resp:
        final = resp.geturl()
    tag = final.rstrip("/").rsplit("/", 1)[-1]
    if not tag.startswith("v"):
        sys.exit(f"Could not determine the latest HPO release from {final}")
    return tag


def chrom_sort_key(chrom: str) -> tuple:
    c = chrom.removeprefix("chr")
    return (0, int(c), "") if c.isdigit() else (1, 0, c)


def gff_attributes(column: str) -> dict[str, str]:
    return dict(kv.split("=", 1) for kv in column.split(";") if "=" in kv)


# ---------------------------------------------------------- contig naming

def ensembl_name(path: Path) -> Path:
    """refs/clinvar.vcf.gz -> refs/clinvar.ensembl.vcf.gz (idempotent)."""
    if ".ensembl.vcf.gz" in path.name:
        return path
    return path.with_name(path.name.replace(".vcf.gz", ".ensembl.vcf.gz"))


def ensembl_rename_map(path: Path) -> dict[str, str]:
    """{chr-name: Ensembl-name} for the primary contigs of a VCF; empty if already Ensembl.

    Only the 25 primary contigs are renamed (chr1->1 ... chrX->X, chrM->MT). Decoys, alt
    and unplaced contigs are deliberately left alone: they cannot match a primary-assembly
    VCF anyway, and renaming them would invent names that exist in no reference file.
    """
    with pysam.VariantFile(str(path)) as vcf:
        contigs = list(vcf.header.contigs)
    mapping: dict[str, str] = {}
    for contig in contigs:
        if not contig.startswith("chr"):
            continue
        short = contig.removeprefix("chr")
        short = "MT" if short == "M" else short
        if short in PRIMARY_CONTIGS:
            mapping[contig] = short
    return mapping


def to_ensembl(src: Path, threads: int = BCFTOOLS_THREADS) -> tuple[Path, bool]:
    """Rewrite src with Ensembl contig naming; return (path to use, converted?).

    A file that is already Ensembl-named yields an empty rename map and is returned
    unchanged, so calling this on every source is idempotent. The rewrite goes through the
    bcftools bundled in pysam - the same binary the annotation pipeline uses.
    """
    dest = ensembl_name(src)
    if dest == src:
        return src, False
    if dest.exists() and Path(f"{dest}.tbi").exists():
        print(f"  \u2714 already present: {dest.name}")
        return dest, True
    mapping = ensembl_rename_map(src)
    if not mapping:
        print(f"  \u00b7 {src.name}: already Ensembl naming, left as is")
        return src, False
    chrmap = dest.with_name(dest.name + ".chrmap.tsv")
    chrmap.write_text("".join(f"{old}\t{new}\n" for old, new in mapping.items()))
    tmp = dest.with_name(dest.name + ".part")
    print(f"  \u21bb {src.name} -> {dest.name} ({len(mapping)} contigs renamed)")
    bcftools.annotate("--rename-chrs", str(chrmap), "--threads", str(threads),
                      "-Oz", "-o", str(tmp), str(src), catch_stdout=False)
    if not tmp.exists() or tmp.stat().st_size == 0:
        tmp.unlink(missing_ok=True)
        sys.exit(f"Renaming the contigs of {src.name} produced no output.")
    tmp.replace(dest)
    pysam.tabix_index(str(dest), preset="vcf", force=True)
    return dest, True


def unify_naming(files: dict[str, Path], naming: dict[str, str],
                 keep_original: bool) -> tuple[dict[str, Path], dict[str, str], dict[str, str]]:
    """Step 7: every annotation source in Ensembl naming, measured rather than assumed."""
    print("\n[7/7] Ensembl contig naming for the annotation sources")
    converted: dict[str, str] = {}
    for key in ("clinvar", "gnomad_af", "gnomad_v4_exomes", "gnomad_v4_genomes"):
        src = files.get(key)
        if not src or not src.exists():
            continue
        path, did_convert = to_ensembl(src)
        files[key] = path
        style, _ = vcf_naming(path)        # re-read the file instead of trusting the rename
        naming[key] = style
        if style != "no_chr":
            print(f"  ! {path.name} still uses '{style}' naming after step 7")
        if did_convert:
            converted[key] = str(path)
        # The 3 GB chr-named original is the one copy worth reclaiming. Deleting it is safe
        # only because step_gnomad_af() looks for the converted file first, so a re-run does
        # not start the 3 GB download all over again.
        if did_convert and key == "gnomad_af" and not keep_original and src != path:
            size = src.stat().st_size
            src.unlink(missing_ok=True)
            Path(f"{src}.tbi").unlink(missing_ok=True)
            print(f"  \u00b7 removed the chr-named original {src.name} ({human(size)} reclaimed)")
    return files, naming, converted


# ---------------------------------------------------------------- steps

def step_fasta(out: Path, rel: int, keep_gz: bool) -> Path:
    print("\n[1/7] Reference genome (Ensembl primary assembly)")
    fasta = out / "GRCh38.fa"
    if fasta.exists():
        if not Path(f"{fasta}.fai").exists():
            print("  the FASTA is present but has no .fai - indexing...")
            pysam.faidx(str(fasta))
        print(f"  ✔ already present: {fasta.name}")
        return fasta
    gz = download(ENSEMBL_FASTA.format(rel=rel), out / "GRCh38.fa.gz")
    print("  decompressing (≈3.1 GB)...")
    tmp = fasta.with_name(fasta.name + ".part")
    with gzip.open(gz, "rb") as src, tmp.open("wb") as dst:
        shutil.copyfileobj(src, dst, length=1 << 22)     # also verifies the gzip CRC
    tmp.replace(fasta)
    if not keep_gz:
        gz.unlink()
    print("  indexing...")
    pysam.faidx(str(fasta))
    return fasta


def step_gff(out: Path, rel: int) -> tuple[Path, Path, dict[str, dict]]:
    """Download the GFF3, write mane.tsv and return coordinates of every gene by name."""
    print("\n[2/7] Gene models (Ensembl GFF3, chromosomes only)")
    gff = download(ENSEMBL_GFF.format(rel=rel), out / f"Homo_sapiens.GRCh38.{rel}.chr.gff3.gz")

    gene_names: dict[str, str] = {}
    coords: dict[str, dict] = {}
    mane_rows = []
    with gzip.open(gff, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attrs = gff_attributes(f[8])
            if attrs.get("ID", "").startswith("gene:"):
                name = attrs.get("Name", attrs["ID"].removeprefix("gene:"))
                gene_names[attrs["ID"]] = name
                coords.setdefault(name, {"chrom": f[0], "start": int(f[3]), "end": int(f[4]),
                                         "strand": f[6], "gene_id": attrs["ID"].removeprefix("gene:")})
            elif "MANE_Select" in attrs.get("tag", ""):
                transcript = attrs.get("ID", "").removeprefix("transcript:")
                mane_rows.append(f"{transcript}\t{gene_names.get(attrs.get('Parent', ''), '')}\n")

    mane = out / "mane.tsv"
    mane.write_text("transcript\tgene\n" + "".join(mane_rows))

    # Every gene with its coordinates, in Ensembl contig naming. The panel step needs the
    # coordinates of a few hundred genes picked from HPO; without this table it would have
    # to re-parse the 103 MB GFF3 to get them, or ask a remote API - and asking a remote
    # API for the coordinates of a phenotype-derived gene list would disclose the
    # phenotype indirectly.
    coords_tsv = out / "gene_coords.tsv"
    coords_tsv.write_text(
        "gene\tchrom\tstart\tend\tstrand\tgene_id\n"
        + "".join(f"{g}\t{c['chrom']}\t{c['start']}\t{c['end']}\t{c['strand']}\t{c['gene_id']}\n"
                  for g, c in sorted(coords.items(),
                                     key=lambda kv: (chrom_sort_key(kv[1]["chrom"]), kv[1]["start"]))))
    print(f"  {len(gene_names):,} genes, {len(mane_rows):,} MANE Select transcripts -> {mane.name}")
    print(f"  gene coordinates -> {coords_tsv.name}")
    return gff, mane, coords_tsv, coords


def step_hpo(out: Path, update: bool, known_tag: str | None) -> tuple[dict[str, Path], str, dict[str, list[str]]]:
    """Download HPO files from one pinned release and find genes linked to MVA."""
    print("\n[3/7] HPO annotations and ontology")
    hpo_dir = out / "hpo"
    present = all((hpo_dir / f).exists() for f in HPO_FILES)
    if present and known_tag and not update:
        tag = known_tag
        print(f"  ✔ HPO release {tag} already present (use --update-hpo for the latest)")
    else:
        tag = latest_hpo_release()
        print(f"  HPO release {tag}")
    paths = {f: download(f"{HPO_REPO}/download/{tag}/{f}", hpo_dir / f,
                         overwrite=update or tag != known_tag) for f in HPO_FILES}

    # sanity check: header and version of phenotype.hpoa
    with paths["phenotype.hpoa"].open() as fh:
        version = next((ln.split(":", 1)[1].strip() for ln in fh if ln.startswith("#version")), "unknown")
    print(f"  phenotype.hpoa version {version}")

    # genes linked to an MVA disease, selected by ID allow-list OR by name
    diseases: dict[str, str] = {}
    with paths["phenotype.hpoa"].open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) > 1 and (f[0] in MVA_DISEASE_IDS
                               or MVA_DISEASE_PATTERN.search(f[1])):
                diseases[f[0]] = f[1]
    # An allow-listed id that is absent must be LOUD: HPO renaming an entry is
    # how SLF2 and SMC5 went missing in the first place, and a quietly shorter
    # gene list is indistinguishable from a correct one.
    if missing := sorted(MVA_DISEASE_IDS - set(diseases)):
        print(f"  ! allow-listed MVA diseases absent from HPO {tag}: {missing}")
        print("    check their current titles before trusting the gene list")
    mva_genes: dict[str, set[str]] = {}
    with paths["genes_to_phenotype.txt"].open() as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > 5 and f[5] in diseases:
                mva_genes.setdefault(f[1], set()).add(f[5])
    print(f"  MVA diseases in HPO: {len(diseases)}")
    for did, name in sorted(diseases.items()):
        genes = sorted(g for g, ds in mva_genes.items() if did in ds)
        print(f"    {did:<13} {name}  →  {', '.join(genes) or '-'}")
    return paths, tag, {g: sorted(ds) for g, ds in mva_genes.items()}


def read_extra_genes(path: Path | None) -> dict[str, dict]:
    """Manually added genes: TSV with columns gene, source, note (lines starting with # are comments)."""
    if path is None or not path.exists():
        return {}
    extras: dict[str, dict] = {}
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t")
        if f[0] == "gene":
            continue
        extras[f[0].strip()] = {"source": f[1].strip() if len(f) > 1 else "",
                                "note": f[2].strip() if len(f) > 2 else ""}
    return extras


def write_mva_genes(out: Path, mva: dict[str, list[str]], coords: dict[str, dict],
                    extras: dict[str, dict]) -> dict[str, dict]:
    """mva_genes.tsv: genes linked to MVA in HPO plus manually added genes, with coordinates.

    Falls back to the classic three genes if HPO found none."""
    if not mva:
        print("  ! No MVA genes found in HPO – using the classic BUB1B, CEP57, TRIP13.")
        mva = {g: [] for g in MVA_FALLBACK_GENES}

    entries: dict[str, dict] = {g: {"origin": "hpo", "diseases": ds, "source": "", "note": ""}
                                for g, ds in mva.items()}
    for g, e in extras.items():
        if g in entries:
            print(f"  · {g} is already in the HPO list – manual entry ignored")
            continue
        entries[g] = {"origin": "manual", "diseases": [], **e}

    table: dict[str, dict] = {}
    for gene in sorted(entries):
        if gene not in coords:
            print(f"  ! {gene} not found in the GFF3 – skipped (check the official gene symbol)")
            continue
        table[gene] = {**coords[gene], **entries[gene]}

    header = "gene\tchrom\tstart\tend\tstrand\tgene_id\torigin\tdiseases\tsource\tnote\n"
    rows = [f"{g}\t{c['chrom']}\t{c['start']}\t{c['end']}\t{c['strand']}\t{c['gene_id']}\t"
            f"{c['origin']}\t{','.join(c['diseases'])}\t{c['source']}\t{c['note']}\n"
            for g, c in sorted(table.items(), key=lambda kv: (chrom_sort_key(kv[1]["chrom"]), kv[1]["start"]))]
    (out / "mva_genes.tsv").write_text(header + "".join(rows))
    n_manual = sum(c["origin"] == "manual" for c in table.values())
    print(f"  {len(table)} MVA genes with coordinates ({len(table) - n_manual} from HPO, "
          f"{n_manual} added manually) -> mva_genes.tsv")
    for g, c in table.items():
        tag = "" if c["origin"] == "hpo" else f"  [manual: {c['source']}]"
        print(f"    {g:<8} {c['chrom']}:{c['start']:,}-{c['end']:,} ({c['strand']}){tag}")
    return table


def step_clinvar(out: Path) -> Path:
    print("\n[4/7] ClinVar")
    vcf = out / "clinvar.vcf.gz"
    fresh = not vcf.exists()
    download(CLINVAR, vcf)
    download(CLINVAR + ".tbi", out / "clinvar.vcf.gz.tbi")
    if fresh:
        md5_file = download(CLINVAR + ".md5", out / "clinvar.vcf.gz.md5")
        expected = md5_file.read_text().split()[0]
        if md5sum(vcf) != expected:
            vcf.unlink()
            sys.exit("ClinVar MD5 mismatch – the file was deleted, please re-run.")
        print("  ✔ MD5 verified")
    return vcf


def step_gnomad_af(out: Path) -> Path:
    print("\n[5/7] gnomAD allele frequencies, genome-wide (GATK af-only resource)")
    converted = ensembl_name(out / "af-only-gnomad.hg38.vcf.gz")
    if converted.exists() and Path(f"{converted}.tbi").exists():
        # Step 7 already rewrote this file and may have deleted the chr-named original;
        # without this guard a re-run would download the 3 GB again.
        print(f"  ✔ already present with Ensembl naming: {converted.name}")
        return converted
    vcf = download(GNOMAD_AF, out / "af-only-gnomad.hg38.vcf.gz")
    download(GNOMAD_AF + ".tbi", out / "af-only-gnomad.hg38.vcf.gz.tbi")
    return vcf


def merge_padded_regions(mva: dict[str, dict]) -> list[list]:
    """Padded gene regions, merged where they overlap. 1-based and inclusive.

    Merging matters twice: no gnomAD record is written twice into the extract, and the BED
    written from the same function covers exactly the intervals that were extracted.
    """
    padded = sorted(((g["chrom"], max(1, g["start"] - MVA_PADDING), g["end"] + MVA_PADDING, name)
                     for name, g in mva.items()), key=lambda r: (chrom_sort_key(r[0]), r[1]))
    regions: list[list] = []
    for chrom, start, end, name in padded:
        if regions and regions[-1][0] == chrom and start <= regions[-1][2]:
            regions[-1][2] = max(regions[-1][2], end)
            regions[-1][3] += f"+{name}"
        else:
            regions.append([chrom, start, end, name])
    return regions


def write_checklist_bed(out: Path, mva: dict[str, dict]) -> Path:
    """BED of the padded checklist regions, for the annotation step's fill check.

    It has to be the SAME intervals the gnomAD v4.1 extracts were cut from: those files
    hold records only inside these regions, so a fill fraction measured anywhere else
    would be a fraction of a percent by construction and would test nothing. BED is
    0-based and half-open, the regions above are 1-based and inclusive, hence start - 1.
    """
    dest = out / "mva_checklist.bed"
    regions = merge_padded_regions(mva)
    dest.write_text("".join(f"{c}\t{s - 1}\t{e}\t{n}\n" for c, s, e, n in regions))
    span = sum(e - (s - 1) for _, s, e, _ in regions)
    print(f"  {len(regions)} merged regions, {span / 1e6:.2f} Mb -> {dest.name} "
          f"(padding {MVA_PADDING:,} bp)")
    return dest


def step_mva_gnomad(out: Path, mva: dict[str, dict], refresh: bool) -> dict[str, str]:
    print(f"\n[6/7] gnomAD v4.1 for the MVA genes ({', '.join(mva)}; gene regions only)")
    status: dict[str, str] = {}
    regions = merge_padded_regions(mva)
    for dataset, template in GNOMAD_V4.items():
        dest = (out / f"gnomad_v4.1_{dataset}.mva.vcf.gz").resolve()
        if dest.exists() and Path(f"{dest}.tbi").exists() and not refresh:
            print(f"  ✔ already present: {dest.name}")
            status[dataset] = "ok"
            continue
        try:
            # htslib writes remote index files to the working directory, so use a temporary one
            with tempfile.TemporaryDirectory() as tmp, contextlib.chdir(tmp):
                # Open every chromosome we need FIRST and merge their headers.
                # Three things make this necessary, all of them verified by test:
                #   * a per-chromosome header may declare only its own contig;
                #   * headers may declare different INFO sets, and writing a
                #     record carrying an INFO tag the output header lacks fails
                #     with "Can't write record: Invalid argument";
                #   * header.merge() unions contigs AND INFO/FORMAT/FILTER.
                # Opening once per chromosome instead of once per region also
                # avoids re-opening the same remote file for neighbouring genes.
                sources: dict[str, pysam.VariantFile] = {}
                for chrom, _, _, _ in regions:
                    short = chrom.removeprefix("chr")
                    if short not in sources:
                        url = template.format(chrom=short)
                        check_url(url)
                        sources[short] = pysam.VariantFile(url)
                header = next(iter(sources.values())).header.copy()
                for short, src in list(sources.items())[1:]:
                    header.merge(src.header)
                print(f"  {dataset}: merged headers of {len(sources)} chromosome "
                      f"file(s) -> {len(list(header.contigs))} contigs, "
                      f"{len(list(header.info))} INFO fields")
                writer = pysam.VariantFile(str(dest), "wz", header=header)
                try:
                    for chrom, start, end, name in regions:
                        short = chrom.removeprefix("chr")
                        region = f"chr{short}:{start}-{end}"
                        n = 0
                        for rec in sources[short].fetch(region=region):
                            # translate() rebinds the record to the output
                            # header IN PLACE and returns None, so it cannot be
                            # inlined into write(): write(rec.translate(h))
                            # passes None and raises "record must not be None".
                            rec.translate(writer.header)
                            writer.write(rec)
                            n += 1
                        print(f"  {dataset}: {name} {region} -> {n:,} records")
                finally:
                    writer.close()
                    for src in sources.values():
                        src.close()
            pysam.tabix_index(str(dest), preset="vcf", force=True)
            status[dataset] = "ok"
        except (OSError, ValueError) as err:
            dest.unlink(missing_ok=True)
            print(f"  ! {dataset}: extraction failed ({err}). Continuing without it.")
            status[dataset] = f"failed: {err}"
    return status


# ---------------------------------------------------------------- checks

def vcf_naming(path: Path) -> tuple[str, str | None]:
    """('chr' | 'no_chr', first record) for a VCF, judged from header contigs or the first record."""
    with pysam.VariantFile(str(path)) as vcf:
        contigs = list(vcf.header.contigs)
        first = next(iter(vcf), None)
    names = contigs or ([first.chrom] if first else [])
    naming = "chr" if any(c.startswith("chr") for c in names) else "no_chr"
    return naming, f"{first.chrom}:{first.pos}" if first else None


def verify(files: dict[str, Path]) -> tuple[dict[str, str], list[str]]:
    """Genome build, readability and chromosome naming of every file."""
    print("\nChecks")
    naming: dict[str, str] = {}
    problems: list[str] = []

    if fasta := files.get("fasta"):
        with pysam.FastaFile(str(fasta)) as fa:
            refs = list(fa.references)
            # try both namings, so a chr-prefixed FASTA reports its real length
            # instead of a misleading "chromosome 1 length None"
            chr1 = next((r for r in ("1", "chr1") if r in refs), None)
            length = fa.get_reference_length(chr1) if chr1 else None
        naming["fasta"] = "chr" if any(r.startswith("chr") for r in refs) else "no_chr"
        if length == GRCH38_CHR1:
            print(f"  ✔ FASTA: GRCh38 ({chr1} = {length:,} bp), naming '{naming['fasta']}'")
        elif chr1 is None:
            print(f"  ✘ FASTA: no contig named '1' or 'chr1' among {len(refs)} references")
            problems.append("FASTA has no chromosome 1 contig")
        else:
            # The case that matters most: a GRCh37 FASTA also HAS a '1' contig,
            # just a different length (249,250,621 bp). This branch is the whole
            # point of the check, so it must never be an unhandled fall-through.
            print(f"  ✘ FASTA: {chr1} is {length:,} bp – expected {GRCH38_CHR1:,} "
                  f"for GRCh38 (GRCh37 chromosome 1 is 249,250,621 bp)")
            problems.append(f"FASTA is not GRCh38: {chr1} = {length:,} bp")

    if gff := files.get("gff3"):
        with gzip.open(gff, "rt") as fh:
            seqid = next((ln.split("\t", 1)[0] for ln in fh if not ln.startswith("#")), "")
        naming["gff3"] = "chr" if seqid.startswith("chr") else "no_chr"
        print(f"  ✔ GFF3: naming '{naming['gff3']}'")

    for key in ("clinvar", "gnomad_af", "gnomad_v4_exomes", "gnomad_v4_genomes"):
        path = files.get(key)
        if not path or not path.exists():
            continue
        name_style, first = vcf_naming(path)
        naming[key] = name_style
        if first:
            print(f"  ✔ {path.name}: readable, naming '{name_style}', first record {first}")
        else:
            print(f"  ✘ {path.name}: no records")
            problems.append(f"{path.name} has no records")

    styles = set(naming.values())
    if len(styles) > 1:
        groups = {s: [k for k, v in naming.items() if v == s] for s in styles}
        print("  ! Mixed chromosome naming: "
              + "; ".join(f"'{s}': {', '.join(ks)}" for s, ks in groups.items())
              + "\n    The annotation step must rename chromosomes per source and check that each annotation landed.")
    return naming, problems


# ---------------------------------------------------------------- main

def main() -> None:
    p = argparse.ArgumentParser(description="Download reference files for GRCh38 VCF annotation.")
    p.add_argument("--out", type=Path, default=ROOT / "refs",
                   help="output directory (default <repo>/refs)")
    p.add_argument("--ensembl-release", type=int, help="Ensembl release for FASTA + GFF3 (default: current)")
    p.add_argument("--update-hpo", action="store_true", help="replace HPO files with the latest release")
    p.add_argument("--skip-gnomad-af", action="store_true", help="skip the genome-wide gnomAD AF file")
    p.add_argument("--skip-mva-gnomad", action="store_true", help="skip gnomAD v4.1 for the MVA genes")
    p.add_argument("--keep-gz", action="store_true", help="keep the compressed FASTA after decompressing")
    p.add_argument("--keep-chr-naming", action="store_true",
                   help="skip step 7: leave every source in its own contig naming")
    p.add_argument("--keep-original-naming", action="store_true",
                   help="in step 7 also keep the chr-named copy of the 3 GB gnomAD AF file")
    p.add_argument("--extra-genes", type=Path, default=ROOT / "config" / "mva_genes_extra.tsv",
                   help="TSV of manually added MVA-related genes "
                        "(default <repo>/config/mva_genes_extra.tsv)")
    args = p.parse_args()

    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    free_gb = shutil.disk_usage(out).free / 1024**3
    print(f"Output: {out.resolve()}  (free space {free_gb:,.0f} GB)")
    if free_gb < MIN_FREE_GB:
        print(f"! Less than {MIN_FREE_GB} GB free – the downloads may not fit.")

    manifest_path = out / "manifest.json"
    old = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    rel = args.ensembl_release or old.get("ensembl_release") or latest_ensembl_release()
    if old.get("ensembl_release") and old["ensembl_release"] != rel:
        sys.exit(f"{out}/ already holds Ensembl release {old['ensembl_release']}; "
                 f"use another --out for release {rel}.")
    print(f"Ensembl release: {rel}")

    files: dict[str, Path] = {}
    files["fasta"] = step_fasta(out, rel, args.keep_gz)
    files["gff3"], files["mane"], files["gene_coords"], coords = step_gff(out, rel)

    hpo_paths, hpo_tag, mva_hpo = step_hpo(out, args.update_hpo, old.get("hpo_release"))
    for name, path in hpo_paths.items():
        files[f"hpo_{name}"] = path
    extras = read_extra_genes(args.extra_genes)
    if args.extra_genes and not args.extra_genes.exists():
        print(f"  ! {args.extra_genes} not found – no manually added genes")
    mva = write_mva_genes(out, mva_hpo, coords, extras)
    files["mva_genes"] = out / "mva_genes.tsv"
    if mva:
        files["mva_checklist_bed"] = write_checklist_bed(out, mva)
    genes_changed = sorted(mva) != sorted(old.get("mva_genes", {}))

    files["clinvar"] = step_clinvar(out)
    if not args.skip_gnomad_af:
        files["gnomad_af"] = step_gnomad_af(out)
    mva_status: dict[str, str] = {}
    if not args.skip_mva_gnomad and mva:
        mva_status = step_mva_gnomad(out, mva, refresh=genes_changed)
        for dataset, state in mva_status.items():
            if state == "ok":
                files[f"gnomad_v4_{dataset}"] = out / f"gnomad_v4.1_{dataset}.mva.vcf.gz"

    naming, problems = verify(files)
    ensembl_converted: dict[str, str] = {}
    if args.keep_chr_naming:
        print("\n[7/7] Ensembl contig naming - skipped (--keep-chr-naming)")
    else:
        files, naming, ensembl_converted = unify_naming(files, naming, args.keep_original_naming)

    print("\nComputing MD5 sums for the manifest...")
    old_files = old.get("files", {})
    manifest = {
        "created": old.get("created", date.today().isoformat()),
        "updated": date.today().isoformat(),
        "ensembl_release": rel,
        "hpo_release": hpo_tag,
        "sources": {
            "fasta": ENSEMBL_FASTA.format(rel=rel), "gff3": ENSEMBL_GFF.format(rel=rel),
            "hpo": f"{HPO_REPO}/tag/{hpo_tag}", "clinvar": CLINVAR, "gnomad_af": GNOMAD_AF,
            "gnomad_v4_mva": GNOMAD_V4,
        },
        "chrom_naming": naming,
        "ensembl_converted": ensembl_converted,
        "mva_genes": mva, "mva_padding": MVA_PADDING, "mva_gnomad_status": mva_status,
        "files": {k: file_info(v, old_files.get(k)) for k, v in files.items() if v.exists()},
    }
    manifest_path.write_text(json.dumps(manifest, indent=2))

    print(f"\nDone. Manifest: {manifest_path}")
    for key, info in manifest["files"].items():
        print(f"  {key:<28} {Path(info['path']).name:<40} {human(info['size'])}")
    if problems:
        print("\n! Problems:\n  " + "\n  ".join(problems))
        sys.exit(1)


if __name__ == "__main__":
    main()