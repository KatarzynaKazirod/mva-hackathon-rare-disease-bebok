#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12,<3.14"
# dependencies = ["pysam>=0.24,<0.25"]
# ///
"""
annotate_vcf.py

1. Checks the chromosome naming convention and the genome build of every input file.
2. Harmonises EVERY VCF that takes part in a coordinate join - the sample and both
   annotation sources - along the same path: keep only contigs present in the FASTA,
   rename to Ensembl-style names (chr1 -> 1, chrM -> MT), then `norm -m -any -f`.
3. Runs the annotation pipeline with bcftools:
     csq -l -> per-variant functional consequences (BCSQ) from GFF3 + FASTA
     annotate (ClinVar) -> CLNSIG, CLNDN, CLNREVSTAT, ALLELEID
     annotate (gnomAD)  -> gnomAD_AF
4. Verifies that every new INFO field actually landed and FAILS if it did not.
5. Writes annotated.vcf.gz (+ index) and annotated.tsv into --workdir, and
   annotation_summary.json - counts, fill fractions, tool and data versions - into
   --outdir. The summary deliberately contains no patient file names.

Uses the system `bcftools` if it is on PATH, otherwise the copy bundled in pysam.

Index policy: whatever index the data came with is used as-is. Its modification time is
never examined, because mtime cannot tell a genuinely stale index apart from the order the
files happened to be copied in. A file with no index at all gets one built next to it.

Step caching is NOT mtime-based. Each step writes a .stamp.json next to its output holding
the resolved path, size and mtime of every input plus the exact argv, and the stamp is
written last, so an interrupted write leaves no usable cache. Changing the input to an
older file, or swapping in another sample copied with `cp -p`, invalidates the stamp -
which plain mtime comparison would not.

Run a smoke test on a single chromosome before committing hours to the whole genome:

    ./src/annotate_vcf.py --regions 21 --workdir PATIENT_DATA/smoke --no-tsv ...

Example:
    ./src/annotate_vcf.py \\
        --vcf PATIENT_DATA/raw/<your_sample>.vcf.gz \\
        --fasta refs/GRCh38.fa \\
        --gff refs/Homo_sapiens.GRCh38.116.chr.gff3.gz \\
        --clinvar refs/clinvar.vcf.gz \\
        --gnomad refs/af-only-gnomad.hg38.ensembl.vcf.gz \\
        --workdir PATIENT_DATA/work --outdir results --cache refs
"""

import argparse
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pysam

PRIMARY = re.compile(r"^(chr)?([0-9]{1,2}|X|Y|M|MT)$", re.IGNORECASE)
GRCH38_CHR1 = 248956422
GRCH37_CHR1 = 249250621
#: ALLELEID is carried along so a candidate can be traced back to its ClinVar record.
CLINVAR_FIELDS = ["CLNSIG", "CLNDN", "CLNREVSTAT", "ALLELEID"]
BCFTOOLS_BIN = shutil.which("bcftools")

#: Repository root. The script lives in src/, so the default output directories are
#: anchored here rather than to the current working directory.
ROOT = Path(__file__).resolve().parents[1]

#: Minimum fraction of records that must carry each new INFO field. `bcftools annotate`
#: exits 0 having matched ZERO records, so these are the only thing between a silent
#: naming mismatch and an answer read off empty columns.
#:
#: Calibrated against a measured run, not guessed: on chromosome 21 of a single-sample WGS
#: (87,337 records after normalisation) the observed fill was BCSQ 61.5%, CLNSIG 0.7%,
#: gnomAD_AF 71.7%. The BCSQ floor started at 5% on the assumption that csq reports only
#: coding consequences - it does not: `intron` and `non_coding` count too, and gene bodies
#: plus lncRNA/pseudogene transcripts cover most of the genome. A floor of 5% against a
#: normal value of 60% can only catch a total failure, so it is 30%: low enough for a
#: gene-poor region, high enough to catch a GFF3 that fits only part of the data.
#: One chromosome of one sample is calibration, not validation - revisit if a legitimate
#: run trips it.
MIN_FILL = {"BCSQ": 0.30, "CLNSIG": 0.002, "gnomAD_AF": 0.10}

#: Records dropped by `norm --check-ref x` as a fraction of the input. With a hard build
#: check upstream this should be zero; anything above noise means the VCF and the FASTA
#: are not the same assembly and no coordinate join downstream can be trusted.
REF_MISMATCH_MAX_FRAC = 1e-4

#: ##INFO lines for fields we RENAME. A renamed tag is declared in no source header, so
#: annotate has to be handed the definition with -h. Fields copied under their own name
#: (CLNSIG and friends) need no entry - their definition travels with them from ClinVar.
INFO_HEADERS = {
    "gnomAD_AF": '##INFO=<ID=gnomAD_AF,Number=A,Type=Float,Description="Alternate allele '
                 'frequency in gnomAD (af-only resource; GLOBAL AF only, no popmax and no '
                 'nhomalt; a missing value means absent from THIS resource, not AF=0)">',
}


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def log(msg):
    print(f"[annotate] {msg}", file=sys.stderr, flush=True)


def warn(msg):
    print(f"[annotate] WARNING: {msg}", file=sys.stderr, flush=True)


def die(msg):
    print(f"[annotate] ERROR: {msg}", file=sys.stderr, flush=True)
    sys.exit(1)


def bcftools(*args, out=None):
    """Run a bcftools subcommand, via the binary or via pysam.

    `out`, when given, is checked afterwards. That check is not decoration: while pysam
    captures stdout (its default) a `-o FILE` argument is silently ignored, no file appears
    and no exception is raised. Every call here passes catch_stdout=False, and the check
    makes a regression impossible to miss.
    """
    args = [str(a) for a in args]
    log("bcftools " + " ".join(args))
    if BCFTOOLS_BIN:
        subprocess.run([BCFTOOLS_BIN, *args], check=True)
    else:
        from pysam import bcftools as pb
        getattr(pb, args[0])(*args[1:], catch_stdout=False)
    if out is not None:
        out = Path(out)
        if not out.exists() or out.stat().st_size == 0:
            die(f"bcftools {args[0]} exited 0 but produced no output at {out}")


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def is_bgzf(path):
    with open(path, "rb") as fh:
        return fh.read(4) == b"\x1f\x8b\x08\x04"


def has_index(path):
    return any(Path(f"{path}{ext}").exists() for ext in (".tbi", ".csi"))


def index_vcf(path):
    """Index a bgzipped VCF (.tbi) or a BCF (.csi). -t is tabix-only, hence the split."""
    path = Path(path)
    flags = ["-f"] if path.name.endswith(".bcf") else ["-f", "-t"]
    bcftools("index", *flags, path)


def fingerprint(paths, argv):
    """What a cached step was built FROM, not merely when.

    mtime alone is not an answer to "is this output current": pointing an argument at an
    older file, or swapping in another sample copied with `cp -p` / `rsync -a`, leaves the
    output looking fresh. Resolved path catches a changed argument, size catches a
    different file under the same name, mtime_ns catches an edit in place, and argv catches
    a changed flag - including --regions, where a stale cache would silently turn one
    chromosome into a "genome-wide" result.
    """
    return {
        "argv": [str(a) for a in argv],
        "inputs": sorted(
            [str(Path(p).resolve()), Path(p).stat().st_size, Path(p).stat().st_mtime_ns]
            for p in paths
        ),
    }


def cached_step(label, out, inputs, build_args, index=True):
    """Run one bcftools step unless its stamp proves the output is already current.

    The output is written to <out>.part and moved into place with os.replace, and the
    stamp is written LAST, so it acts as a commit marker: an interrupted or killed write
    leaves a stale .part and no stamp, never a truncated file that the next run would
    accept as cached. That matters because htslib reads a truncated BGZF with nothing worse
    than a `no BGZF EOF marker` warning.
    """
    out, stamp = Path(out), Path(f"{out}.stamp.json")
    want = fingerprint(inputs, build_args(out))
    if out.exists() and stamp.exists() and (not index or has_index(out)):
        try:
            if json.loads(stamp.read_text()) == want:
                log(f"{label}: using cached {out.name}")
                return out
        except (json.JSONDecodeError, OSError):
            pass
    stamp.unlink(missing_ok=True)
    tmp = Path(f"{out}.part")
    tmp.unlink(missing_ok=True)
    bcftools(*build_args(tmp), out=tmp)
    for ext in (".tbi", ".csi"):
        Path(f"{out}{ext}").unlink(missing_ok=True)
    os.replace(tmp, out)
    if index:
        index_vcf(out)
    stamp.write_text(json.dumps(want))          # commit marker, written last
    return out


def write_atomic(dest, writer):
    """Same commit discipline for the files we write from Python (FASTA, GFF3)."""
    dest = Path(dest)
    tmp = Path(f"{dest}.part")
    writer(tmp)
    os.replace(tmp, dest)
    return dest


# --------------------------------------------------------------------------- #
# counting (for the fill checks)
# --------------------------------------------------------------------------- #
def _count_query(vcf, expr=None):
    """Stream one line per matching record into a temp file and count the lines.

    Counted from a file rather than from memory so that five million hits never have to
    fit in RAM, and without `bcftools index -n`, which needs count metadata that indexes
    written by older tabix versions do not carry.
    """
    args = ["query"]
    if expr:
        args += ["-i", expr]
    args += ["-f", "1\n"]
    with tempfile.TemporaryDirectory() as tmp:
        hits = Path(tmp) / "hits.txt"
        bcftools(*args, "-o", hits, vcf, out=hits)
        total = 0
        with hits.open("rb") as fh:
            while chunk := fh.read(1 << 22):
                total += chunk.count(b"\n")
    return total


def bcftools_stdout(*args):
    """Run a bcftools subcommand and return its stdout (small outputs only)."""
    args = [str(a) for a in args]
    if BCFTOOLS_BIN:
        return subprocess.run([BCFTOOLS_BIN, *args], check=True,
                              capture_output=True, text=True).stdout
    from pysam import bcftools as pb
    return getattr(pb, args[0])(*args[1:]) or ""


def count_records(vcf):
    """Record count: from the index when it can tell, otherwise by streaming.

    Every file this script writes gets a fresh index that carries count metadata, so the
    cheap path is taken for all intermediates; an input index written by an older tabix
    has no counts and forces the streaming count. The fallback must trigger on ANY failure
    of `index -n`: subprocess.CalledProcessError stringifies to the exit status alone and
    never carries the stderr text, so matching on the message is not enough to notice.
    """
    try:
        if BCFTOOLS_BIN:
            proc = subprocess.run([BCFTOOLS_BIN, "index", "-n", str(vcf)],
                                  capture_output=True, text=True)
            if proc.returncode == 0:
                return int(proc.stdout.strip())
        else:
            from pysam import bcftools as pb
            return int((pb.index("-n", str(vcf)) or "").strip())
    except Exception:
        pass
    return _count_query(vcf)


def count_info_present(vcf, key):
    """Records that actually carry INFO/<key>.

    In a bcftools expression a missing value equals ".", so `INFO/KEY!="."` is the
    documented way of asking whether a tag is present at all.
    """
    return _count_query(vcf, f'INFO/{key}!="."')


def assert_landed(vcf, key, floor, total):
    """Fail unless at least `floor` of the records carry INFO/<key>."""
    present = count_info_present(vcf, key)
    frac = present / total if total else 0.0
    log(f"  {key}: {present:,} / {total:,} ({frac:.1%})")
    if total == 0:
        die("the annotated file has no records at all")
    if frac < floor:
        die(f"INFO/{key} landed on only {frac:.2%} of records, expected at least "
            f"{floor:.1%}.\n"
            "  Check, in this order: contig naming of the annotation source, whether BOTH "
            "sides were normalised, and the -c column map.")
    return {"present": present, "total": total, "fraction": round(frac, 4)}


# --------------------------------------------------------------------------- #
# chromosome naming
# --------------------------------------------------------------------------- #
def to_numeric(name):
    """chr1 -> 1, chrX -> X, chrM/M -> MT. Non-standard contigs are left as-is."""
    m = PRIMARY.match(name)
    if not m:
        return name
    core = m.group(2).upper()
    return "MT" if core in ("M", "MT") else core


def needs_rename(names):
    return any(to_numeric(n) != n for n in names)


def describe(names):
    primary = [n for n in names if PRIMARY.match(n)]
    if not primary:
        return "no standard chromosome names found"
    n_chr = sum(n.lower().startswith("chr") for n in primary)
    if n_chr == 0:
        return "numeric (1, 2, ..., X, Y, MT)"
    if n_chr == len(primary):
        return "UCSC-style (chr1, chr2, ..., chrM)"
    return "mixed"


def report(label, names):
    log(f"{label}: {describe(names)} | e.g. {', '.join(names[:5])}")
    if not any(PRIMARY.match(n) for n in names):
        warn(f"{label}: could not recognise any chromosome names")


def write_map(path):
    lines = [f"chr{c}\t{c}" for c in list(map(str, range(1, 23))) + ["X", "Y"]]
    lines += ["chrM\tMT", "chrMT\tMT", "M\tMT"]
    Path(path).write_text("\n".join(lines) + "\n")


def write_header_file(keys, dest):
    """##INFO definitions for the renamed fields, handed to `annotate -h`."""
    Path(dest).parent.mkdir(parents=True, exist_ok=True)
    Path(dest).write_text("".join(INFO_HEADERS[k] + "\n" for k in keys))
    return dest


# --------------------------------------------------------------------------- #
# contig discovery
# --------------------------------------------------------------------------- #
def vcf_contigs(path, max_records=200_000):
    """Contigs with data (from the index) if indexed, else header / record scan.

    Asking the index rather than the header matters twice: the header of a GATK VCF lists
    thousands of contigs that carry no records, and contig names live in the core of the
    index, not in the optional count metadata that older tabix versions omit.
    """
    if has_index(path):
        try:
            tbx = pysam.TabixFile(str(path))
            names = list(tbx.contigs)
            tbx.close()
            if names:
                return names
        except (OSError, ValueError):
            pass
    with pysam.VariantFile(str(path)) as vf:
        names = list(vf.header.contigs)
        if names:
            return names
        seen = []
        for i, rec in enumerate(vf):
            if rec.chrom not in seen:
                seen.append(rec.chrom)
            if i >= max_records:
                break
        return seen


def fasta_contigs(path):
    fai = Path(f"{path}.fai")
    if not fai.exists():
        log(f"indexing {path}")
        pysam.faidx(str(path))
    with open(fai) as fh:
        return [line.split("\t", 1)[0] for line in fh if line.strip()]


def fai_lengths(fasta):
    lengths = {}
    with open(f"{fasta}.fai") as fh:
        for line in fh:
            if line.strip():
                name, length = line.split("\t")[:2]
                lengths[name] = int(length)
    return lengths


def gff_contigs(path, max_lines=500_000):
    names, seen = [], set()
    with open_text(path) as fh:
        for i, line in enumerate(fh):
            if i >= max_lines:
                break
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) < 2:
                    continue
                name = parts[1]
            elif line.startswith("#") or not line.strip():
                continue
            else:
                name = line.split("\t", 1)[0]
            if name not in seen:
                seen.add(name)
                names.append(name)
    return names


def gff_sequence_regions(path, max_lines=100_000):
    """{contig: length} from the GFF3's ##sequence-region pragmas, for the build check."""
    out = {}
    with open_text(path) as fh:
        for i, line in enumerate(fh):
            if i >= max_lines or (line.startswith("#") is False and out):
                break
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        out[parts[1]] = int(parts[3])
                    except ValueError:
                        pass
    return out


# --------------------------------------------------------------------------- #
# harmonising every VCF that takes part in a coordinate join
# --------------------------------------------------------------------------- #
def harmonise_vcf(path, label, fasta, dest_dir, mapfile, threads, wanted=None,
                  do_norm=True):
    """One path for the sample AND for every annotation source.

    `bcftools annotate -a` joins literally on CHROM, POS, REF and ALT, so both sides have
    to be in the same representation. Harmonising only the target is a silent data loss:
    a multiallelic source record (ACT -> A,ACTCT) never meets the split, left-aligned
    target record it corresponds to. It is a fraction of a percent of records, invisible to
    any fill threshold - and it is exactly the indels where a missing frequency then reads
    as "rare".

    Three stages, each cached and atomic:
      1. bgzip + index if needed;
      2. keep only contigs the FASTA has (plus --regions), and rename to Ensembl naming -
         `annotate -R <bed> --rename-chrs <map>` is the only subcommand that does both in
         one pass, and `norm -f` aborts on a contig the FASTA lacks;
      3. `norm -m -any -f` - split multiallelics, left-align indels, drop REF mismatches.
    """
    path = ensure_bgzf_indexed(path, label, dest_dir, threads)
    names = vcf_contigs(path)
    report(label, names)

    fasta_len = fai_lengths(fasta)
    keep = [n for n in names if to_numeric(n) in fasta_len]
    drop = [n for n in names if to_numeric(n) not in fasta_len]
    deselected = []
    if wanted is not None:
        unknown = wanted - {to_numeric(n) for n in names}
        if unknown:
            warn(f"{label}: --regions names contig(s) this file does not have: "
                 f"{', '.join(sorted(unknown))} - check for a typo, they are ignored")
        deselected = [n for n in keep if to_numeric(n) not in wanted]
        keep = [n for n in keep if to_numeric(n) in wanted]

    if not keep:
        die(f"{label}: none of its contigs are present in the FASTA (or survived "
            "--regions) - check the genome build / reference file")
    if drop:
        shown = ", ".join(drop[:10]) + (" ..." if len(drop) > 10 else "")
        log(f"{label}: dropping records on {len(drop)} contig(s) absent from the FASTA: "
            f"{shown}")
    if deselected:
        shown = ", ".join(deselected[:10]) + (" ..." if len(deselected) > 10 else "")
        log(f"{label}: {len(deselected)} contig(s) excluded by --regions: {shown}")
    log(f"{label}: keeping {len(keep)} contig(s): {', '.join(keep[:10])}"
        f"{' ...' if len(keep) > 10 else ''}")

    current = path
    if drop or deselected or needs_rename(keep):
        bed = dest_dir / f"{label}.keep.bed"
        content = "".join(f"{n}\t0\t{fasta_len[to_numeric(n)]}\n" for n in keep)
        if not (bed.exists() and bed.read_text() == content):
            bed.write_text(content)
        subset = dest_dir / f"{label}.subset.vcf.gz"
        rename = ["--rename-chrs", str(mapfile)] if needs_rename(keep) else []
        current = cached_step(
            label, subset, [current, bed],
            lambda o, src=current, rn=rename: ["annotate", "--threads", threads, "-R", bed,
                                               *rn, "-Oz", "-o", o, src])

    stats = {}
    if do_norm:
        before = count_records(current)
        normed = dest_dir / f"{label}.norm.bcf"
        src = current
        current = cached_step(
            label + " norm", normed, [src, fasta],
            lambda o, s=src: ["norm", "--threads", threads, "-m", "-any", "-f", fasta,
                              "--check-ref", "x", "-Ob", "-o", o, s])
        after = count_records(current)
        # --check-ref x EXCLUDES a record whose REF disagrees with the reference, instead
        # of passing it through with a warning that pysam then discards. Excluded records
        # are therefore counted here rather than hoped about.
        dropped = max(0, before - after)
        stats = {"records_in": before, "records_out": after,
                 "ref_mismatch_dropped": dropped}
        log(f"{label}: norm {before:,} -> {after:,} records "
            f"({dropped:,} dropped by --check-ref x)")
        if before and dropped / before > REF_MISMATCH_MAX_FRAC:
            die(f"{label}: {dropped:,} records ({dropped / before:.3%}) disagree with the "
                "reference base. Above noise this means a different assembly or FASTA "
                "flavour, and no coordinate join can be trusted.")
    return current, stats


def ensure_bgzf_indexed(path, label, cache, threads):
    path = Path(path)
    if not is_bgzf(path):
        bgz = cache / f"{label}.bgz.vcf.gz"
        current = cached_step(label, bgz, [path],
                              lambda o, s=path: ["view", "--threads", threads,
                                                 "-Oz", "-o", o, s])
        return current
    if not has_index(path):
        index_vcf(path)
    return path


def prepare_fasta(path, cache):
    path = Path(path)
    names = fasta_contigs(path)
    report("fasta", names)
    if not needs_rename(names):
        return path

    renamed = cache / "reference.numeric.fa"
    if renamed.exists() and Path(f"{renamed}.fai").exists():
        log("fasta: using cached reference.numeric.fa")
        return renamed

    log("fasta: converting chromosome names to numeric (this can take a few minutes)")

    def _write(tmp):
        with open_text(path) as fi, open(tmp, "w") as fo:
            for line in fi:
                if line.startswith(">"):
                    name, _, rest = line[1:].rstrip("\n").partition(" ")
                    fo.write(">" + to_numeric(name) + (f" {rest}" if rest else "") + "\n")
                else:
                    fo.write(line)

    write_atomic(renamed, _write)
    pysam.faidx(str(renamed))
    return renamed


def prepare_gff(path, cache):
    path = Path(path)
    names = gff_contigs(path)
    report("gff3", names)
    if not needs_rename(names):
        return path

    renamed = cache / "genes.numeric.gff3.gz"
    if renamed.exists():
        log("gff3: using cached genes.numeric.gff3.gz")
        return renamed

    log("gff3: converting chromosome names to numeric")

    def _write(tmp):
        with open_text(path) as fi, gzip.open(tmp, "wt") as fo:
            for line in fi:
                if line.startswith("##sequence-region"):
                    parts = line.split()
                    if len(parts) > 1:
                        parts[1] = to_numeric(parts[1])
                    fo.write(" ".join(parts) + "\n")
                elif line.startswith("#") or not line.strip():
                    fo.write(line)
                else:
                    seqid, sep, rest = line.partition("\t")
                    fo.write(to_numeric(seqid) + sep + rest)

    return write_atomic(renamed, _write)


# --------------------------------------------------------------------------- #
# sanity checks
# --------------------------------------------------------------------------- #
def check_build(fasta, gff, vcfs):
    """A GRCh37 FASTA also has a contig named '1' - it is just 249,250,621 bp long."""
    chr1 = fai_lengths(fasta).get("1")
    if chr1 == GRCH37_CHR1:
        die(f"the FASTA is GRCh37/hg19, not GRCh38 (chromosome 1 is {GRCH37_CHR1:,} bp; "
            f"GRCh38 is {GRCH38_CHR1:,} bp). Nothing downstream could be trusted.")
    if chr1 is None:
        die("the FASTA has no contig named '1' after harmonisation")
    if chr1 != GRCH38_CHR1:
        die(f"unexpected length of chromosome 1 in the FASTA ({chr1:,}); GRCh38 is "
            f"{GRCH38_CHR1:,} bp")

    # The GFF3 states its own assembly in ##sequence-region, and a GFF3 from the wrong
    # release is the one build mismatch the FASTA check cannot see.
    gff_len = gff_sequence_regions(gff).get("1")
    if gff_len and gff_len != chr1:
        die(f"the GFF3 says chromosome 1 is {gff_len:,} bp, the FASTA says {chr1:,} bp - "
            "they are different assemblies")
    if gff_len:
        log(f"gff3: ##sequence-region 1 = {gff_len:,} bp, matches the FASTA")

    for label, path in vcfs:
        with pysam.VariantFile(str(path)) as vf:
            if "1" in vf.header.contigs:
                length = vf.header.contigs["1"].length
                if length and length != chr1:
                    die(f"{label}: chromosome 1 is {length:,} bp against the FASTA's "
                        f"{chr1:,} bp - a different genome build")


def check_not_gvcf(path):
    """A gVCF would fail the fill thresholds with a misleading message about naming."""
    with pysam.VariantFile(str(path)) as vf:
        if "NON_REF" in str(vf.header) or "<NON_REF>" in str(vf.header):
            die("this looks like a gVCF (<NON_REF> is declared in the header). Reference "
                "blocks cannot be annotated by a coordinate join - genotype it first "
                "(GenotypeGVCFs) and annotate the resulting VCF.")
        for i, rec in enumerate(vf):
            if rec.alts and any(a.startswith("<") for a in rec.alts):
                warn(f"symbolic allele at {rec.chrom}:{rec.pos} ({','.join(rec.alts)}) - "
                     "structural calls cannot be annotated by coordinate join and need "
                     "separate handling")
                break
            if i > 1000:
                break


def clinvar_fields(path):
    with pysam.VariantFile(str(path)) as vf:
        present = [f for f in CLINVAR_FIELDS if f in vf.header.info]
    missing = set(CLINVAR_FIELDS) - set(present)
    if missing:
        warn(f"ClinVar file lacks: {', '.join(sorted(missing))}")
    if not present:
        die("ClinVar file has none of the expected INFO fields")
    return present


def check_source_fields(path, label, required):
    with pysam.VariantFile(str(path)) as vf:
        absent = [f for f in required if f not in vf.header.info]
    if absent:
        die(f"{label} has no INFO/{', INFO/'.join(absent)} field")


def source_provenance(label, path):
    """Which release of a database a classification came from, recorded at run time.

    Paths alone are not provenance: refs/clinvar.vcf.gz means nothing in a month. ClinVar
    states its own ##fileDate; for the others the size plus the file name carries the
    release.
    """
    path = Path(path)
    info = {"file": path.name, "size_bytes": path.stat().st_size}
    try:
        with pysam.VariantFile(str(path)) as vf:
            for line in str(vf.header).splitlines():
                for key in ("fileDate", "reference", "source", "dbDate"):
                    if line.startswith(f"##{key}="):
                        info[key] = line.split("=", 1)[1]
    except (OSError, ValueError):
        pass
    log(f"{label}: {info.get('fileDate', info.get('source', 'no date in header'))}")
    return info


def bcftools_version():
    if BCFTOOLS_BIN:
        out = subprocess.run([BCFTOOLS_BIN, "--version"], check=True,
                             capture_output=True, text=True).stdout
        return out.splitlines()[0] if out else "unknown"
    # pysam exposes __bcftools_version__ separately from __samtools_version__; the two
    # happen to agree today, but this field exists to make provenance certain, so ask for
    # the right one and fall back only if this pysam build does not have it.
    bcf = getattr(pysam.version, "__bcftools_version__", None) \
        or getattr(pysam.version, "__samtools_version__", "unknown")
    return (f"bcftools {bcf} / htslib {pysam.version.__htslib_version__} "
            f"(bundled in pysam {pysam.__version__})")


# --------------------------------------------------------------------------- #
# pipeline
# --------------------------------------------------------------------------- #
def run_pipeline(sample, fasta, gff, clinvar, gnomad, fields, workdir, threads, want_tsv):
    """Three annotation passes over an already-normalised sample, BCF between them."""
    csq = workdir / "step1.csq.bcf"
    # -l: local (per-record) consequences, no haplotype merging. -p a would replace the
    # 2nd variant's consequence with an @<pos> back-reference, hiding compound-het pairs
    # from any consequence filter.
    csq = cached_step("csq", csq, [sample, fasta, gff],
                      lambda o: ["csq", "--threads", threads, "-f", fasta, "-g", gff,
                                 "-l", "-Ob", "-o", o, sample])

    # No -h here: the fields keep their names, so their ##INFO definitions travel with
    # them from the ClinVar header.
    cv = cached_step("clinvar", workdir / "step2.clinvar.bcf", [csq, clinvar],
                     lambda o: ["annotate", "--threads", threads, "-a", clinvar,
                                "-c", ",".join(f"INFO/{f}" for f in fields),
                                "-Ob", "-o", o, csq])

    # gnomAD_AF:=INFO/AF is mandatory, not cosmetic: a caller-produced VCF carries its own
    # INFO/AF meaning "frequency in this sample", and writing a population frequency into
    # the same tag would destroy both. The renamed tag is declared in no header, hence -h.
    hdr = write_header_file(["gnomAD_AF"], workdir / "hdr" / "gnomad_af.hdr")
    gn = cached_step("gnomad", workdir / "step3.gnomad.bcf", [cv, gnomad, hdr],
                     lambda o: ["annotate", "--threads", threads, "-a", gnomad,
                                "-c", "INFO/gnomAD_AF:=INFO/AF", "-h", hdr,
                                "-Ob", "-o", o, cv])

    final = workdir / "annotated.vcf.gz"
    final = cached_step("final", final, [gn],
                        lambda o: ["view", "--threads", threads, "-Oz", "-o", o, gn])

    tsv = None
    if want_tsv:
        tsv = workdir / "annotated.tsv"
        fmt = ("%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%GT]\t%INFO/BCSQ\t"
               + "\t".join(f"%INFO/{f}" for f in fields) + "\t%INFO/gnomAD_AF\n")
        bcftools("query", "-H", "-f", fmt, "-o", tsv, final, out=tsv)

    return final, tsv


def verify(final, fields):
    """Count what actually landed. `bcftools annotate` exits 0 having matched nothing."""
    log("-" * 50)
    total = count_records(final)
    fill = {"BCSQ": assert_landed(final, "BCSQ", MIN_FILL["BCSQ"], total)}
    clin_key = "CLNSIG" if "CLNSIG" in fields else fields[0]
    fill[clin_key] = assert_landed(final, clin_key, MIN_FILL.get(clin_key, 0.002), total)
    fill["gnomAD_AF"] = assert_landed(final, "gnomAD_AF", MIN_FILL["gnomAD_AF"], total)
    return total, fill


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True, help="sample VCF (.vcf or .vcf.gz)")
    ap.add_argument("--fasta", required=True, help="GRCh38 reference FASTA")
    ap.add_argument("--gff", required=True, help="Ensembl GFF3 (.gff3 or .gff3.gz)")
    ap.add_argument("--clinvar", required=True, help="ClinVar VCF (.vcf.gz)")
    ap.add_argument("--gnomad", required=True, help="gnomAD AF-only VCF (.vcf.gz)")
    ap.add_argument("--workdir", default=None,
                    help="patient data: intermediates, the annotated VCF and TSV "
                         "(default <repo>/PATIENT_DATA/work)")
    ap.add_argument("--outdir", default=None,
                    help="numbers only: annotation_summary.json, no patient file names "
                         "(default <repo>/results)")
    ap.add_argument("--cache", default=None,
                    help="harmonised copies of the PUBLIC reference files "
                         "(default <workdir>/harmonised; point it at refs/ to keep "
                         "gigabytes of public data out of the patient work directory)")
    ap.add_argument("--regions", default=None,
                    help="comma-separated contigs to keep, in either naming (21 or chr21). "
                         "Two uses: a smoke test on one chromosome, which exercises the "
                         "whole chain in minutes instead of hours, and the documented "
                         "fallback when the genome-wide run does not fit the time budget - "
                         "if it ships, the methods section has to say so")
    ap.add_argument("--skip-source-norm", action="store_true",
                    help="do not normalise ClinVar and gnomAD. Saves roughly half an hour "
                         "once, at the cost of silently missing indels written "
                         "differently on the two sides - mostly at multiallelic sites. "
                         "If you use it, say so in the methods section")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--no-tsv", action="store_true",
                    help="skip the flat table (it is multi-GB at WGS scale)")
    args = ap.parse_args()

    workdir = Path(args.workdir) if args.workdir else ROOT / "PATIENT_DATA" / "work"
    outdir = Path(args.outdir) if args.outdir else ROOT / "results"
    cache = Path(args.cache) if args.cache else workdir / "harmonised"
    for directory in (workdir, outdir, cache):
        directory.mkdir(parents=True, exist_ok=True)
    mapfile = cache / "chr_to_numeric.txt"
    write_map(mapfile)

    log(f"bcftools: {bcftools_version()}")
    log(f"using bcftools from: {BCFTOOLS_BIN or 'pysam (bundled)'}")

    wanted = None
    if args.regions:
        wanted = {to_numeric(c.strip()) for c in args.regions.split(",") if c.strip()}
        log(f"restricted to {len(wanted)} contig(s): {', '.join(sorted(wanted))}")

    fasta = prepare_fasta(args.fasta, cache)
    gff = prepare_gff(args.gff, cache)
    check_not_gvcf(args.vcf)

    sample, sample_stats = harmonise_vcf(args.vcf, "sample", fasta, workdir, mapfile,
                                         args.threads, wanted=wanted)
    do_norm = not args.skip_source_norm
    if not do_norm:
        warn("--skip-source-norm: ClinVar and gnomAD are left unnormalised, so indels "
             "written differently on the two sides will not match. This has to go in the "
             "methods section.")
    clinvar, _ = harmonise_vcf(args.clinvar, "clinvar", fasta, cache, mapfile,
                               args.threads, do_norm=do_norm)
    gnomad, _ = harmonise_vcf(args.gnomad, "gnomad", fasta, cache, mapfile,
                              args.threads, do_norm=do_norm)

    check_build(fasta, gff, [("sample", sample), ("clinvar", clinvar), ("gnomad", gnomad)])
    fields = clinvar_fields(clinvar)
    check_source_fields(gnomad, "gnomAD", ["AF"])
    provenance = {"clinvar": source_provenance("clinvar", args.clinvar),
                  "gnomad": source_provenance("gnomad", args.gnomad),
                  "gff3": {"file": Path(args.gff).name,
                           "size_bytes": Path(args.gff).stat().st_size},
                  "fasta": {"file": Path(args.fasta).name,
                            "size_bytes": Path(args.fasta).stat().st_size}}

    final, tsv = run_pipeline(sample, fasta, gff, clinvar, gnomad, fields,
                              workdir, args.threads, want_tsv=not args.no_tsv)
    total, fill = verify(final, fields)

    # No patient file names here: this file is meant to be safe to commit, and a sample
    # file name normally carries the sample or patient identifier.
    summary = {
        "tools": bcftools_version(),
        "dispatch": BCFTOOLS_BIN or "pysam (bundled)",
        "regions": args.regions or "genome-wide",
        "source_normalisation": "skipped" if args.skip_source_norm else "norm -m -any -f",
        "n_records": total,
        "fill": fill,
        "sample": sample_stats,
        "clinvar_fields": fields,
        "reference_provenance": provenance,
        "caveats": [
            "gnomAD af-only carries GLOBAL AF only: no popmax/grpmax and no nhomalt, so "
            "BA1/BS1/PM2 cannot be evaluated from it alone.",
            "A missing gnomAD_AF means absent from THIS resource, not AF=0 - the resource "
            "is a liftover of an older gnomAD release and its non-coding coverage is "
            "thinner than the current genomes.",
            "bcftools csq emits no HGVS and has no NMD-escape logic; both are done by hand "
            "for the final candidates.",
        ],
    }
    summary_path = outdir / "annotation_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))

    log(f"done:    {final}")
    if tsv:
        log(f"table:   {tsv}")
    log(f"summary: {summary_path}")
    log("Intermediates in the work directory can be deleted once the final file verifies.")


if __name__ == "__main__":
    main()