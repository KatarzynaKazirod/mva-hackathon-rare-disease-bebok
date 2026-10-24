# rare-disease-bebok

## About

Submission for **"Rare Disease, Real Kid: The MVA Hackathon 2026"** (Sage Bionetworks with the
MVA Society, Hugging Face and BEACON), for a child with suspected mosaic variegated aneuploidy.
Two tracks, both in this repository:

| track | question | where |
| --- | --- | --- |
| **Track 1** | which variant explains the phenotype? | this repository root |
| **Track 2** | does an existing medicine act on that biology? | [`track2/`](track2/) |

Inputs are a single-sample GRCh38 WGS VCF and seven HPO terms - no BAM/CRAM, no parental data.

## Environment setup

Requires [uv](https://docs.astral.sh/uv/), Python 3.12 and `bcftools`.

```bash
uv sync                           # dependencies
brew install bcftools             # macOS; required by the annotation step
```

Track 2 needs none of this: its scripts are standard-library only and run under plain
`python3`.

## Track 1

### Where things live

```
PATIENT_DATA/        input VCF and all intermediates            NOT tracked
  raw/                 put the delivered .vcf.gz + .tbi HERE
  work/                annotated.vcf.gz, annotated.tsv, intermediates
refs/                reference files from download_refs.py      NOT tracked, ~14 GB
config/              hand-edited inputs
  patient_hpo.example  TEMPLATE: format, the .docx it comes     tracked
    .txt               from, and why relatives' terms are out
  patient_hpo.txt      HPO terms read from the challenge's      NOT tracked
                       Challenge_Clinical_Phenotype_1.docx;
                       7 active + 1 commented out (maternal)
  family_history.txt   relatives' HPO terms - a note to the     NOT tracked
                       reader; no script reads this file
  mva_genes_extra.tsv  genes HPO cannot supply, one PMID each
results/             generated tables                           NOT tracked, two exceptions
src/                 the pipeline, one executable per step
submission/track_1/  what is actually submitted
  bebok_hpo_ordered_   the ranked answer, organiser's column layout
    cascade.csv
  METHODS.md           the methods-description form
  bebok_track1_        the full write-up
    report.md
track2/              drug repurposing - self-contained, see track2/README.md
```

### The data that matters

| file | what it is |
| --- | --- |
| **`submission/track_1/bebok_hpo_ordered_cascade.csv`** | **the answer.** Row 1 is the primary claim; the `epcr` column *is* the ranking |
| `results/final_candidates.tsv` | the same table as the pipeline writes it, before CSV conversion |
| `results/finalists_evidence.tsv` | one row per allele: frequency, NMD status, ClinVar |
| `results/ranking.json` | the tier rules behind the order, plus the provenance of the run |
| `results/cascade_counts.json` | how many variants survived each filter |

Everything else in `results/` exists to explain how row 1 got there.

## How to run the pipeline

Run from the repository root, in this order.

### Step A - get the patient data

The dataset is **gated**: request access on the Hugging Face dataset page, accept the
Hackathon rules, and wait for approval before anything below will work. Then download **two
files, from the same source and in the same pass**:

```bash
mkdir -p PATIENT_DATA/raw
# from the gated HF dataset page, into PATIENT_DATA/raw/ :
#   <sample>.vcf.gz        the variant calls
#   <sample>.vcf.gz.tbi    the tabix index - NOT optional
```

The index matters more than it looks. A `.tbi` holds byte offsets and **no link back to its
data file** - no filename, no checksum - so htslib simply trusts that `X.vcf.gz.tbi` describes
`X.vcf.gz`. A mismatched or stale index does not error; region queries silently return wrong or
empty results. Download the pair together and do not substitute one from another copy.

BAM/CRAM are offered too and are **not needed here** - this pipeline reads only the VCF. The
full dataset is ~85 GB; the VCF and its index are ~1.7 GB.

Two rules attach to this file from the moment you have it: it **may not be reshared through any
channel**, and it must be **deleted within 30 days of the Hackathon close** - from every
machine and private repository - with confirmation emailed to the organisers.

### Step B onward - the pipeline

Point `VCF` at the downloaded file:

```bash
export VCF=PATIENT_DATA/raw/<your-file>.vcf.gz

# 0. references - once per clone, a few hours, needs 20 GB free
./src/download_refs.py

# 1. what IS this file?                                            minutes
./src/audit_vcf.py --vcf "$VCF"

# 2. consequence + ClinVar + gnomAD                                hours, genome-wide
./src/annotate_vcf.py --vcf "$VCF" \
    --fasta refs/GRCh38.fa \
    --gff   refs/Homo_sapiens.GRCh38.116.chr.gff3.gz \
    --clinvar refs/clinvar.vcf.gz \
    --gnomad  refs/af-only-gnomad.hg38.ensembl.vcf.gz \
    --cache refs/

# 3. which GENE fits the phenotype - no variants involved          seconds
#    the template explains the format and where the terms come from
cp config/patient_hpo.example.txt config/patient_hpo.txt
$EDITOR config/patient_hpo.txt   # transcribe Challenge_Clinical_Phenotype_1.docx
./src/build_panel.py --leave-one-out

# 4. which VARIANTS qualify - genome-wide, panel not consulted     minutes
./src/cascade.py

# 5. which CLAIM is most likely                                    seconds
./src/rank_candidates.py --top 5 --proband-id <id-from-the-dataset>
```

Three things are worth knowing before you run it:

**The panel and the cascade are logically independent** - the panel needs no variants, the
cascade needs no panel - and that independence is demonstrated rather than asserted. Run
`./src/rank_candidates.py --panel-mode agnostic` and the phenotype is discarded entirely: the
candidate set is unchanged and *BUB1B* still comes first. An earlier version of this pipeline
instead committed a hash of the term list before the cascade, as a pre-registration record.
That was removed - the challenge does not ask for it, and the agnostic run is better evidence
than a commitment anyway, because it tests the claim instead of attesting to it.

**Everything is harmonised to Ensembl chromosome naming** (`1`, not `chr1`). That is the single
invariant the rest of the pipeline relies on. The submission CSV is the one place it is
converted back, because the organisers' template is `chr`-prefixed.

**`src/` is a directory of executables, not a package.** Each anchors its paths to the
repository root via `Path(__file__).resolve().parents[1]`, so moving a script elsewhere
silently relocates its outputs - it will run without error and write to the wrong place.

Steps 0 and 1 are independent of everything else, so the audit can run while the references
download.

## Track 2

Drug repurposing, taking the Track 1 mechanism as its input: *does an existing medicine act on
this biology?* Self-contained in [`track2/`](track2/) with its own `config/`, `data/`,
`results/` and `src/`; standard-library only, no API key, one command to reproduce:

```bash
cd track2 && python3 src/run_pipeline.py --mode check
```

**The central result is negative, and that is the finding.** Of 65 spindle-checkpoint and
kinetochore nodes, 60 carry no drug of any kind and none of the rest has an approved one; the
17 cohesion-layer targets return zero approved drugs, separase included. Asked of the disease
instead of the gene, the nine MVA entries return no medicine at all. Every candidate that
survives therefore acts **downstream** of the lesion - 14 of them, in four arms, each with a
ChEMBL identifier, a verified mechanism citation and product-level approval checked against
openFDA, alongside an explicit list of what must **not** be given to this child.

Everything is in one document:
[`track2/bebok_track2_report.md`](track2/bebok_track2_report.md) - mechanism, the negative
result, the candidate arms, the contraindications, how it was produced, and the organisers'
methods form as Appendix A. [`track2/README.md`](track2/README.md) is the directory map.

## Data privacy

Everything derived from the VCF or the phenotype is patient material. The question for each
file is not whether it carries an identifier, but whether the patient can be inferred from it.

| not committed | why |
| --- | --- |
| `PATIENT_DATA/**` | the genotypes themselves |
| `results/panel.meta.json` | contains the HPO terms outright |
| `results/panel.tsv` | a ranking of ~1400 genes is a function of those terms, so the phenotype can be substantially inferred |
| `config/patient_hpo.txt`, `config/family_history.txt` | the phenotype, and a relative's |

### The phenotype files

The proband's clinical features arrive as `Challenge_Clinical_Phenotype_1.docx` in the gated
dataset and are transcribed by hand into `config/patient_hpo.txt`, one HPO identifier per line.
A relative's feature goes into `config/family_history.txt` instead, which no script reads: it
exists so that a term left out of the scoring is visibly set aside rather than apparently
forgotten. `config/patient_hpo.example.txt` is the tracked template that documents the
format; the filled-in files are not tracked and must not become tracked.
