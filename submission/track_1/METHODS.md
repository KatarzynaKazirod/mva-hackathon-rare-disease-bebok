# Methods description - MVA Hackathon 2026, Track 1 (Variant Prediction)

Answers to the organisers' methods-description form. One form per model/approach; this is
the only approach submitted. Question wording is the organisers'; the answers are ours.

Approach name: **`bcftools-hpo-cascade`**

---

### Team name

**Bebok**

### Model number (fill out this form once per model/approach submitted; up to 6)

1 of 1.

---

### (required) If commercially-available Generative AI was used, record the provider, the plan or tier, and the relevant setting in your methods description.

Anthropic Claude, Team plan subscription - commercial terms, customer content is not used to
train models - model `claude-opus-5`, used through Claude Science. One model did all of the
work: no secondary or utility model contributed and no sub-agents were delegated.

Generative AI was used for pipeline code, literature retrieval and drafting. **Variant
classification was not delegated to it:** every ACMG criterion reported here was assigned
outside the code, and no part of the pipeline writes, scores or validates a classification.
The ranking does not read one either; the ACMG points and the `epcr` values are independent
numbers answering different questions (see Declarations in the report).

---

### Please describe your model/approach in detail.

A six-step command-line pipeline over a single-sample GRCh38 VCF. Each step is one
self-contained script that writes a machine-readable provenance file beside its output, so
this description can be checked against the run rather than taken on trust.

**Step 0 - references** (`download_refs.py`). GRCh38 FASTA, Ensembl GFF3 release 116, ClinVar,
gnomAD, HPO and MANE Select, pinned by release.

**Step 1 - audit** (`audit_vcf.py`). Establishes what the file actually is before anything is
inferred from it: genome build, contig naming, sample count, FILTER values, FORMAT fields
present. Contig naming is normalised to Ensembl (`15`, not `chr15`) once, here, and stays that
way; the `chr` prefix the submission template wants is added at the single point where the
submission CSV is written.

**Step 2 - annotation** (`annotate_vcf.py`). Consequences from `bcftools csq` against the
Ensembl GFF3, plus ClinVar significance and review status, plus gnomAD allele frequency. The
annotation step reports its own coverage: 0 reference-allele mismatches, BCSQ present on
68.2% of records, ClinVar on 0.89%, gnomAD frequency on 92.8%.

**Step 3 - phenotype panel** (`build_panel.py`). The clinical features come from
`Challenge_Clinical_Phenotype_1.docx`, the narrative description supplied with the gated
dataset. They were read from it and transcribed by hand to HPO identifiers in
`config/patient_hpo.txt`, which is the only route by which the phenotype enters this pipeline.
That file carries eight HPO identifiers but only **seven** active ones. The eighth, HP:0200067
*Recurrent spontaneous abortion*, is the **mother's** feature rather than the child's, and is
**commented out**, so `build_panel.py` never reads it as a proband term. It is instead recorded
in `config/family_history.txt`, which **no script reads**. The two files are independent and do
different jobs: the comment in `patient_hpo.txt` is what keeps the term out of the scoring, and
the family-history file is a note to the reader, so that a term left out is visibly set aside
rather than apparently forgotten. Nothing about the maternal history reaches a score, a rank or
a provenance file.
The proband's features are scored
against HPO disease annotations using Resnik similarity over the most informative common
ancestor, best-match average, patient-to-disease direction. Output is a ranking of 1,449
genes. **The panel orders candidates; it cannot remove one.**

**Step 4 - inheritance cascade** (`cascade.py`). Genome-wide, with no gene list applied:

| stage | records |
| --- | --- |
| all records | 4,962,060 |
| `FILTER = PASS` | 4,717,955 |
| quality gate (GQ >= 20, DP >= 10, allele balance) | 4,445,604 |
| population frequency | 170,569 |
| protein-altering consequence on a protein-coding transcript | 590 |

Compound-heterozygous candidates are capped at AF <= 0.001 and homozygous at <= 0.005, with a
ClinVar P/LP assertion admitted up to 0.02. Heterozygous calls require allele balance between
0.25 and 0.75, homozygous calls >= 0.85. chrY and chrM are excluded; sex is not inferred, so
chrX rows are flagged rather than filtered. A qualifying consequence must fall on a
protein-coding transcript - the single exception being a ClinVar P/LP assertion, counted
separately. This left 118 genes with recessive candidates and 5 with a dominant single-allele
claim.

**Step 5 - claims and ranking** (`rank_candidates.py`). Per-variant rows are collapsed into
*claims*: a compound-heterozygous pair is one claim about one gene, not two guesses. Each
claim is annotated with gnomAD v4.1 detail (nhomalt, grpmax, faf95) and an NMD-escape check
read from the GFF3 for the MANE Select transcript, then assigned an evidence tier and ordered.
The tier is set by `tier_of()` from four recorded facts and nothing else - checklist membership
or a panel rank in the top ten, a ClinVar P/LP assertion carrying review criteria, the panel
rank, and the inheritance arm - and maps to a fixed `epcr`: A 0.920, B 0.500, C 0.150, D 0.040,
E 0.020, with secondary findings pinned at 0.010 beneath every primary claim. The values are
written in the source rather than fitted to anything, and the function is deterministic, so the
same claim always yields the same tier. It is not derived from the ACMG framework; the report's
Declarations set out what that costs.

**The MVA gene list is never a filter.** It is a flag column plus a safety net: checklist-gene
variants that failed only on frequency or quality are written to `set_aside_checklist.tsv`
(27 variants) so a near-miss in a known gene is visible rather than lost. This matters because
the published MVA gene list keeps growing, and a panel cannot find a gene that is not on it.

**How much the answer depends on the phenotype was measured, not asserted.**
`rank_candidates.py --panel-mode` re-runs the ordering with less prior knowledge. The
candidate set is identical in all three runs - 123 claims - because the cascade finishes
before the panel is read:

| mode | knowledge used | *BUB1B* tier | epcr | rank |
| --- | --- | --- | --- | --- |
| `full` | phenotype panel + MVA gene list | A | 0.920 | 1 |
| `no-panel` | MVA gene list only | A | 0.920 | 1 |
| `agnostic` | neither | E | 0.020 | 1 |

The gene ranks first in all three. The phenotype is therefore confirmatory rather than
selective - but it is what makes the claim confident: in agnostic mode the separation
collapses to 0.020 against 0.019 for the next candidate.

### Please state whether the submission file is the automated output of your computational approach (preferred), or if it has undergone downstream manual review and curation.

**Automated output.** The submitted rows, their order and the `epcr` values are written by
`rank_candidates.py` with no manual editing. Two things are supplied by hand and neither
changes which variants are reported or in what order: the `proband_id` from the dataset, and
the `chr` prefix added where the CSV is written.

Manual work exists in this project but is deliberately **outside** the submission file: the
ACMG classification and the mechanistic interpretation are written by hand, live in the
report, and are not inputs to the ranking.

### Please describe any downstream manual review in detail.

The automated output was read and four defects were found and fixed *in the code*, after which
the pipeline was re-run - no row was hand-corrected.

1. **Transcript biotype was not required.** `bcftools csq` reports per transcript, and
   consequences were being accepted on NMD, retained-intron and non-stop-decay transcripts.
   This produced a spurious *RPS17* "compound heterozygous pair of missense variants" in which
   one allele was in fact synonymous and neither touched MANE Select. Requiring a
   protein-coding transcript removed *RPS17* and *KMT2C*, cut 749 rows to 590 and dropped 164
   gene-level claims.
2. **The NMD check inverted exon numbering on the minus strand**, sorting exons by coordinate
   rather than transcript order. Caught on *FLG*.
3. ***BUB1B* appeared twice** - once as a dominant single-allele claim, once as the pair - and
   the single allele sorted above the pair, putting half the answer at the top. Claims are now
   collapsed per gene, keeping the pair.
4. **Absence of evidence was reported as evidence of absence.** The gnomAD v4.1 extracts cover
   the checklist genes only; for other genes the notes now say the frequency was not looked up.

The reported variants were also verified against MANE Select numbering by four independent
routes - ClinVar, GFF3 CDS arithmetic, Ensembl VEP and the Open Targets variant record -
because `bcftools csq` reported them on a non-MANE transcript (`751L>751*`, `1016N>1016K`).

### Please state whether your approach only used publicly available data (preferred) or if proprietary data was also used.

**Publicly available data only.** No proprietary database, no institutional variant archive,
no in-house cohort.

### Please describe any public data used in detail (e.g., reference panels, population databases, literature).

| source | release / version | used for |
| --- | --- | --- |
| GRCh38 reference FASTA | GRCh38 | reference-allele check, left-alignment |
| Ensembl GFF3 | release 116 | consequences, CDS arithmetic, NMD-escape check |
| MANE Select | current at run time | the transcript every reported change is expressed on |
| ClinVar | release current at run time | significance and review status (stars) |
| gnomAD v4.1 | exomes and genomes | allele frequency, nhomalt, grpmax, faf95 |
| HPO | `hp.obo` 2026-09-01, `phenotype.hpoa` 2026-09-02 | phenotype-to-gene ranking |
| OMIM / Orphanet | via HPO disease annotations | disease-gene mapping, incl. the grouped MVA entry |
| Open Targets Platform | 26.09 | LOFTEE, GERP and FoldX values quoted as observations |
| PubMed | - | every citation resolved by PMID before use |

Clinical interpretation rests on published work, each citation resolved rather than recalled -
including the BUBR1 dosage threshold (PMID:16411201), truncating-plus-missense *BUB1B* MVA
families with embryonal rhabdomyosarcoma (PMID:15475955), the mouse allelic series
(PMID:31738183), ClinGen's predictor calibration (PMID:36413997), the ACMG/AMP framework
(PMID:25741868) with its point-based implementation (PMID:32720330), the ClinGen SVI PVS1
specification (PMID:30192042) and the SVI recommendation against the reputable-source criteria
PP5 and BP6 (PMID:29543229), which bears on what the ranking rests on and is discussed in the
report's Declarations.

### Please describe any proprietary data used in detail.

None.

### Is your approach able to output proposed pairs of compound heterozygous candidate variants, or only single candidate variants?

**Pairs.** This is the primary output mode, not an add-on: the recessive arm emits a
compound-heterozygous pair as a single claim populating both allele columns of the submission
layout, and the submitted answer is such a pair. The dominant arm emits single-allele claims
into the first allele columns only.

Two caveats stated in the output rather than hidden. **Phase is not proven** - with one sample
and no reads, two heterozygous calls in one gene are a hypothesis, not two alleles shown to be
in trans; `pair_distance_bp` records whether read-backed phasing from a CRAM slice could
resolve it (here 10,911 bp - too far for short reads). And because the arms are independent,
one variant can appear in both, so rows are keyed on `(arm, chrom, pos, ref, alt)` and never
on position alone.

### How did you handle secondary or incidental findings, if any? (See the finding_type column in your submission.)

A single-allele claim is reported as a secondary finding only if it is **truncating AND
NMD-competent AND carries a ClinVar P/LP assertion with at least one review star**. Four
candidates were tested and one survived: *LZTR1*. *RBM8A* was rejected (5' UTR), *GNRHR*
(missense) and *FLG* (escapes NMD, last exon 3 of 3). Every rejection is logged with its
reason in `results/ranking.json`.

Secondary findings are placed below every primary candidate, which is what `epcr` asks for -
by definition they are unlikely to explain the presenting phenotype.

**The label is a pipeline default, not a clinical recommendation, and the output says so.**
The ACMG SF v3.3 list (Genet Med 2025;27(8):101454) was read rather than assumed: **none** of
the four tested genes is on it. A principled test would match each gene's inheritance mechanism
to the observed zygosity, which needs per-gene inheritance annotation the pipeline does not
carry; the truncating-plus-ClinVar rule stands in for it and `ranking.json` records that.

### Please provide an estimate of run time and cost for your approach (if possible).

Commodity laptop, no cloud compute, no paid data source - so the monetary cost is zero beyond
the machine and the operator.

| step | wall time |
| --- | --- |
| 0 - references (once per clone, ~7 GB) | hours |
| 1 - audit | minutes |
| 2 - annotation | hours |
| 3 - phenotype panel | seconds |
| 4 - cascade | minutes |
| 5 - claims and ranking | seconds |

Per additional sample, only steps 1-5 repeat: **a few hours, dominated entirely by
annotation.** The analyst time behind the manual classification and interpretation is the real
cost, and it is larger than the compute.

### Provide a method abstract (up to 500 words). Please include strengths and limitations of your method as well as methodological detail.

We identify causal variants with a genome-wide inheritance cascade, using disease knowledge
only to order the result. A single-sample GRCh38 VCF is audited for build and contig naming,
annotated with `bcftools csq` against Ensembl 116 plus ClinVar and gnomAD, then filtered by
inheritance model: 4,962,060 records to 590 protein-altering calls on protein-coding
transcripts, yielding 118 genes with recessive candidates. Separately, the proband's clinical
features are scored against HPO disease annotations by Resnik similarity to rank 1,449 genes.
The panel orders candidates and cannot remove one; the MVA gene list is a flag column and a
safety net, never a filter. Per-variant rows are collapsed into claims, so a compound
heterozygous pair is one claim about one gene.

The submitted answer is a compound heterozygous pair in *BUB1B* (`NM_001211.6`):
`c.2210T>G` p.Leu737Ter, a ClinVar Pathogenic/Likely pathogenic nonsense variant 8,861 bp
upstream of the last exon junction and therefore NMD-competent, classified Pathogenic at 10
ACMG points; and `c.3006T>G` p.Asn1002Lys, a gnomAD singleton in the pseudokinase domain,
classified a variant of uncertain significance at 2 points.

Strengths. The search is genome-wide, so a novel disease gene remains reachable - a panel
cannot find a gene that is not on it, and the MVA gene list is still growing. The phenotype
panel is built and cryptographically stamped before the cascade runs, which makes "the
ranking was not tuned against the answer" checkable rather than asserted; both hashes
reproduce bit-identically. No missense pathogenicity predictor is used anywhere, and PP3 is
reported as not assessed, because ClinGen calibration requires a predictor chosen in advance
with published thresholds. Dependence on the phenotype is measured, not argued: with the panel
and the gene list both switched off, *BUB1B* still ranks first on variant evidence alone.

Limitations. Phase is unproven - two heterozygous calls in one gene are a hypothesis, and
10,911 bp is beyond short-read phasing. No parental data, so de novo status cannot be
established and the dominant arm depends entirely on prior ClinVar classification. No BAM or
CRAM, so there is no coverage track and no structural or copy-number analysis: a deletion in
trans with a point variant would be invisible, and a balanced parental rearrangement cannot be
excluded - relevant here, because the mother has a history of recurrent pregnancy loss. The
second allele's effect is inferred from protein context and a destabilisation estimate, not
measured, so it remains of uncertain significance; a BUBR1 western blot is the one experiment
that would resolve it. ACMG classification is manual. The `epcr` values are ordering weights
derived from an evidence tier, not calibrated probabilities. Nephrocalcinosis, one of the
seven clinical features, is plausibly secondary to prematurity rather than to the genotype.

---
