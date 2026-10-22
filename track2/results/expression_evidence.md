# Expression evidence from public cell and animal data

Question asked: is there cell or animal expression data for the shortlisted drugs that
would support or undermine the mechanism, by showing what happens to SAC-related pathways?

Answer in one line: **the one usable dataset argues against the leading candidate at the
exposure tested, and the reason it argues against it is instructive.**

## The caveat that governs the interpretation

**Spindle-assembly-checkpoint function is post-translational.** It runs on KARD
phosphorylation by PLK1, APC/C-mediated degradation, and assembly of the mitotic checkpoint
complex. mRNA levels of SAC genes are not a measure of checkpoint competence: sodium
phenylbutyrate could raise BUBR1 protein without moving a single transcript, and the
chaperone hypothesis is invisible to RNA-seq by construction.

What expression data *can* establish is **target engagement and the direction of pathway
perturbation**. That is the question answered below, and it is not the same question as
"does the checkpoint work again".

A second confounder dominates everything here: **SAC and kinetochore genes are a
proliferation program.** Anything that stops cells dividing lowers all of them. A drop in
the SAC module is therefore uninformative unless proliferation is controlled for.

## GSE134781 - the right model, the wrong tissue

Sieben 2020 (PMID:31738183), the mouse allelic series behind this project's mechanism
section. Design is exactly what the question needs: **3 genotypes x 3 replicates** -
wild-type, `hetL1002P` (the missense allele, equivalent to human p.Leu1012Pro) and
`hetX753` (the truncating allele) - in skeletal muscle at 3 months. Skeletal muscle is
phenotypically relevant: the proband has muscle atrophy and Sieben reported sarcopenia.

**It cannot answer the question.** Skeletal muscle is post-mitotic and the mitotic module is
not expressed in it:

- `Mki67` mean CPM **0.32** - the proliferation marker is effectively absent
- median CPM across the 62 mapped SAC nodes: **0.63**
- only **28 of 62** nodes clear a 1-CPM-in-3-samples filter
- raw counts in the core: `Bub1` 0-1, `Ttk` 0-2, `Aurkb` 0-4, `Plk1` 1-5, `Birc5` 1-6,
  `Trip13` 0-5, `Cdc20` 1-11, on libraries of 20-34 M reads

A competitive test of the SAC set against background found no shift in any contrast
(L1002P vs WT p=0.22; X753 vs WT p=0.37; L1002P vs X753 p=0.90). That is a null result from
an underpowered assay, not evidence of no effect.

**A result retracted during this analysis.** `Bub1b` initially appeared reduced ~2-fold in
both mutant genotypes (X753 0.46x, p=0.002; L1002P 0.45x, p=0.033), which was tempting
because NMD predicts exactly that for the truncating allele. The raw counts are **21, 11, 16
(WT) vs 10, 9, 2 (L1002P) vs 8, 8, 5 (X753)**. At that depth the ratio is noise, and the
equal reduction in the *missense* genotype - which should not trigger NMD - confirms it is
not an NMD signal. The claim is withdrawn. It is recorded here because it is the kind of
result that survives into a report if nobody looks at the counts.

*Aside, not pursued:* `Bub3`, `Mad1l1` and `Aurka` are robustly expressed in muscle while
the rest of the module is not, consistent with non-mitotic moonlighting roles.

## GSE86928 - sodium phenylbutyrate, and the answer is a caution

A549 cells treated with 5 mM sodium butyrate or 5 mM **sodium phenylbutyrate** against
untreated control (PMID:28033672). Proliferating human cells, so the mitotic module *is*
expressed. The submitters frame phenylbutyrate as an **HDAC inhibitor**, which is a
different proposed mechanism from the chemical-chaperone one in `bebok_track2_report.md`.

**What the data show.** Of the SAC and kinetochore nodes present in the submitters' filtered
probe list, **32 of 33 move down**, including:

| gene | fold change | | gene | fold change |
|---|---|---|---|---|
| KNTC1 | -4.03 | | TTK | -2.54 |
| NUF2 | -3.97 | | BUB1 | -2.50 |
| MAD2L1 | -3.85 | | **BUB1B** | **-2.40** |
| NDC80 | -3.30 | | BIRC5 | -2.25 |
| TRIP13 | -3.14 | | AURKB | -2.21 |
| CDC20 | -3.01 | | **USP44** | **+3.39** (the only one up) |

**Why the drop is not a finding about the checkpoint.** The proliferation program collapses
alongside it: `MKI67` -3.65/-3.56/-2.20/-2.01, `CCNA2` -3.29/-4.58, `FOXM1` -2.99,
`TOP2A` -2.22, `CCNB1` -2.18, with `HIST1H1C` +2.08 consistent with HDAC inhibition. The
cells have arrested, and every mitotic gene follows. The SAC module did not move *as a
module*; it moved as part of the cell cycle.

**Why it is nonetheless a caution against the candidate.** Our own signed network
(`data/signed_network.tsv`) gives BUB1B a required direction of **increase**. At this
exposure sodium phenylbutyrate produces a **2.4-fold decrease in BUB1B transcript**. Even
granting that the chaperone effect is post-translational and invisible here, a proteostatic
gain would have to overcome a transcriptional loss in any proliferating tissue. That is the
wrong direction, and it belongs in the report rather than in a footnote.

**Limitations, which are severe.** n=1 per condition with no replicates, so no statistics are
possible and only fold changes are interpretable. 5 mM is supra-therapeutic relative to
plasma levels in urea-cycle treatment. A549 is a lung adenocarcinoma line, not a patient
cell. The supplementary file is a **filtered list of changed probes**, so a gene's absence
means "not in the filtered list", not "unchanged" - `CDKN1A`, `PLK1`, `E2F1` and `MYC` are
absent for that reason, and a proper competitive gene-set test would require renormalising
the raw CEL files.

## Datasets found for the other candidates

| candidate | GEO survey result |
|---|---|
| hydroxychloroquine | 56 series, mostly autophagy-inhibition combinations in oncology (e.g. GSE275078 n=136 human, GSE274516 PDAC GEMM); none in a chromosomal-instability context |
| ataluren | 3 series, none a readthrough experiment in a relevant model; GSE307222 (aniridia cornea organoids) is the nearest analogue |
| acetylcysteine | 557 series, overwhelmingly as an antioxidant co-treatment rather than the perturbation under test |
| metformin, dasatinib | abundant, none addressing aneuploidy tolerance |

Related BubR1 datasets worth reading but not drug experiments: GSE277997 (BubR1
insufficiency, cardiac, n=12), GSE270226 (BubR1-knockout embryonic heart, single cell),
GSE338948 (human BUB1B dependency, identifies CENP-E inhibition as the tractable approach).

## What would actually test the hypothesis

Expression profiling is the wrong assay for the leading candidate, and this exercise shows
why. The hypothesis is that a chemical chaperone raises the steady-state **protein** level of
a destabilised missense product past Matsuura's >50%-reduction threshold. The experiment is
a **BUBR1 western blot in a patient-derived or p.Asn1002Lys knock-in line, plus or minus
drug**, with a PCS rate as the functional secondary readout - as set out in
`validation_plan.md`.

If transcriptomics is used at all, it must be in a **proliferating** cell at a
**therapeutically achievable** concentration, with proliferation explicitly controlled - for
example by normalising the SAC module against a general cell-cycle signature, so that a
specific effect can be separated from arrest. Neither public dataset meets that
specification.
