# Track 1 - variant identification

Suspected mosaic variegated aneuploidy in a child. Single-sample GRCh38 whole-genome VCF,
seven clinical features coded in HPO, no BAM or CRAM, no parental samples.

**Answer submitted: a compound heterozygous pair in *BUB1B* (MANE Select `NM_001211.6`),
the gene for mosaic variegated aneuploidy syndrome 1.**

| | allele 1 | allele 2 |
| --- | --- | --- |
| GRCh38 | 15:40209701 T>G | 15:40220612 T>G |
| coding | `c.2210T>G` | `c.3006T>G` |
| protein | p.Leu737Ter | p.Asn1002Lys |
| consequence | stop_gained, exon 17/23 | missense, exon 23/23 |
| ClinVar | Pathogenic/Likely pathogenic, 2 stars | not present |
| gnomAD v4.1 exomes | 115/1,461,846 (no homozygotes) | 1/1,461,878 - singleton |
| ACMG | PVS1 + PM2_supporting + PM3_supporting = 10 points, **Pathogenic** | PM2_supporting + PM3_supporting = 2 points, **Uncertain significance** |

Numbering on the MANE Select transcript was confirmed four independent ways before anything was
reported on it: ClinVar, CDS arithmetic from the Ensembl GFF3, an Ensembl REST VEP query, and
the Open Targets variant record. This matters because `bcftools csq` reports per transcript and
picked a non-MANE one, where the same two alleles read `751L>751*` and `1016N>1016K` - those
numbers appear in `results/candidates.tsv` and are **not** the ones to quote.

The alleles are 10,911 bp apart. Phase is **not** proven: with one sample and no reads, two
heterozygous calls in one gene are a compound-heterozygous *hypothesis*. What raises it above
a coincidence is that the gene is the one the phenotype points at and one allele is an
independently asserted null.

---

## 1. What this case needed, and why a gene panel was not the way to get it

MVA is caused by biallelic variants in any of about eleven genes, and the published list keeps
growing. A panel of those eleven genes would have found this answer - and would have been the
wrong instrument, because it cannot find anything else. A child with this presentation and a
variant in a twelfth gene returns a negative panel, and a negative panel is read as "not MVA"
rather than "not on our list".

So the design decision underneath everything here is: **the search is genome-wide and the
disease knowledge is only ever used to order the output.** The eleven-gene list exists in the
repository, but it is a flag column and a safety net, never a filter. The phenotype is used
the same way. Section 6 measures how much the answer actually depended on either.

### Four constraints fixed before any data was read

1. **Patient data never leaves the machine.** Everything under `results/` is git-ignored, and
   the two config files holding the phenotype are checked against the git index on every
   `build_panel.py` run, which refuses to proceed if either is tracked.
2. **The phenotype panel orders, it never filters.** It cannot remove a candidate.
3. **No missense pathogenicity predictor anywhere** - not in filtering, not in
   classification. PP3 is reported as *not assessed*. ClinGen's calibration work
   (PMID:36413997) requires a predictor chosen in advance with published thresholds; picking
   one after seeing the variant is how a predictor becomes a rationalisation.
4. **Every stage records its own limitations**, so the methods description can be assembled
   from provenance files rather than from memory.

---

## 2. The clinical picture and what it implies about inheritance

Seven features were provided: nephrocalcinosis, failure to thrive, small for gestational age,
premature birth, **rhabdomyosarcoma**, skeletal muscle atrophy and short stature. The mother
has a history of recurrent pregnancy loss.

All of this was read from `Challenge_Clinical_Phenotype_1.docx`, the clinical description
supplied with the dataset, and transcribed by hand into HPO identifiers. The transcription
carries eight HPO identifiers but only seven active ones: the maternal recurrent pregnancy
loss (HP:0200067) is present as a commented-out line, so the panel never scores it, and is
recorded separately in `config/family_history.txt` for audit. The reason it is held out is set
out below.

Two of these carry most of the information. **Rhabdomyosarcoma in early childhood** is a
cancer-predisposition signal, and embryonal rhabdomyosarcoma is specifically reported in MVA
with biallelic *BUB1B* (PMID:15475955). The combination of intrauterine growth restriction,
prematurity and postnatal growth failure is the MVA growth phenotype.

Nephrocalcinosis is the feature that fits least well. It is common in preterm infants for
reasons unrelated to the primary diagnosis - in one series, 23 of 152 infants born before 33
weeks, with respiratory distress and parenteral nutrition as the strongest associations, and
resolving in 76.5% by three months (PMID:40660023). It is treated here as plausibly secondary
to prematurity rather than as evidence for any gene.

**The maternal recurrent pregnancy loss is treated as phenotypic input, and it is doing work
in this report** - it is the reason chromosomal instability was the starting hypothesis rather
than a conclusion. Recurrent aneuploid conception is reported with maternal *BUB1B* variants
(PMID:42434306, PMID:42124334), which is consistent with a carrier parent; it is also what a
**balanced parental rearrangement** looks like, a differential that single-sample WGS cannot
exclude, and it is weak support for the two alleles lying in *trans* (section 8).

What it is deliberately **not** is a term in the similarity scoring, and the distinction is
mechanical rather than a judgement about relevance. HPO disease matching treats every term in
the input as a feature of the *proband*, so HP:0200067 *Recurrent spontaneous abortion* entered
as a proband term makes the algorithm search for diseases in which the affected child
miscarries - antiphospholipid syndrome, thrombophilias, uterine anomalies - and those displace
real candidates from a ranking in which the answer sits at rank 2. The term is therefore held
out of the scoring and recorded in `config/family_history.txt`, which no script reads - it is
there so the omission is visible to a reader rather than looking like an oversight. No ACMG
criterion follows from it either way; PP1 needs genotyped
relatives.

---

## 3. Pipeline

Six steps, each a single-file script with its own provenance output, preceded by the data
acquisition itself. Full commands in [README.md](../../README.md).

**Input.** The proband's GRCh38 VCF and its tabix index were obtained from the gated Hugging
Face dataset after access approval, and placed in `PATIENT_DATA/raw/`. Only the VCF was used:
BAM/CRAM are available in the dataset but this pipeline never reads them, which is why no
step here can inspect read-level evidence - the reason phase could not be resolved. 

| step | script | what it establishes |
| --- | --- | --- |
| 0 | `download_refs.py` | reference FASTA, Ensembl GFF3, HPO, ClinVar, gnomAD extracts, MANE |
| 1 | `audit_vcf.py` | build, contig naming, sample count, FILTER and FORMAT inventory |
| 2 | `annotate_vcf.py` | consequences via `bcftools csq`, ClinVar, gnomAD frequency |
| 3 | `build_panel.py` | phenotype-driven gene ranking, stamped before any candidate exists |
| 4 | `cascade.py` | genome-wide inheritance-model filtering |
| 5 | `rank_candidates.py` | claims, evidence, ordering, submission layout |

### How the ordering is assigned

The `epcr` column the challenge asks for is produced by one function, `tier_of()`, which reads
four things and nothing else: whether the gene is on the MVA checklist or in the phenotype
panel's top ten, whether at least one allele carries a ClinVar P/LP assertion *with review
criteria*, the gene's panel rank, and which inheritance arm the claim came from. It is a pure
function of recorded evidence - the same claim always produces the same tier - and the tiers
are fixed in the source, not fitted:

| tier | `epcr` | recessive arm: two alleles in a gene that is... |
| --- | --- | --- |
| A | 0.920 | on the checklist or top-10 in the panel, **and** has a ClinVar P/LP assertion |
| B | 0.500 | on the checklist or top-10 in the panel, with no such assertion |
| C | 0.150 | ranked in the panel's top 100 |
| D | 0.040 | present in the panel at any rank |
| E | 0.020 | outside the panel entirely |

The dominant arm is deliberately harder: a single heterozygous allele reaches tier A only if
the gene both fits the phenotype and carries the assertion, because for one allele in a
dominant gene the assertion *is* the evidence. Otherwise it lands at C if the gene ranks in
the panel's top 100, and E below that. Secondary findings bypass the tiers altogether and are
pinned at `epcr = 0.010`, beneath every primary candidate, because `epcr` asks for the
probability that the variant explains *this* phenotype and a secondary finding is by
definition not that.

Two properties of this scheme matter for reading the rest of the report. It never looks at the
NMD verdict, the gnomAD counts or the ClinVar star count as quantities - the assertion is read
as a yes/no - and it is **not** derived from the ACMG framework. What that costs is set out in
section 11.

### The annotation, checked

Splitting multiallelics took 4,895,719 input records to **4,962,060**, and **zero** were
dropped for a reference-allele mismatch against the FASTA - the check that catches a build
mix-up before it becomes a wrong answer. Field fill after annotation: `BCSQ` 68.2%, `gnomAD_AF`
92.8%, `CLNSIG` 0.89%. The low ClinVar figure is expected rather than a failure; most positions
in a genome have no clinical assertion. These figures are recorded in
`results/annotation_summary.json`, which is committed.

### The cascade, measured

Genome-wide, no region restriction:

| stage | records |
| --- | --- |
| all records | 4,962,060 |
| `FILTER = PASS` | 4,717,955 |
| quality gate (GQ ≥ 20, DP ≥ 10, allele balance) | 4,445,604 |
| population frequency | 170,569 |
| protein-altering consequence on a protein-coding transcript | **590** |

Thresholds: compound-heterozygous candidates at AF ≤ 0.001, homozygous at ≤ 0.005, and a
ClinVar P/LP assertion permitted up to 0.02. Heterozygous calls require allele balance between
0.25 and 0.75; homozygous calls require ≥ 0.85. chrY and chrM are excluded; sex is not
inferred, and chrX rows are flagged rather than filtered.

That left **118 genes with recessive candidates** and 5 with a dominant single-allele claim.
Of the eleven MVA genes, exactly one carried candidates: ***BUB1B***. A further 27 variants in
checklist genes that failed only on frequency or quality were written to
`set_aside_checklist.tsv`, so that a near-miss in a known gene is visible rather than lost.

### The phenotype panel

The seven features were scored against HPO disease annotations (release 2026-09-01) by Resnik
similarity over the most informative common ancestor, best-match average, patient-to-disease
direction. This produced a ranking of **1,449 genes**, in which ***BUB1B* ranked 2** - from
phenotype alone, before any variant was examined.

The genes immediately behind it are informative about how the ranking behaves. *BUB1*, *BUB3*,
*CEP57* and *TRIP13* follow at ranks **49-52**, and the route differs: *BUB1B* enters through
`OMIM:257300`, its own gene-specific MVA type 1 entry, while all four of the others enter only
through `ORPHA:1052`, the grouped MVA entry that carries no gene-level specificity. The panel
therefore recovers the MVA gene family as a block, and *BUB1B*'s lead over the rest of that
block comes from having a disease entry of its own rather than from anything about the child's
variants.

The obvious objection to any phenotype panel is that one built after seeing the candidate gene
can be made to rank it highly. An earlier version of this pipeline answered that by committing
a SHA-256 of the sorted term list before the cascade ran, as a pre-registration record. That
was dropped: the challenge does not require it, and hashing a set of seven terms drawn from a
small range of plausible paediatric features offers thin protection, since such a digest can be
enumerated. The objection is instead answered empirically in section 6 - discard the phenotype
altogether and the candidate set is unchanged and *BUB1B* still ranks first - which tests the
claim rather than attesting to it.

---

## 4. Evidence for the answer

### Allele 1, p.Leu737Ter - the null

A nonsense variant in exon 17 of 23, **8,861 bp upstream of the last exon junction**, so
nonsense-mediated decay applies and PVS1 holds at full strength. The NMD check reads exon
structure from the Ensembl GFF3 for the MANE Select transcript rather than assuming it;
an earlier version of that check sorted exons by coordinate, which inverts exon numbering on
the minus strand, and was corrected.

ClinVar records it as Pathogenic/Likely pathogenic with two review stars, for mosaic
variegated aneuploidy syndrome 1. Open Targets reports LOFTEE **high-confidence**. It is
present in gnomAD v4.1 at 115/1,461,846 with **no homozygotes**, consistent with a recessive
null allele carried in the population.

The absent homozygotes are not themselves evidence, and the table above should be read with
that in mind. At this frequency Hardy-Weinberg predicts **0.0045** homozygotes across the
730,923 individuals in the exome call set, so gnomAD has no power to observe one: `nhomalt = 0`
removes a potential refutation rather than supplying support. What the frequency data does
settle is the benign side. `AF_grpmax` stays at or below 1e-04 in **every** population group,
which closes BA1 and BS1 - something a global allele frequency on its own could not do, since a
variant can be globally rare and common in one group.

### Allele 2, p.Asn1002Lys - the uncertain one

A singleton: one allele in 1,461,878, absent from genomes. It is **not in ClinVar** and is not
indexed in Open Targets.

Codon 1002 is `AAT`, and `c.3006` is its third base - which explains the three ClinVar entries
at that position: T>C gives `AAC`, a synonymous change classified likely benign; T>A gives
`AAA`, p.Asn1002Lys, classified uncertain; and T>G, this child's allele, gives `AAG` - the
same protein change as the T>A allele, but a different nucleotide change and therefore not in
ClinVar. **PS1 and PM5 are both unmet**: PS1 needs an established pathogenic variant with the
same protein change, and the only same-change allele is itself uncertain.

**The call is therefore uncertain on two grounds**: the ACMG criteria reach 2 points, below the
6 needed for likely pathogenic, and the only other missense ClinVar holds at this codon agrees
- `c.3006T>A`, the same p.Asn1002Lys reached from a different nucleotide, is itself classified
uncertain. That agreement supplies no criterion, since PS1 would need it to be pathogenic, but
it does mean the uncertainty here matches a diagnostic laboratory's rather than following from
a thinner analysis.

It sits in the **pseudokinase domain (UniProt O60566, residues 766-1050)**, six residues from
the PLK1 phosphosite T1008. The domain is catalytically dead and functions as a structural
and protein-protein interaction module, so the available mechanistic question is stability,
not activity.

**Direction of effect.** Three independent lines say loss of function rather than gain:
recessive inheritance with a high-confidence null in trans; position inside a non-catalytic
stability domain; and a FoldX ddG of **+2.40 kcal/mol**, destabilising. For comparison, the
characterised MVA missense p.Leu1012Pro scores **+8.01** - 3.3-fold more destabilising - which
matches the mouse survival data independently: mice modelling the severe allele with a null in
trans died prematurely, while the equivalent of a milder allele was viable with MVA features
(PMID:31738183). Two calibrations belong beside that comparison. p.Leu1012Pro is the only
functionally characterised MVA missense in this region - ten residues away, same last exon,
same domain - it has a dedicated mouse model, and ClinVar still classifies it **uncertain**,
with criteria provided by multiple submitters and no conflicts. That is how far a missense in
this gene gets on present evidence. And PubMed returns no paper mentioning residue 1002 of
*BUB1B* at all, so the child's residue is uncharacterised rather than contested. In this
genotype a hypomorph and an antimorph are indistinguishable anyway,
because there is no wild-type protein left to poison.

**None of that changes the classification.** FoldX is not a ClinGen-calibrated pathogenicity
predictor, so it supports no ACMG criterion. Allele 2 remains a **variant of uncertain
significance at 2 points**, and the report says so.

### Why the pair is still the answer

Matsuura and colleagues established that MVA requires BUBR1 to fall below roughly half of
normal (PMID:16411201), and Hanks and colleagues described truncating *and* missense *BUB1B*
variants across five MVA families, two of them with embryonal rhabdomyosarcoma
(PMID:15475955). A null plus a destabilising missense is the configuration those papers
describe. The genotype is also *not* homozygous null, and that matters: complete BUBR1 loss is
not seen in living patients, which is why nearly every reported case pairs a null with a
hypomorphic allele.

**The two sources do not agree about this genotype, and saying so is worth more than picking
the convenient one.** Hanks shows truncating-plus-missense recurring in living human MVA
patients. Sieben's mouse allelic series shows the analogous combination - `X753/L1002P`, a null
in trans with this missense, which is precisely the child's architecture - dying prematurely,
while the viable combination was hypomorph-plus-missense. The reconciliation that fits both is
a hypothesis and not a result: that the particular missense decides whether the pairing is
survivable, and that p.Asn1002Lys retains more residual function than p.Leu1012Pro does. On
that reading the child being alive is evidence about the allele rather than against the
architecture. The mouse abstract does not give the timing of death, so this report does not
call it embryonic or perinatal.

---

## 5. Secondary findings

A single-allele finding is reported only if it is truncating **and** NMD-competent **and**
carries a ClinVar P/LP assertion with at least one review star. Four candidates were tested
and one survived: *LZTR1*. *RBM8A* was rejected (5' UTR), *GNRHR* (missense) and *FLG*
(escapes NMD - last exon 3 of 3). Each rejection is logged with its reason.

**LZTR1 is not on the ACMG SF v3.3 list** (Genet Med 2025;27(8):101454), which was read
rather than assumed - none of the four candidates is, and v3.3 added exactly three genes while
declining three others. It is reported anyway, and the reason is worth stating because the two
facts together look contradictory until the rule behind them is named.

The filter above is mechanical: truncating, NMD-competent, P/LP with review stars. What
actually separates *LZTR1* from the other three is that **its zygosity matches its inheritance
mechanism.** *LZTR1*-related schwannomatosis is autosomal dominant, so one truncating allele
*is* the disease configuration. The others fail that test: *GNRHR* hypogonadotropic
hypogonadism is recessive, making a single heterozygote a carrier - a category ACMG SF
deliberately excludes; *RBM8A* TAR syndrome needs the hypomorph in trans with a 1q21.1
deletion; *FLG* ichthyosis vulgaris is semi-dominant with a minimal heterozygous phenotype.

That the gene is absent from ACMG SF is a statement about a *reporting list*, not about the
variant. The list carries *NF2* as NF2-related schwannomatosis; *LZTR1* is the other
schwannomatosis gene and was simply not added. A truncating P/LP allele in a
tumour-predisposition gene, in a child who already has a rhabdomyosarcoma, is worth surfacing
on those grounds alone - and the Hackathon rules invite incidental findings and score them as
harmless.

So the `finding_type` label is a pipeline default rather than a clinical recommendation, and
`ranking.json` says so explicitly. In a clinical report under the current professional
standard none of these four would be returned as a secondary finding; this is a research
submission, and the row is flagged `secondary` at `epcr = 0.010` precisely so it cannot
compete with the primary claim.

---

## 6. How much does the answer depend on the phenotype?

This is the question a phenotype-driven pipeline has to answer about itself, so it is built
into the tool as a switch rather than argued in prose. `rank_candidates.py --panel-mode`
re-runs the ordering with progressively less prior knowledge. **The candidate set is identical
in all three runs - 123 claims - because the cascade has already finished before any of this
is read.** Only the order changes.

| mode | knowledge used | *BUB1B* tier | *BUB1B* epcr | rank |
| --- | --- | --- | --- | --- |
| `full` | phenotype panel + MVA gene list | A | 0.920 | **1** |
| `no-panel` | MVA gene list only, phenotype discarded | A | 0.920 | **1** |
| `agnostic` | neither | E | 0.020 | **1** |

***BUB1B* comes first in all three.** Stripped of the phenotype *and* the disease gene list,
it still wins, on variant evidence alone: a ClinVar-asserted null paired with a rare missense
in the same gene.

The honest reading of the third row is that the phenotype is **confirmatory, not selective** -
but it is what makes the claim *confident*. In agnostic mode the separation collapses: 0.020
against 0.019 for the next candidate. The panel does not choose the answer; it is what turns a
narrow lead into a defensible one.

---

## 7. What the pipeline got wrong

These were found by running on real data, not by reading the code, and they are listed because
a methods section that reports no errors has usually not looked.

1. **Transcript biotype was not required.** `bcftools csq` reports per transcript, and
   consequences were being accepted on NMD, retained-intron and non-stop-decay transcripts.
   This produced a spurious *RPS17* "compound heterozygous pair of missense variants" in which
   one allele was actually **synonymous** and neither touched MANE Select. Requiring a
   protein-coding transcript removed *RPS17* and *KMT2C*, cut 749 rows to 590 and dropped 164
   gene-level claims. ClinVar P/LP assertions are the single declared exception and are
   counted separately.
2. **The NMD check inverted exon numbering on the minus strand**, because it sorted exons by
   coordinate rather than by transcript order. Caught on *FLG*.
3. ***BUB1B* appeared twice** - once as a dominant single-allele claim, once as the pair - and
   the single allele sorted above the pair, putting half the answer at the top of the
   submission. Claims are now collapsed per gene, keeping the pair.
4. **Absence of evidence was reported as evidence of absence.** The gnomAD v4.1 extracts cover
   the checklist genes only; for any other gene the notes now say the frequency was not looked
   up, rather than reporting it as absent. `v41_region_covered` records which.

---

## 8. Limitations

- **Phase is unproven, and *trans* is an assumption.** The two alleles are 10,911 bp apart -
  too far for short-read phasing from a single library, well within range for a long read or
  for parental testing. Reporting the pair as compound heterozygous asserts that the alleles
  sit on opposite copies; in *cis* the second copy is intact and the genotype explains nothing.
  The maternal pregnancy-loss history is weakly consistent with *trans* (section 2), but that
  is an argument about which parent carries what rather than phase evidence. Parental
  genotyping settles it, and is next step 3.
- **A balanced parental rearrangement is not excluded** by single-sample WGS, and the maternal
  pregnancy-loss history is compatible with one.
- **No structural or copy-number analysis.** Without a BAM or CRAM there is no coverage track,
  so a deletion in trans with a point variant would be invisible.
- **PP3 is not assessed** - by design, see constraint 3.
- **PP4 is not claimed** for either allele. The phenotype is consistent with MVA, but the panel
  already uses the phenotype to order candidates, so counting it again as independent evidence
  for the same gene would be circular.

---

## 9. What should be done next

In order of value, and the first two are measurements rather than analyses:

1. **BUBR1 western blot** on patient cells against a control. This tests the premise the whole
   mechanistic argument rests on - that the missense allele reduces protein - and is the only
   realistic route to PS3, which is what would move allele 2 out of uncertain significance.
2. **Premature chromatid separation rate** by G-banding with centromere FISH. This is the
   cytogenetic hallmark of MVA, it is the cheapest confirmatory test available, and it is
   absent from the work-up.
3. **Parental genotyping** - resolves phase, converts PM3 from supporting to a stronger tier,
   and addresses the balanced-rearrangement differential.
4. **Tumour surveillance** per a cancer-predisposition protocol, given the rhabdomyosarcoma
   and the Wilms tumour risk reported in MVA.

---

## 10. Reproducing this

```bash
uv sync
./src/download_refs.py                          # ~7 GB into refs/

export VCF=PATIENT_DATA/raw/<your-file>.vcf.gz
./src/audit_vcf.py     --vcf "$VCF"
./src/annotate_vcf.py  --vcf "$VCF" \
    --fasta refs/GRCh38.fa \
    --gff   refs/Homo_sapiens.GRCh38.116.chr.gff3.gz \
    --clinvar refs/clinvar.vcf.gz \
    --gnomad  refs/af-only-gnomad.hg38.ensembl.vcf.gz \
    --cache refs/

$EDITOR config/patient_hpo.txt                  # the proband's features
./src/build_panel.py --leave-one-out

./src/cascade.py
./src/rank_candidates.py --top 5 --proband-id <id>
./src/rank_candidates.py --panel-mode agnostic  # the robustness run in section 6
```

The panel step comes before the cascade because the two are independent, not because the order
is enforced: `--panel-mode agnostic` re-runs the ordering with the phenotype discarded and is
what makes "the ranking was not tuned against the answer" a checkable claim.

Every step writes a provenance file beside its output: `vcf_audit.json`,
`annotation_summary.json`, `panel.meta.json`, `cascade_counts.json`, `ranking.json`. Every
number in this report is read from one of them. Only `annotation_summary.json` is committed,
being a non-identifying aggregate; the other four are derived from the proband's genotype or
phenotype and stay out of version control, which `build_panel.py` enforces for
`panel.meta.json` on every run. A reviewer therefore reproduces the figures by re-running the
pipeline on the dataset rather than by reading them out of this repository.

---

## 11. Declarations

**ACMG classification is a manual annotation.** The criteria, their strengths and the resulting
class in section 4 were assigned outside the code; nothing in the pipeline produces or scores a
classification. Both verdicts are checkable against ClinVar and both check out: allele 1 agrees
with a two-star entry carrying criteria from multiple submitters and no conflicts, and the
uncertain call on allele 2 agrees with ClinVar's assessment of the same protein change reached
from a different nucleotide. The classification deliberately does not use PP5, the
reputable-source criterion ClinGen recommends against on the grounds that interpretation should
rest on primary data rather than on assertions decoupled from the evidence behind them
(PMID:29543229). The Pathogenic call is carried by PVS1 - the NMD verdict computed here from
the Ensembl GFF3 - together with population frequency, and the ClinVar agreement is reported as
corroboration rather than scored.

**The `epcr` column is a prioritisation weight, not a classification.** `tier_of()` reads four
recorded facts - checklist membership or a panel rank in the top ten, a ClinVar P/LP assertion
carrying review criteria, the gene's panel rank, and the inheritance arm - and is
deterministic, so the same claim always yields the same tier. These are the ordinary grounds
for deciding which candidate to examine first. The two numbers in this report answer different
questions and neither feeds the other: ACMG points measure how strong the evidence for an
allele is, `epcr` which claim to look at first. Because the tier reads the ClinVar assertion as
a yes/no rather than scoring the NMD and frequency evidence quantitatively, `epcr` is an
ordering weight and not a calibrated probability. Making the tier read the NMD verdict directly
would align the two, and would change every `epcr` value in the submission, so it is left as
future work.

**AI use.** Anthropic Claude, Team plan subscription - commercial terms, customer content is not
used to train models - model `claude-opus-5`, used through Claude Science for pipeline code,
literature retrieval and drafting.
