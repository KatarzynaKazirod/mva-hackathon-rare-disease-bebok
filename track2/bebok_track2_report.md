# Track 2 - drug repurposing in BUB1B-related mosaic variegated aneuploidy

**Team Bebok** - Rare Disease, Real Kid: The MVA Hackathon 2026, Track 2.

Hypotheses for follow-up, not evidence that any medicine works.

Every drug named below carries a ChEMBL identifier retrieved from Open Targets, and every
mechanistic claim carries a PMID verified through NCBI E-utilities. No candidate was
admitted from recollection. Supporting tables are in `data/` and `results/`.

**The search was re-examined after a first pass and the corrections are folded into the tables
below**, so this document is the current state rather than a first draft. What the re-examination
changed, and why, is section 10. The methods description required by the organisers is
Appendix A.

**Supporting material in the repository**

| file | what it adds |
| --- | --- |
| [results/ranked_candidates.tsv](results/ranked_candidates.tsv) | the shortlist, scoring components beside the total |
| [results/validation_plan.md](results/validation_plan.md) | the experiment that would refute each hypothesis |
| [results/expression_evidence.md](results/expression_evidence.md) | public cell and animal data, including a result retracted on inspection |
| [data/axes_considered.tsv](data/axes_considered.tsv) | every therapeutic axis examined, with its verdict and reason |

**Every count in this report is falsifiable in one command.**

```bash
python3 src/run_pipeline.py --mode check      # add --offline to skip the network stages
```

That rebuilds the tables from the live Open Targets API into a temporary directory, diffs them
against what is committed without touching the repository, and re-resolves every cited PMID
against PubMed. It exits non-zero on any difference. This matters more here than it would for a
positive result: the central finding is a **negative** one - that no approved medicine acts on
any node of the broken machinery - and a negative claim rots quietly when the database beneath
it moves to a new release. The check turns that drift into a failure instead of a silence.

## 1. What is broken

The proband carries a compound-heterozygous BUB1B genotype on MANE Select `NM_001211.6`:
`c.2210T>G` **p.Leu737Ter** (pathogenic, 10 ACMG points) in presumed *trans* with
`c.3006T>G` **p.Asn1002Lys** (uncertain significance, 2 points). The first is NMD-competent -
exon 17 of 23, 8,861 bp upstream of the last junction - so it contributes no protein. The
net state is a **reduced dose of partially functional BUBR1**, not its absence; complete
loss is embryonically lethal and the patient is alive.

UniProt O60566 places the protein kinase (pseudokinase) domain at **766-1050** and the
PLK1-phosphorylated KARD serine S676 at 676. Two consequences follow, both checkable:

- p.Leu737Ter stops at 737, **before** the domain begins at 766, so the truncated product
  would lose the entire pseudokinase module while retaining KARD.
- p.Asn1002Lys lies **inside** that domain, six residues from the PLK1 phosphosite T1008.

Both alleles therefore converge on the C-terminal stability function rather than on
checkpoint signalling.

### Is the second allele loss- or gain-of-function?

Missense pathogenicity predictors cannot answer this: REVEL, SIFT and PolyPhen are trained on
pathogenic-versus-benign labels and a gain-of-function allele scores as readily as a
loss-of-function one. The question needs a **mechanism-specific** method. Three lines point to
**hypomorph - partial loss of function through reduced protein stability** - but they are not
three independent lines, and the first pass of this report claimed they were. Line 1 is an
expectation that presupposes the diagnosis in order to establish the allele's mechanism; lines
2 and 3 are the same observation at two levels of description, since the folding-energy
estimate is computed on the structure whose domain assignment is line 2. Read as **two partly
independent lines plus an inheritance-based expectation**:

1. **Inheritance.** MVA1 is autosomal recessive and allele 1 is an amorph (NMD-competent
   nonsense; Open Targets reports **LOFTEE high-confidence LoF**). For a recessive
   loss-of-function phenotype the second allele must also reduce function. A hypermorph or
   neomorph in *trans* with a null does not produce recessive disease.
2. **Location.** Residue 1002 lies inside the pseudokinase domain (766-1050, UniProt
   O60566), which is catalytically dead in humans and contributes protein stability. Damage
   there predicts *less protein*, not altered activity.
3. **Biophysics.** Open Targets carries a **FoldX** stability prediction of
   **ddG = +2.40** for p.Asn1002Lys - destabilising. The proband's `c.3006T>G` is itself a
   singleton and not indexed, but the catalogued `c.3006T>A` produces the identical amino-acid
   substitution, and predictors act on the protein change rather than the nucleotide; the
   record confirms `aa=N1002K` on canonical ENST00000287598, a fourth independent check on the
   MANE numbering. A synonymous control at the same base (`c.3006T>C`) returns no FoldX score,
   as it should.

**The direction calibrates against the one characterised neighbour; the magnitude does not.**
Human **p.Leu1012Pro** - ten residues away, and the allele modelled in mouse, where it reduces
total BUBR1 - scores **ddG = +8.01** against +2.40 here. The defensible reading is **"both
destabilising, the proband's allele less so"**. The first pass quoted the ratio as roughly
3.3x, which over-interprets the number: the scale is non-linear, the method's error is on the
order of +/-0.5 kcal/mol, and the comparator is a substitution to proline, which the method
penalises predictably on backbone grounds. A ratio computed from two such values is not a
3.3-fold difference in anything measurable.

> **Numbering warning, because two residue numbers for one allele appear here.** The Sieben
> abstract states verbatim: *"we mimicked the BUBR1L1012P mutation in mice (BubR1L1002P)"*.
> **Human L1012P and mouse L1002P are the same allele** under a 10-residue species offset,
> and ClinVar confirms the human form as `NM_001211.6:c.3035T>C` p.Leu1012Pro. Mouse genotypes
> are quoted in mouse numbering (`X753/L1002P`, `H/L1002P`) and must not be read as touching
> the proband's residue: **p.Asn1002Lys is residue 1002 in HUMAN numbering**, a different
> residue from mouse L1002P. The two numbers merely collide.

The FoldX comparison independently supports the reconciliation proposed from the mouse
survival data: `X753/L1002P` animals (null plus missense, the proband's genotype architecture)
died prematurely while this child is alive, which requires N1002K to retain more residual
function than L1012P. Mouse phenotype and biophysical prediction agree without being told to.

Two limits on this. FoldX is **not** a ClinGen-calibrated PP3 tool, so the ACMG arithmetic is
unchanged at 2 points, uncertain significance - what improves is the *mechanistic assumption*,
not the classification. And in this genotype **hypomorph and antimorph are experimentally
indistinguishable**: allele 1 produces no protein, so every BUBR1 molecule in the cell is
already N1002K and there is no wild-type partner for a dominant-negative to poison. That
distinction could only be resolved in a heterozygous carrier - which is where the recurrent
pregnancy loss phenotype of Wei 2026 (PMID:42434306) sits. See
`data/variant_effect_predictions.tsv`. Matsuura 2006 (PMID:16411201) placed the
pathogenic threshold at a **>50% reduction** in BUBR1 - a number that matters below, because
it means the therapeutic goal is to cross a threshold, not to restore wild-type protein.

BUBR1 does two jobs and both are impaired: it contributes the inhibitory MCC signal that
restrains APC/C-CDC20, and through KARD it recruits PP2A-B56 to counterbalance Aurora B at
kinetochores (PMID:36441015). A checkpoint defect alone explains at most half the phenotype.

## 2. The therapeutic goal has to be reframed

A germline hypomorph cannot be pharmacologically corrected in every cell of a child. The
question is not "which drug fixes BUB1B" but "which layer of the consequence is tractable".
We treated that as an empirical question rather than a rhetorical one, and the answer is
unusually clear-cut.

### Two alleles, two different interventions

A consequence of the compound-heterozygous state that a single-drug answer cannot accommodate:
the two alleles are broken in different ways and would need different corrections. Allele 1 is a
premature stop, so the only allele-directed route is translational readthrough - which requires
NMD inhibition to produce anything at all, and whose approved representatives are excluded here.
Allele 2 is a destabilised missense protein, so its route is proteostatic: raise the fraction
that survives quality control. Nothing acts on both. A proposal naming one medicine for this
child is therefore either addressing one allele and silent about the other, or acting downstream
of both - and the report says which, for every candidate.

## 3. Direction of action, and why it is a gate

Network proximity alone recommends harmful drugs here. MPS1/TTK **activates** the spindle
assembly checkpoint, so an MPS1 inhibitor - well-evidenced, tractable, abundantly
characterised - would rank highly on any unsigned method and would further weaken a
checkpoint that is already hypomorphic. The same applies to PLK1 inhibitors: PLK1
phosphorylates the BubR1 KARD (PMID:39516273, PMID:33248027), so inhibiting it removes what
remains of PP2A-B56 recruitment.

We therefore curated a signed table (`data/signed_network.tsv`, one PMID per row) giving the
required **direction** of perturbation per node, and used it as a **hard filter** rather than
a scoring dimension - because a drug with strong evidence and the wrong sign would otherwise
float to the top.

**The finding that direction is tissue-dependent.** Somatic tissue needs segregation fidelity
restored. The rhabdomyosarcoma is already chromosomally unstable and the therapeutic aim
there is the opposite: push it past the tolerable threshold. The same table carries two
direction columns, and drugs sort into two non-overlapping groups accordingly. MPS1
inhibition is systemically contraindicated and tumour-directed rational - simultaneously.

STRING supplied topology; it is unsigned and undirected and we do not pretend otherwise. The
sign is curated. Reactome's REST API returned Cloudflare 1010 to our client, so pathway
membership came from the project's own five-layer catalogue validated gene-by-gene against
Open Targets - 65 of 78 candidate symbols confirmed, four mis-resolved tokens rejected.

## 4. The target-based route is empirically closed

Of **65 spindle-checkpoint and kinetochore nodes queried, 60 have no drug of any kind** -
including every core component: BUB1B itself, BUB1, BUB3, CDC20, MAD1L1, MAD2L1, MAD2L1BP,
TRIP13, the entire APC/C, the entire CENP platform, the entire KMN network.

The five exceptions carry 49 drug records between them: AURKA 19, AURKB 18, PLK1 8, BIRC5 2,
TTK 2. Their maximum clinical stage distribution is PHASE_1 29, PHASE_2 14, PHASE_2_3 1,
PHASE_3 5. **Approved: zero.**

Applying the direction gate leaves **19 of the 49 compounds acting correctively in somatic
tissue, across two nodes only**: 18 AURKB inhibitors and one BIRC5 antisense inhibitor
(gataparsen, CHEMBL2219774). Both work through the same argument - partial reduction restores
the kinase-phosphatase balance that reduced KARD-mediated PP2A-B56 recruitment has tipped,
whether the reduction is applied to Aurora B itself or to the survivin that localises the CPC.
Not one of the 19 is approved. The most advanced reach phase 3 - ibcasertib (CHEMBL5095033) and AT-9283 (CHEMBL495727); barasertib, often named first in this class, is recorded at PHASE_2_3, which sits below phase 3 in trial-phase order.

The remaining 30 compounds fail the somatic gate: 19 AURKA inhibitors (that node is signed
intact, so the required direction is neutral), 8 PLK1 and 2 TTK inhibitors (both nodes need
their function preserved or raised, and an inhibitor lowers it), and one survivin peptide
vaccine, which fails in both contexts because its action type raises an immune response rather
than lowering target function.

So: **no approved medicine acts on any node of the machinery this disease breaks.** That is a
result, not a gap in the search, and it forces every viable hypothesis downstream.

An important corollary for the therapeutic window: the aim at AURKB would be *partial*
rebalancing, not maximal inhibition, because over-correcting the checkpoint arrests cells in
mitosis. A hypothesis that defines its own ceiling is more credible than one that does not.

### The cohesion layer, checked separately

The node set above covers the checkpoint, its silencing, the kinetochore platform and error
correction. Sister-chromatid cohesion was deliberately left out as upstream machinery - but
**premature chromatid separation is the cytogenetic hallmark of BUB1B-MVA**, so the exclusion
had to be tested rather than assumed. It was, and it costs nothing:

**All 17 cohesion-layer targets return zero approved drugs, and 16 of the 17 return no
compound of any kind** (`data/cohesion_layer_census.tsv`) - cohesin itself (RAD21, SMC1A,
SMC3, STAG1, STAG2), its loader and regulators (NIPBL, PDS5A, PDS5B, CDCA5, WAPL), the
protection machinery (SGO1, SGO2, PPP2R5A, PPP2R5B), and **ESPL1 (separase)**.

Separase deserves naming. It is the protease that cleaves cohesin, so its premature action
*is* the phenotype, and proteases are ordinarily a tractable class - yet Open Targets records
no compound against it at all. The single most specific conceivable intervention against
premature chromatid separation has no chemical matter behind it.

A precision point that must not be fudged: **BUB1**, not BUBR1, is the kinase that
phosphorylates H2A-T120 to recruit shugoshin and protect centromeric cohesion. BUBR1 is a
pseudokinase and does not perform that step, so BUB1's cohesion role cannot be transferred to
it. How reduced BUBR1 produces premature chromatid separation is therefore **not settled** in
what this project verified; plausible routes are competition for shared PP2A-B56 pools, the
SAC-independent attachment-stabilisation defect, or cohesion fatigue secondary to abnormal
mitotic timing. Those are hypotheses and are labelled as such.

### Asking the disease rather than the gene

The census above interrogates the gene's machinery. The same question put to the disease is an
independent route to the same place, and it agrees.

| MVA entries in MONDO (the grouped entry and every subtype) | with any medicine recorded | approved |
| --- | --- | --- |
| **9** | **0** | **0** |

Walking from disease to genetically associated target to drug returns **11 genes** above the
association-score floor, **9** of them on the curated eleven-gene MVA checklist, and exactly
**one** carrying any compound at all (`data/mva_disease_drugs.tsv`, `data/mva_group_genes.tsv`).

That walk also contains a trap, which the pipeline reports per row. **PAK6** associates with MVA
at 0.52 and is not an MVA gene: it overlaps *BUB1B* at chr15:40,217,428-40,277,489, with a
readthrough locus annotated across both, so the association is positional rather than biological.
It carries no drugs and nothing propagated from it - but a pipeline that walked target to drug
without checking coordinates would have proposed a PAK inhibitor on the strength of a transcript
annotation.

## 5. The consequence arm

Because these drugs do not target network nodes, they cannot be generated from the network.
They are literature-derived hypotheses **verified against databases**. Of eight candidates
carried forward, seven propose a mechanism that is not the drug's registered mechanism of action.
**That is not in itself a criticism, and framing it as one was an error of ours.** Using a medicine
outside its registered indication, on the strength of an activity documented elsewhere, is precisely
what drug repurposing means - a urea-cycle drug proposed here for chaperone activity is the intended
logic, not a deviation from it.

What the check is actually good for is separating three different evidential situations that the
word "mismatch" collapses into one: the registered mechanism IS the lever invoked; the activity is
**established pharmacology but not the registered mechanism** (phenylbutyrate as a chaperone,
hydroxychloroquine blocking autophagic flux - both well documented, neither registered); or the
activity is **asserted for that drug** without being part of its established pharmacology. Only the
last is a defect, and it is what this check removed: fingolimod and perphenazine were both nominated
as PP2A activators on review-literature grounds that neither drug's primary pharmacology supports -
and in fingolimod's case the dominant registered action, S1P-receptor agonism, is immunosuppression
in a child with an active sarcoma.

| candidate | ChEMBL | registered mechanism | proposed mechanism | arm |
|---|---|---|---|---|
| sodium phenylbutyrate | CHEMBL1746 | **none recorded** | chemical chaperone for the destabilised missense | somatic |
| hydroxychloroquine | CHEMBL1535 | TLR7/TLR9 antagonist | autophagy inhibition in aneuploid cells | **tumour only** |
| ataluren | CHEMBL256997 | 80S ribosome modulator | readthrough of the TGA stop | somatic, blocked |

> **Status update after product-level verification.** openFDA returns one record for ataluren and it is `BULK INGREDIENT`: no NDA, no ANDA, so there is no approved product in the United States, and the EU conditional authorisation was reportedly not renewed in March 2025. This candidate may fail the rule that only medicines already in human use qualify, and is flagged rather than deleted pending an EMA check.

| acetylcysteine | CHEMBL600 | none recorded | antioxidant buffering | somatic |
| metformin | CHEMBL1431 | complex I inhibitor | metabolic stress in aneuploid cells | tumour only |
| dasatinib | CHEMBL5416410 | BCR-ABL inhibitor | senolytic clearance | **contraindicated** |
| fingolimod | CHEMBL314854 | S1P receptor agonist | PP2A activation | **rejected** |
| perphenazine | CHEMBL567 | dopamine D2 antagonist | PP2A activation | **rejected** |

> **Moved to contraindicated on re-examination.** Dasatinib was first recorded as
> tumour-directed. Delayed epiphyseal closure, osteopenia and growth suppression are recorded
> paediatric harms, which puts the drug's principal harm and this child's principal complaint on
> the same axis. Rapalogs were considered in the same arm and are contraindicated for a
> mechanistic reason: TORC1 inactivation promotes **APC/C-dependent mitotic slippage**
> (PMID:35141499), and the lesion here is already a weakened MCC inhibitory signal toward
> APC/C-CDC20, so the drug pushes the lever the disease is already pushing. Both are listed in
> section 6 and recorded in `data/contraindicated.tsv`.

Fingolimod and perphenazine were the mechanistically most attractive idea in this arm -
boosting PP2A to compensate for lost KARD recruitment. Both are rejected: the PP2A activity
is literature-only and absent from the registered mechanism, fingolimod is immunosuppressive
in a child with an active sarcoma, and perphenazine carries a black box warning.

### Shortlist

**Arimoclomol** is the highest-ranked candidate and the best mechanistic fit on this arm. It is a
heat-shock response co-inducer: rather than chaperoning one protein, it amplifies the cell's own
chaperone response under stress, which is the behaviour a destabilised missense allele needs if
the goal is to raise the steady-state level of a partly functional protein past a dosage
threshold. It is an approved medicine - FDA NDA, marketed as a capsule since September 2024,
with paediatric experience from 2 years of age - so it clears this project's rule that only
medicines already in human use qualify.

Three objections are recorded against it. **First, its registration is for use in combination
with miglustat**, in a different disease; using it here would be off-label both in indication and
in regimen, which is a further step than off-label indication alone. Second, no mechanism of
action is recorded for it in the target database, so the proposed mechanism is literature-based,
not database-verified - and the database still records it at PHASE_3, which is a release-lag
artefact rather than a disagreement about its status. Third, and shared with every candidate on
this arm, it has never been tested in a chromosomal-instability system, so the chaperone
hypothesis rests on an inference about protein stability that the missing protein-level
measurement would settle.

**Sodium phenylbutyrate** is the second candidate aimed at the causal lesion. If p.Asn1002Lys
destabilises the protein, raising its steady-state level past Matsuura's 50% threshold is a
lower bar than restoring wild-type protein. Approved, oral, with chronic paediatric
experience.

Its weakness is stated precisely, because the loose version of it was wrong. The drug is not
without a mechanism - it has several, and that is the problem. ChEMBL records exactly one
*registered* mechanism, **glutamine sequestering**, the nitrogen-scavenging route behind its
urea-cycle indication; Open Targets surfaces no mechanism row at all, which is a coverage gap
rather than an absence. The two activities we care about are both well documented in the
literature and neither is registered: **chemical chaperone** activity, which is the lever we
propose, and **non-competitive HDAC inhibition**, which is not.

**Those two activities pull in opposite directions on this target, and that sharpens the
counter-evidence into an explanation.** The chaperone activity should raise the fraction of
destabilised p.Asn1002Lys surviving quality control. The HDAC activity alters transcription -
and an HDAC-inhibition signature is precisely what we observed in GSE86928 (`HIST1H1C` +2.08
with `MKI67` -3.7, `CCNA2` -4.6, `FOXM1` -3.0) alongside *BUB1B* falling 2.4-fold, while our
signed network requires *BUB1B* to rise. So the expression result is not an unexplained
coincidence: the same molecule plausibly causes it through its second pharmacology, and a
proteostatic gain would have to overcome a transcriptional loss the drug itself produces.

The limits on that counter-evidence stand unchanged - n=1 per condition, 5 mM is
supra-therapeutic, A549 is a cancer line, and expression is the wrong assay for a
post-translational hypothesis. It does not refute the candidate; it relocates the burden of
proof onto the protein-level measurement. The drug has also never been tested in any
chromosomal-instability system.

**And public expression data argue against it at the exposure tested.** In GSE86928
(A549 cells, 5 mM sodium phenylbutyrate, PMID:28033672) **BUB1B transcript falls 2.4-fold**,
alongside 31 of 32 other spindle-checkpoint and kinetochore genes. The drop is secondary to
cell-cycle arrest - `MKI67` -3.7, `CCNA2` -4.6, `FOXM1` -3.0, with `HIST1H1C` up, the
signature of HDAC inhibition - so it is not a specific effect on the checkpoint. But our own
signed network requires BUB1B to go **up**, and a proteostatic gain would have to overcome a
transcriptional loss in any proliferating tissue. The study is n=1 per condition at a
supra-therapeutic concentration in a cancer line, so it cannot refute the hypothesis; it
does relocate the burden of proof onto the protein-level experiment. See
`results/expression_evidence.md`.

**Hydroxychloroquine**, tumour arm only. Our own direction logic disqualifies it
systemically - somatic tissue needs autophagy to tolerate aneuploidy.

**Ataluren** is the most genotype-specific hypothesis and currently unexecutable. p.Leu737Ter
is `TTA` -> `TGA`, and TGA is the most readthrough-permissive stop codon - but the transcript
is NMD-degraded, so readthrough needs an NMD inhibitor and none is approved. The obvious
approved readthrough route, aminoglycosides, is excluded here on nephrotoxicity grounds.

## 6. What must not be given

Derived from the same direction logic (`data/contraindicated.tsv`):

- **aminoglycosides** - nephrotoxic against existing nephrocalcinosis
- **growth hormone** - growth-promoting in a cancer-predisposition syndrome with an active
  sarcoma; the intuitive answer to short stature and failure to thrive, and the wrong one
- **MPS1/TTK and PLK1 inhibitors** - systemically contraindicated, tumour-directed only
- **vinca alkaloids and taxanes** - they prolong mitosis *by activating the SAC*, which is
  the broken mechanism, and are aneugenic in a constitutionally unstable genome
- **fingolimod** - immunosuppressive under cancer predisposition
- **rapalogs (sirolimus, everolimus)** - TORC1 inactivation promotes APC/C-dependent mitotic
  slippage (PMID:35141499), the same lever the weakened checkpoint already fails to restrain;
  any muscle or growth benefit would be bought by worsening the primary defect
- **dasatinib** - delayed epiphyseal closure, osteopenia and growth suppression are recorded
  paediatric harms, on the same axis as this child's growth failure

## 6b. Approaches outside the repurposing scope

The track asks for existing approved medicines, and sections 1-6 answer that. Three further
approaches fall outside that scope while bearing directly on this child, and the first two
are arguably more consequential than anything in the shortlist.

### Exon skipping - frame-compatible, and the arithmetic is favourable

p.Leu737Ter sits in **exon 17 of 23**. From the project's own GFF3 (Ensembl 116, MANE Select
ENST00000287598):

```
CDS 17: 40209635-40209775    length 141 nt    141 % 3 = 0
coding bases before it: 2143      total CDS: 3153 nt = 1050 residues
```

**141 is divisible by 3, so exon 17 is in-frame skippable.** Skipping it removes 47 codons and
yields an internally shortened but in-frame product rather than a premature stop - which
**abolishes NMD** and converts a null allele into one that makes protein. That is precisely
the direction this entire report searched for and did not find: *raising* BUBR1 dose rather
than inhibiting something.

The deleted segment is approximately **residues 715-761**. KARD ends at 682; the pseudokinase
domain begins at 766. Exon 17 therefore encodes the **linker between the two modules the
mechanism depends on**, and touches neither.

Limits, stated so this does not read better than it is: exon-skipping antisense
oligonucleotides are approved technology (Duchenne) but are always gene- and exon-specific,
and nothing exists for BUB1B. Skipping efficiency would have to be sufficient to matter. The
function of an internally deleted BUBR1 is **unknown** - the absence of an annotated module in
715-761 does not mean the linker is dispensable. And it addresses only allele 1.

### Tumour surveillance - the intervention with actual evidence of benefit

In cancer-predisposition syndromes the proven intervention is a **surveillance protocol, not a
medicine**, and this report said nothing about it until now.

The specific case for renal surveillance in this child: a girl with **MVA syndrome 3** had a
Wilms tumour resected at age 3 and presented with **orbital rhabdomyosarcoma nine years
later**, at 12 (PMID:42595739). Wilms tumour is documented bilaterally in premature-chromatid-
separation / MVA (PMID:31081598), and biallelic TRIP13 - another MVA gene - predisposes to
Wilms tumour directly (PMID:28553959). The tumour spectrum tracks the mechanism rather than
the individual gene.

For a child who already has a rhabdomyosarcoma, the implication is that **it may not be the
only primary**, and renal ultrasound is cheap and non-invasive. General framework for
paediatric predisposition surveillance: PMID:42784470. This is probably the most useful item
in Track 2, and it is not a drug.

### Other modalities, named with their blockers

- **Base or prime editing** - both alleles are single-nucleotide substitutions and are
  conceptually correctable; `c.2210T>G` is a transversion, so prime editing rather than a
  classical base editor. Not approved, realistically *ex vivo*, and the requirement here is to
  reach every dividing cell in the body.
- **Gene addition (AAV)** - the 3,153 nt CDS fits the capsid, but a mitotic gene needs
  cell-cycle-regulated expression and systemic delivery to all proliferating tissue is not
  realistic.
- **Allogeneic haematopoietic stem cell transplant** - replaces the marrow compartment with
  genetically normal cells and is considered in some chromosomal-instability syndromes. It
  would address the leukaemia risk in the MVA tumour triad but **does nothing for solid
  tumours**, and this child's malignancy is solid.
- **Family-level**: parental carrier testing (needed anyway for phase and PM3) and
  preimplantation genetic testing for future pregnancies, given a 25% recurrence risk and the
  mother's recurrent pregnancy loss.

### Considered and rejected: immune checkpoint inhibition

The tempting logic is that aneuploid cells generate micronuclei, cytosolic DNA and cGAS-STING
activation, making chromosomally unstable tumours immunogenic and therefore responsive to
approved checkpoint inhibitors. **The evidence runs the other way.** In glioblastoma,
chromosomal instability is associated with *paradoxical cGAS upregulation together with
impaired STING signalling*, shaping an immunologically cold microenvironment
(PMID:42540549). High instability does not translate into immunogenicity, so this is recorded
as considered and declined rather than omitted.

## 6c. The symptomatic arm - approved for other diseases, applied to these features

Sections 1-6 ask what acts on the *mechanism*. This section asks a different and more
tractable question: **which medicines already approved for other diseases address the features
this child actually has?** It is the only part of the analysis where approved drugs meet the
proband's own phenotype, and the only part that could affect daily function. Every entry was
resolved to a ChEMBL identifier and its approval status and indication list read from Open
Targets (`data/symptomatic_candidates.tsv`).

This is **symptomatic management, not treatment of MVA** - none of it improves chromosome
segregation. Every item is a question for the treating paediatric nephrologist, oncologist and
clinical dietitian, who hold context this analysis does not.

### Nephrocalcinosis - the strongest position, and our own analysis strengthens it

| drug | ChEMBL | registered indication match |
|---|---|---|
| **potassium citrate** | CHEMBL3989822 | **calcium nephrolithiasis; renal tubular acidosis** |
| hydrochlorothiazide | CHEMBL435 | none - listed indications are cardiovascular |
| chlorthalidone | CHEMBL1055 | none - cardiovascular and renal |

Potassium citrate matches the indication class directly: urinary alkalinisation inhibiting
calcium crystallisation. Thiazides inhibit the sodium-chloride cotransporter, which raises
distal calcium reabsorption and lowers calciuria - the **mechanism** is registered but the
**hypocalciuric indication is not listed**, so that use is off-label in the database's terms
however established it is in paediatric nephrology practice. Stated rather than glossed.

This is the one candidate our own findings make *more* appropriate rather than less. Section 1
concluded that the nephrocalcinosis most likely follows from prematurity and its management
rather than from BUB1B. A feature whose cause lies outside the gene is precisely the feature
that symptomatic management should address - and the prognosis is favourable, since in the
preterm cohort it resolved in 76.5% by three months (PMID:40660023).

### Skeletal muscle atrophy - one approved drug matches, and it collides with our own data

**Givinostat** (CHEMBL1213492) is a histone deacetylase inhibitor approved for, among others,
**Duchenne muscular dystrophy** - oral, with paediatric experience.

It is not uniquely matched to this child's phenotype. **Seven** of the entries in
`data/symptomatic_candidates.tsv` carry a registered indication matching an HPO feature -
givinostat plus potassium citrate, deflazacort, cyproheptadine, megestrol acetate, somatropin
and metreleptin - against two that match on mechanism only (the two thiazides). What
distinguishes givinostat is narrower: it is the only muscle-directed
option here that is **not a glucocorticoid**, so it escapes the immunosuppression and
growth-suppression objections that disqualify deflazacort. It carries a different objection
instead.

It collides with section 5: GSE86928 shows HDAC inhibition arresting the cell cycle and
lowering BUB1B 2.4-fold. Skeletal muscle is post-mitotic so that matters less there, but in a
cancer-predisposed child with a proliferating tumour it is a real concern. Interesting and
requiring caution, not obvious.

A second caveat on transfer: muscle atrophy in MVA is probably neither dystrophic nor
inflammatory, so borrowing from Duchenne rests on similarity of *symptom*, not of mechanism.

**Glucocorticoids rejected.** Deflazacort (CHEMBL1201891) and prednisone are approved in
muscular dystrophies but are immunosuppressive under cancer predisposition and
growth-suppressive in a child who already has growth failure.

### Failure to thrive

Nutritional support first; pharmacology is secondary. **Cyproheptadine** (CHEMBL516) carries
*malnutrition* among its indications and **megestrol acetate** (CHEMBL1201139) carries
*cachexia*. Both approved, both with a weak evidence base in children; megestrol is hormonal
and carries adrenal-suppression risk. Held rather than recommended.

### Short stature - the obvious approved answer is the forbidden one

**Somatropin** (CHEMBL1201621) is approved for short stature in Turner syndrome, growth delay
from IGF deficiency and others - the exact "approved for another disease" pattern. It is
**contraindicated here**: growth-promoting in a cancer-predisposition syndrome with an active
sarcoma. It appears in `data/contraindicated.tsv`, and it is listed precisely because it is
the intuitive move and the wrong one.

### No pharmacological option

Small for gestational age and premature birth are historical findings describing the past, not
states to treat. Lipodystrophy and sarcopenia from the Sieben mouse model have not manifested;
if generalised lipodystrophy ever did, metreleptin (CHEMBL2107857) is approved for it, but
anticipating it serves nothing and the drug carries a black box warning.

### Why this arm carries the scalability argument

Nothing above is BUB1B-specific. The procedure is: take a phenotype description, resolve each
feature to drugs whose registered indications match it, separate indication-level from
mechanism-level matches, and screen the result against the patient's own contraindications.
That generalises to any rare disease with an HPO description, which is what the scalability
criterion asks for.

## 6d. The shortlist, in rank order

Sections 4 to 6c each examine one arm. This section is the consolidated output: every candidate
the search produced, in the order it ranks them. It is generated from
`results/ranked_candidates.tsv` rather than retyped, so the table and the data cannot drift
apart.

Four arms, by what the candidate addresses. **A** acts on the causal mechanism, **B** on a
cellular consequence of aneuploidy, **C** on a clinical feature the child actually has, and
**D** on tumour tissue only. The symptomatic arm carries an explicit **x0.5** multiplier,
because relieving a feature is worth less than acting on the cause; the raw scores are retained
in the table's source so the weighting can be removed and the ordering re-examined.

**Mechanism basis** is the column that matters most for reading the rest. *Registered* means the
proposed action is the one the medicine is licensed for. *Established, unregistered* means the
action is documented in the literature but is not the registered mechanism - which is the premise
of repurposing, not a defect. *Asserted only* means neither, and it is a defect: it applies to
exactly one entry.

| # | candidate | ChEMBL | arm | weighted | mechanism basis | flag |
|---|---|---|---|---|---|---|
| 1 | arimoclomol | `CHEMBL2107726` | A causal | 6.0 | established, unregistered | - |
| 2 | sodium phenylbutyrate | `CHEMBL1746` | A causal | 5.0 | established, unregistered | - |
| 3 | glycerol phenylbutyrate | `CHEMBL2105745` | A causal | 4.0 | established, unregistered | - |
| 4 | acetylcysteine | `CHEMBL600` | B consequence | 4.0 | registered | - |
| 5 | givinostat | `CHEMBL1213492` | C symptomatic | 2.5 | registered | - |
| 6 | potassium citrate | `CHEMBL3989822` | C symptomatic | 2.5 | registered | - |
| 7 | hydrochlorothiazide | `CHEMBL435` | C symptomatic | 2.5 | registered | - |
| 8 | hydroxychloroquine | `CHEMBL1535` | D tumour | 2.0 | established, unregistered | **SYSTEMICALLY CONTRAINDICATED** |
| 9 | ataluren | `CHEMBL256997` | A causal | 2.0 | registered | **APPROVAL IN QUESTION** |
| 10 | chlorthalidone | `CHEMBL1055` | C symptomatic | 2.0 | registered | - |
| 11 | cyproheptadine | `CHEMBL516` | C symptomatic | 1.5 | established, unregistered | - |
| 12 | palbociclib | `CHEMBL189963` | D tumour | 1.0 | established, unregistered | **SYSTEMICALLY CONTRAINDICATED** |
| 13 | megestrol acetate | `CHEMBL1201139` | C symptomatic | 1.0 | registered | - |
| 14 | metformin | `CHEMBL1431` | D tumour | 0.0 | **asserted only** | **CAUTION** |

Three entries carry flags that override their rank. Hydroxychloroquine and palbociclib are
**systemically contraindicated** and belong to tumour-directed use only - their position in this
table is not a recommendation for systemic administration. Ataluren's **approval is in question**:
openFDA returns a single record for it and that record is a bulk ingredient, not a finished
product. Metformin sits last on a weighted score of 0.0 and is the one candidate whose proposed
mechanism is *asserted* rather than documented; it is listed because the search produced it, not
because it is defensible.

Nothing in this table is evidence that any of these medicines helps this child. They are
hypotheses ordered by how much reason there is to examine them, and section 7 names the
measurement that would settle the ones at the top.

---

## 6e. Every axis considered, with its verdict

Sixteen therapeutic axes were examined and each carries a verdict and a reason
(`data/axes_considered.tsv`). This is the record of what was looked at, not of what survived -
the closed and contraindicated rows are the part worth scrutinising.

| axis | lever | verdict |
| --- | --- | --- |
| allele-specific readthrough | read through the premature TGA to make full-length protein | DEGRADED - and the lead compound has no verified product |
| chemical chaperone / proteostasis | raise the amount of the destabilised missense protein that survives quality control | OPEN - best mechanistic fit, and now with an approved product |
| BubR1 abundance via acetylation (NAD+/SIRT2) | block the deacetylation that routes BubR1 to degradation | DEGRADED |
| mTORC1 inhibition (rapalogs) | reduce a damaging secretory/atrophic state | CONTRAINDICATED |
| senolytic clearance | remove cells that have accumulated damage | CONTRAINDICATED |
| autophagy inhibition | remove the buffer that lets aneuploid cells cope | OPEN, CONTESTED |
| metabolic stress on aneuploid cells | exploit the higher metabolic and proteotoxic load aneuploid cells carry | OPEN, CARRIED WITH OBJECTIONS |
| partial Aurora B rebalancing | restore the kinase-phosphatase balance tipped by reduced PP2A-B56 recruitment | OPEN - no approved compound |
| checkpoint reactivation | raise checkpoint strength directly | CLOSED |
| symptomatic management | treat the features the child actually presents with | OPEN - the only arm with approved medicines |
| surveillance, not pharmacology | detect the second tumour early | OPEN - highest expected value |
| measurement, not pharmacology | find out whether BUBR1 abundance is actually reduced | PREREQUISITE |
| arresting division of aneuploid cells | stop aneuploid cells contributing to tissue | CONTRAINDICATED |
| inducing senescence in aneuploid cells | convert aneuploid cells into a stable non-dividing state | CONTRAINDICATED (somatic) / OPEN (tumour, only in sequence) |
| direct pharmacological inhibition of APC/C | supply from outside the restraint on APC/C that the weakened MCC fails to produce | CLOSED - no chemical matter, and aimed at the wrong defect |
| delay separase activation - the cutting step itself | prolong either restraint on separase (securin binding, or cyclin B1-CDK1 binding) so that cohesin is cut later | CLOSED - no chemical matter, and the one druggable target fails on direction |

**One layer was missing from the node set by construction.** The nodes were defined
**structurally** - what the machine is built of - and not **regulatorily**, so nothing that merely
sets the *level* of a component could enter the census. That is exactly where an "up" lever could
live, because a stoichiometric component cannot be up-regulated with an inhibitor, but its
abundance can be raised by blocking its degradation. The acetylation-dependent control of BubR1
abundance is the clearest such lever (PMID:24825348). It is recorded and **degraded**, for three
reasons: a destabilised missense allele is cleared by chaperone-dependent quality control rather
than by that acetylation route, so lever and lesion sit on different pathways for this allele
class; the supporting lifespan result came from transgenic overexpression, not from a drug; and
tumours are NAD+-avid - the entire rationale for NAMPT inhibitors - in a child with an active
sarcoma.

## 7. The test that is missing

Independent of any drug: **premature chromatid separation on the child's own lymphocytes**,
by G-banding with centromere FISH - the assay Wei 2026 (PMID:42434306) applied to carriers.
PCS and mosaic aneuploidy are the cytogenetic hallmark of MVA, this is the cheapest
confirmatory test for the diagnosis, and it appears nowhere in the diagnostic work-up. A
BUBR1 western blot on the same cells would answer whether p.Asn1002Lys actually reduces
protein, and would supply ACMG PS3 - the only realistic route to moving that allele out of
VUS. See `results/validation_plan.md`.

## 8. Generalisation

The structural finding generalises before the pipeline does. The nodes of this machinery are
**stoichiometric components of a complex**: the mechanism needs *more* function, and
small-molecule pharmacology delivers *less*. That holds for any dosage-sensitive complexopathy,
and it predicts the negative result rather than merely recording it.

Nothing in phases 2-6 is BUB1B-specific beyond the curated sign table: given a variant, its
affected module and a phenotype description, the pipeline returns direction-gated candidates
with database identifiers plus an explicit contraindication list. The symptomatic arm, driven
by HPO terms, generalises to any rare disease with a phenotype description. The negative
result is itself portable - "which modules have no approved pharmacology at all" is worth
knowing before a repurposing effort starts.

## 9. Limitations

The *trans* configuration of the two alleles is unproven - no parental data, no phasing, and
10,911 bp is beyond short-read phasing. The functional effect of p.Asn1002Lys is untested.
No missense pathogenicity predictor contributed: SIFT 0.01 and PolyPhen-2 0.997 were returned
incidentally while establishing which calibrated predictors are reachable, and neither
carries an ACMG criterion, because PP3 requires a single predictor chosen in advance with
published calibrated thresholds (PMID:36413997). Protein domain boundaries are UniProt-
verified; the KEN, ABBA and GLEBS positions are from review literature and are flagged as
unverified in `data/domain_defect_map.tsv`. Drug warning text was unavailable from Open
Targets - only the warning type. Approval is checked at two levels: the target database for the substance, and openFDA
`/drug/ndc` for the US product registration category, where 13 of 14 candidates resolve to a
genuine NDA or ANDA (`data/approval_status_openfda.tsv`). **EU product status is unverified for
every row**, because the EMA publishes no comparable API - and that is exactly the gap that
decides ataluren, which has no US product and an EU authorisation reportedly not renewed. The
two levels also disagree for arimoclomol, recorded at PHASE_3 in the target database while
openFDA returns its NDA directly; where they conflict, the regulator is treated as the
authority.

## 10. How this was produced

**Generative AI.** Anthropic Claude, Team plan subscription - commercial terms, customer content
is not used to train models - model `claude-opus-5`, used through Claude Science. One model did
all of the work: no secondary or utility model contributed, and no sub-agents were delegated. The
full methods description, including the boundary between automated and curated steps, is
Appendix A below.

`results/ranked_candidates.tsv` is written by no script. The scores are assigned by hand and the
pipeline *validates* them: `validate_ranking.py` checks the arithmetic, the ordering, the sourcing
of every mechanism claim and the consistency of every flag, and `verify_approval.py` re-checks
approval row by row. Automation established what is available and in which direction it acts;
the ordering is ours, and it is checkable.

**Four rules constrained the search.**

1. **No drug name or mechanism from memory.** Every candidate carries a ChEMBL identifier
   resolved against its source; every mechanistic claim a PMID resolved through NCBI. The rule
   exists because a language model will produce a *plausible* identifier for a claim it has
   invented, and a plausible identifier is indistinguishable from a real one until it is looked up.
2. **Only medicines already in human use.** Approval in any one jurisdiction, for any indication,
   qualifies - but it is verified at product level, not substance level.
3. **Absence is recorded as a row of zeros, never as silence.** A census, not a search.
4. **Direction of action is a gate, not a scoring weight.** A wrong-direction drug is removed with
   its reason recorded rather than down-ranked.

**Corrections made during the work.** Recorded because a method is partly the list of its own
repairs.

- **"No mechanism recorded" was imprecise.** It meant "no row in Open Targets", but read as "this
  drug has no known mechanism". Querying ChEMBL directly showed Open Targets under-reports:
  sodium phenylbutyrate returns no mechanism row there but one in ChEMBL (glutamine sequestering),
  acetylcysteine none and two. Every affected field was rewritten.
- **"Mismatch" was framed as a defect.** Proposing a drug outside its registered indication *is*
  repurposing; treating that as a finding was an error. The check was restructured into three
  evidential levels - registered, established but unregistered, asserted - and only the last is a
  defect. It is what removed fingolimod and perphenazine.
- **Approval was being read at substance level only**, which cannot distinguish a medicine from a
  supplement. openFDA closed this for the US; EU product status remains unverified.
- **The reproducibility claim was not true until it was tested end to end.** The first full online
  check showed three tables differing from their committed copies, and neither cause was data
  drift. *Row order* was inherited from API response order and dict iteration, so tables agreeing
  on every value still compared unequal - `fetch_candidates.py` and `gate_direction.py` now sort
  before writing. And two columns of `candidates_raw.tsv` (`drug_type`, `withdrawn`) **cannot** be
  regenerated, because the public v4 `Drug` type exposes neither, so they are carried as curated
  values named in `CURATED_COLUMNS` rather than silently compared. The lesson generalises: a
  reproducibility check that trips on row order is one a reader quickly learns to ignore, which is
  worse than having no check.
- **An "only node requiring down-regulation" claim was an overstatement.** Three nodes require it;
  APC/C is only the one where the lesion itself is the missing restraint.

**What this repository holds of the patient.** The variant pair and the clinical features, both
of which the organisers require the report to state, and nothing else: no genotype tables, no
sample identifiers, no files derived from the VCF, and no ontology term identifiers - the
clinical feature is named in prose instead. Verified by search across every markdown, table and
script here.

---

## Appendix A. Methods description

Answers to the organisers' methods-description form. Question wording is the organisers'; the
answers are ours. All questions are voluntary except the first and the generative-AI disclosure.

### Team name

**Bebok**

---

### If commercially-available Generative AI was used, record the provider, the plan or tier, and the relevant setting in your methods description.

Anthropic Claude, Team plan subscription - commercial terms, customer content is not used to
train models - model `claude-opus-5`, used through Claude Science. One model did all of the
work: no secondary or utility model contributed and no sub-agents were delegated.

**Skills loaded: one.** `figure-style`, used for the schematic in `docs/track2_flow.png` -
publication-figure conventions, not analysis. The analytical work did not run through domain
skills; it runs through the standard-library scripts in `src/`, which call public APIs directly.

**How the queries are issued.** The scripts call the public GraphQL and REST endpoints directly,
against named URLs, and record the data release beside every count
(`meta{apiVersion dataVersion}` = release 26.09, API 26.9.0), so the same query can be reissued
years later and the drift will be visible. They are standard-library-only with PEP 723 headers,
so `python3 src/<script>.py` runs with no installation, no virtual environment, no API key and no
account.

In summary: generative AI was used for database querying, pipeline code, literature retrieval and
drafting. Exactly one procedural skill was loaded, and it concerned figure conventions, not
analysis; the analytical work runs through standard-library scripts calling public APIs directly,
with the data release recorded beside every count, so no result depends on a client or account
that a reader cannot pin.

It was **not** permitted to supply drug names or mechanisms: every candidate carries a ChEMBL
identifier resolved against Open Targets, and every mechanistic claim a PMID resolved against
PubMed. That rule exists because a language model will produce a plausible identifier for a
claim it has invented, and a plausible identifier is indistinguishable from a real one until it
is looked up. `src/verify_citations.py` re-resolves every identifier in the repository.

---

### Please describe your approach in detail (how you went from variant/mechanism to candidate medication(s)).

**Starting point.** A nonsense allele and a rare missense allele in *BUB1B*. Complete BUBR1 loss
is embryonic-lethal and the child is alive, so the state is a reduced dose of partly functional
protein. The therapeutic target is therefore to cross a dosage threshold, not to reconstitute
wild-type protein - and because a germline hypomorph cannot be corrected in every cell, the
operative question became **"which layer of consequence is reachable"** rather than "which drug
fixes *BUB1B*".

**Step 1 - census, not search.** Every node of the broken machinery was queried and a node with
no drug produced a row saying zero rather than no row at all: 65 spindle-checkpoint and
kinetochore nodes, of which **60 carry no drug of any kind**. The five exceptions hold 49
node-compound pairs (36 distinct molecules): 47 inhibitors, one antisense inhibitor, one vaccine
antigen; phases 1 through 3; **none approved**. The cohesion layer was counted separately,
because premature chromatid separation is the cytogenetic hallmark: **17 targets, 16 with no
compound at all**, separase included. Asked of the disease rather than the gene, nine MONDO
entries across the MVA group return **no medicine at any clinical stage**.

**Step 2 - the structural reading.** These are stoichiometric components of a complex: the
mechanism needs *more* function and small-molecule pharmacology delivers *less*. The obstacle is
also one of coverage, not only direction - three nodes do require the direction to be down
somatically, and at the one where the lesion itself *is* the missing restraint (APC/C), all
seven subunits return zero compounds.

**Step 3 - signed network and a direction gate.** STRING supplies topology but is undirected and
unsigned, so the sign was curated by hand, one PMID per row, with separate columns for what
somatic tissue needs and what the tumour needs. Direction is then applied as a **gate, not a
scoring dimension**: a compound acting on the right protein in the wrong direction is removed
with its reason recorded, because direction treated as one weight among several lets a
well-documented drug with the wrong sign float to the top. 67 rows passed, 31 were rejected,
**none of the passes is approved**. An acceptance test is built in and fails the run if
violated: MPS1/TTK inhibitors must fail for somatic tissue and pass for the tumour.

**Step 4 - consequence arm from the literature, verified against databases.** Because no drug
acts on a network node, remaining hypotheses are literature-derived and were checked against
registered mechanisms. The distinction that matters is NOT registered-versus-not, because proposing a drug outside its
registered indication is what repurposing *is*. It is where the mechanistic claim's evidence sits,
and that splits three ways rather than two: the registered mechanism is itself the lever invoked;
the activity is established pharmacology but is not the registered mechanism, which is ordinary
and legitimate repurposing reasoning; or the activity is asserted for that drug without being part
of its established pharmacology. Only the third is a defect. It is also what the check caught:
fingolimod and perphenazine were both proposed as PP2A activators on review-literature grounds
that neither drug's established pharmacology supports, and both were removed.

> **Caveat on this figure.** The 7-of-8 count comes from the target database's mechanism field, which is NOT a complete mirror of ChEMBL. Demonstrated: sodium phenylbutyrate returns no mechanism row in Open Targets but one in ChEMBL (`SEQUESTERING AGENT: glutamine sequestering agent`); acetylcysteine returns none and two respectively. The count therefore **may overstate mismatch**, since an absent row sits in the same bucket as a genuinely different mechanism. Re-deriving all eight directly against the ChEMBL mechanism endpoint was attempted and did not complete (intermittent HTTP 500), so this is reported as an Open-Targets-derived count. For the candidate that matters most it is unaffected: phenylbutyrate's registered mechanism really is nitrogen scavenging, not chaperoning.


**Step 5 - four arms and two rankings.** Candidates are grouped by what the treatment addresses -
causal mechanism, cellular consequence, symptomatic, tumour-directed - and ranked both within
arm and globally. The global score is
`(benefit ceiling + evidence quality + executability - risk) x category weight`, with the
symptomatic arm at 0.5 because managing a presenting feature answers a different question from
acting on the disease. Components sit beside the total and the raw score is kept, so a reader
who disagrees with a weight can re-rank without recomputing. These are ordering weights, not
calibrated probabilities.

**Step 6 - what must not be given**, as a separate deliverable: eleven entries, each with the
triggering clinical feature.

### Was candidate identification automated (e.g., via a drug-target database query or LLM-assisted search), or primarily manual literature review and expert curation?

**Both, in a deliberate division of labour, and the boundary is recorded per table.**

*Automated:* the census of every node, the compound retrieval with action type and clinical
stage, the direction gate and its self-test, the ranking arithmetic, and the resolution of every
PMID and ChEMBL identifier. These run from `src/` against the public Open Targets GraphQL
endpoint and NCBI E-utilities, and `--mode check` regenerates the tables into a temporary
directory and diffs them against what is committed. The node census reproduces byte-identically.

*Manual:* the **sign** of each network edge (one PMID per row, because the available topology is
unsigned), the scoring components, the arm assignment, and every clinical-safety judgement.

*Hybrid:* candidate generation in the consequence arm began as LLM-assisted literature search
and was then **filtered by database verification**, which rejected most of what it proposed.

The honest summary: **automation decided what is available and in which direction; humans decided
what is advisable.**

### Please describe any manual literature review or expert curation in detail.

Three places where curation was unavoidable.

**Network sign.** STRING edges are undirected and unsigned, so for each node we recorded the
wild-type role, the direction the defect pushes it, and the direction a drug would need to push
it in each tissue - one citation per row. Coverage is **11 of 65 nodes**, and that is stated as
a limitation: this is the component that does not scale, and replacing it with a signed,
directed regulatory resource is the single highest-value engineering change available.

**Mechanistic claims about the protein.** Domain boundaries were taken from UniProt O60566 and
flagged per row as verified or literature-only - the KEN, ABBA and GLEBS positions commonly
quoted in reviews are **not** annotated in UniProt and are marked unverified. The pseudokinase
reading matters for candidate selection: the domain containing the missense residue is
catalytically dead in humans, so the tractable question is protein stability rather than kinase
activity, which is why levers that raise abundance dominate arm A.

**Sequence context read directly from the reference.** For the readthrough hypothesis we read
the stop-codon context ourselves rather than citing it: the premature stop is `TGA`, the most
permissive codon, but the +4 base is `A`, the least permissive position in the C>U>G>A
hierarchy, so local context works *against* readthrough. Near-cognate decoding would insert Trp,
Arg or Cys rather than Leu, so any full-length product is not wild-type BUBR1.

**Sixteen therapeutic axes** carry a written verdict and reason, including the 11 that do not
produce a candidate. That table is the record of what was examined rather than of what survived.

### Please state whether you relied solely on publicly available drug/target databases and literature (preferred), or whether proprietary data sources were also used.

**Publicly available sources only.** No proprietary database, no commercial drug-intelligence
subscription, no in-house compound collection.

Approval is checked at two levels. **Substance level** comes from the target database, which
reports a jurisdiction-agnostic "approved somewhere, for something" status - useful but not
sufficient, as one NAD+ precursor demonstrates by returning approved status while being marketed
as a supplement. **Product level** comes from openFDA `/drug/ndc`, which returns the actual US
registration category per product: 13 of 14 candidates resolve to a genuine NDA or ANDA,
recorded per drug in `data/approval_status_openfda.tsv`.

The qualifying bar is stated so it cannot be read narrowly: **approval in any one jurisdiction,
for any indication, qualifies.** Repurposing needs a medicine with a known safety profile and an
existing supply chain; the indication it was registered for is not the point, and none of these
medicines is registered for this disease.

What remains open, and it is not cosmetic: **EU product status is unverified for every row**,
because the EMA publishes no comparable API. That gap decides the single flagged candidate -
ataluren returns one openFDA record, `BULK INGREDIENT`, so it has no US product, and its EU
conditional authorisation was reportedly not renewed in March 2025. It is retained with an
explicit flag and zero executability rather than deleted, because the EU status needs checking
rather than assuming.

### Please describe the public data sources used in detail (e.g., DrugBank, ChEMBL, ClinicalTrials.gov, OMIM, PubMed literature search).

| source | used for |
| --- | --- |
| Open Targets Platform (release 26.09), public GraphQL endpoint | drug-target associations, action type, mechanism of action, maximum clinical stage, registered indications, disease-to-target association, variant effect predictions, tractability |
| ChEMBL (via Open Targets identifiers) | compound identity; every candidate carries its ChEMBL id |
| STRING | network topology - undirected and unsigned, which is why the sign is curated |
| UniProt (O60566) | BUBR1 domain boundaries, phosphosites |
| MONDO / OMIM / Orphanet | the MVA disease group and its subtypes, including the grouped entry |
| NCBI E-utilities (PubMed) | resolution of every cited PMID |
| GEO (GSE86928, GSE134781) | expression evidence, including a result retracted on inspection |
| Ensembl REST | reference sequence for the stop-codon context |
| **ChEMBL** mechanism endpoint (direct) | registered mechanism of action per molecule, queried directly because the target database does not mirror it completely - sodium phenylbutyrate returns no mechanism row in Open Targets but one in ChEMBL (glutamine sequestering), and acetylcysteine returns none and two respectively |
| **openFDA** `/drug/ndc` and `/drug/label` | **US product-level registration category** (NDA / ANDA / BLA / bulk ingredient), marketing start date and the verbatim label indication |
| ClinicalTrials.gov | checked for MVA and *BUB1B*; returns no registered trials |

**One substitution is declared:** the Reactome REST API returned a Cloudflare block to our
client, so pathway membership was taken from the project's own five-layer catalogue and
validated gene by gene against Open Targets - 65 of 78 candidate symbols confirmed, four
mis-resolved tokens rejected.

### Please describe any proprietary data sources used in detail.

None.

### How did you characterize the variant's mechanism (e.g., loss-of-function vs. gain-of-function, pathway disrupted, downstream biological consequence)?

**Direction of effect: loss of function, specifically a hypomorph acting through reduced protein
stability.** Three lines point that way, and we state plainly that they are **not three
independent lines** - an earlier draft claimed they were. Recessive inheritance with a
high-confidence null in *trans* is an expectation that presupposes the diagnosis in order to
establish the allele's mechanism. Position inside the pseudokinase domain, which is
catalytically dead in humans and contributes stability rather than activity, and a destabilising
folding-energy estimate are the **same observation at two levels of description**, since the
estimate is computed on the structure whose domain assignment is the other line. Read as two
partly independent lines plus an inheritance-based expectation.

Magnitude is reported as **direction only**. The characterised neighbouring allele scores
considerably more destabilising, but quoting a ratio over-interprets a non-linear scale with a
method error of roughly +/-0.5 kcal/mol on a substitution the method penalises predictably. No
calibrated missense predictor contributed, and the folding-energy tool is not ClinGen-calibrated,
so the ACMG classification of the missense allele remains **uncertain significance** and the
mechanistic reading improves the hypothesis, not the classification.

**Pathway disrupted.** BUBR1 has two functions and both are impaired: it is a component of the
mitotic checkpoint complex that restrains APC/C-CDC20, and its KARD recruits PP2A-B56 to oppose
Aurora B. Checkpoint failure alone accounts for at most half the phenotype.

**Downstream consequence.** Reduced BubR1 dose weakens MCC restraint of APC/C, so cyclin B leaks
away faster and the mitotic arrest shortens - more mitotic slippage, more aneuploid and
tetraploid daughters. Two phenomena are kept apart deliberately: premature anaphase **onset** is
a timing defect, whereas premature chromatid **separation** - the cytogenetic hallmark - is a
separation defect in which sisters part before anaphase. That distinction is load-bearing: it is
why direct APC/C inhibition, which addresses timing, may be aimed at the wrong step.

**What would settle it.** A BUBR1 western blot on the child's cells. The abundance premise under
every candidate in arm A is currently unmeasured, and if abundance is normal that premise fails
and we would report it.

### Please provide an estimate of time or effort spent on this analysis (if possible).

Roughly one analyst-week, dominated by verification rather than by generation: drafting
candidates is fast, checking each against registered mechanism, clinical stage, indication list
and primary literature is not. The computational cost is negligible - the full pipeline is a few
hundred public API calls and runs in minutes on a laptop with no paid data source. The
rate-limiting resource was clinical and mechanistic judgement, which is also where the method is
least automatable.

### Provide a method abstract (up to 500 words). Please include strengths and limitations of your approach as well as methodological detail.

We treated the central finding as a result rather than as a gap. Querying all 65 spindle-checkpoint
and kinetochore nodes returns 60 with no drug of any kind; the five exceptions carry 49 node-compound pairs, none approved and none
restorative. The cohesion layer gives 17 targets, 16 with no compound,
separase included. Nine MONDO entries across the MVA group return no medicine at any stage. The reason is structural: these are
stoichiometric components of a complex, so the mechanism needs more function while
small-molecule pharmacology delivers less. The obstacle is also one of coverage - at the one
node where the lesion itself is the missing restraint, all seven subunits return zero compounds.

Method. A census rather than a search, so absence produces a row of zeros instead of silence. A
hand-signed interaction network, because available topology is undirected and unsigned, with
separate required directions for somatic tissue and tumour. Direction applied as a gate that
removes with a recorded reason rather than as a scoring weight, with an acceptance test that fails
the run if MPS1 inhibitors misbehave. Candidates grouped into four
arms by what the treatment addresses, ranked within arm and globally, with the symptomatic arm
down-weighted because managing a feature answers a different question from acting on the
disease. Every identifier re-resolved against its source; every row carrying the strongest
argument against itself.

Strengths. The negative result is reproducible in one command and the census regenerates
byte-identically against a named release - a negative claim rots quietly when the database
beneath it moves. Direction as a filter structurally prevents the
failure mode where a well-documented drug with the wrong sign ranks first. Mechanistic claims were sorted by evidential basis -
registered, established, or asserted - and the two resting on assertion were removed. Sixteen therapeutic axes carry written verdicts including the 11 that
yield no medicine, so the panel sees what was examined rather than what survived.

Limitations. Phase is unproven, so the *trans* configuration is inferred. The missense allele's
functional effect is unmeasured, and the hypomorph reading rests on two partly independent lines
plus an inheritance-based expectation, not three independent lines. Folding-energy estimates
carry direction, not magnitude. Network sign covers 11 of 65 nodes and is hand-curated, which is
the component that does not scale. US product-level approval is verified against openFDA; EU product
status is unverified, which is the gap deciding the one flagged candidate. Safety screens were
applied to the consequence arm but not to the network arm's somatic passes. The top candidate is registered only in combination with another
medicine; the next on that arm carries unfavourable evidence from our own expression analysis.

The conclusion we would defend is that the highest-value next step is a measurement rather than
a medicine: a BUBR1 western blot, a premature-chromatid-separation rate, and tumour surveillance
ahead of any drug on our list.

---
