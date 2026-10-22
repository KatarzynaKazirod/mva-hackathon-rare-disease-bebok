# Validation plan - what would refute each hypothesis

Every candidate below is matched to a model system from `data/model_systems.tsv` and a
readout already published for that system. A hypothesis with no refuting experiment is not
on the shortlist.

These are research proposals for expert review, not treatment recommendations.

## 1. Sodium phenylbutyrate - does it raise BUBR1?

**Hypothesis.** p.Asn1002Lys destabilises the protein; a chemical chaperone raises the
steady-state level of the hypomorphic product past the threshold at which cells function.

**Why a partial effect suffices.** Matsuura 2006 (PMID:16411201) put the pathogenic
threshold at a **>50% reduction** in BUBR1. The goal is not full correction but crossing
back over that threshold - a much lower bar than restoring wild-type protein.

| | |
|---|---|
| model | patient lymphoblastoid line, or a cell line engineered to carry p.Asn1002Lys |
| readout | BUBR1 western blot, +/- drug |
| prediction | measurable increase in total BUBR1 |
| refuted if | BUBR1 unchanged |
| secondary readout | premature chromatid separation rate falls (G-banding + centromere FISH) |
| timeframe | weeks |

**The honest weakness.** Open Targets records **no mechanism of action** for this drug. The
chaperone activity is literature-only, and it has never been tested in any
chromosomal-instability system. This experiment is the whole evidence base.

## 2. Hydroxychloroquine - tumour arm only

**Hypothesis.** Aneuploid cells carry proteotoxic stress and depend on autophagy to survive
it, so blocking autophagy is selectively toxic to them.

**The direction conflict that confines it to the tumour.** Somatic tissue in this child is
already aneuploid and needs autophagy to tolerate that. Blocking autophagy systemically acts
in the wrong direction for the patient's own cells. The hypothesis is therefore about the
rhabdomyosarcoma, not about the syndrome.

| | |
|---|---|
| model | patient-derived aneuploid cells vs euploid control |
| readout | differential viability |
| prediction | selective loss of the aneuploid population |
| refuted if | equal toxicity in both, i.e. no selectivity window |

**Weakness.** The registered mechanism is TLR7/TLR9 antagonism, not autophagy inhibition.
The aneuploidy-stress rationale is generic and not BUB1B-specific.

## 3. Ataluren - genotype-matched but blocked

**Hypothesis.** p.Leu737Ter is `TTA` -> `TGA`, and TGA is the most readthrough-permissive of
the three stop codons, so a readthrough agent should restore full-length protein from that
allele.

**Why it cannot work alone.** The transcript is NMD-competent: the stop sits in exon 17 of
23, 8,861 bp upstream of the last exon-exon junction, so the mRNA is degraded before it can
be read through. Readthrough requires a partner that inhibits NMD, and **no approved NMD
inhibitor exists.** The hypothesis is sound and currently unexecutable.

| | |
|---|---|
| model | patient lymphoblastoid line |
| readout | BUB1B transcript by RT-PCR, then BUBR1 western |
| prediction | neither changes on readthrough agent alone |
| would become testable if | an NMD inhibitor reaches the clinic |

Record also that the obvious approved readthrough route, aminoglycosides, is excluded here
on nephrotoxicity grounds given the nephrocalcinosis - see `data/contraindicated.tsv`.

## 4. The confirmatory test this child does not yet have

Independent of any drug: **premature chromatid separation on the child's own lymphocytes**,
by G-banding with centromere FISH, against healthy controls - the assay Wei 2026
(PMID:42434306) used in carriers. PCS and mosaic aneuploidy are the cytogenetic hallmark of
MVA, this is the cheapest confirmatory test that exists for the diagnosis, and it appears
nowhere in the Track 1 work-up. It should be proposed before any therapeutic hypothesis.

A BUBR1 western blot on the same cells answers the second open question from Track 1 -
whether p.Asn1002Lys actually reduces protein - and would supply PS3 at strong level, which
is the only realistic route to moving that allele out of VUS.
