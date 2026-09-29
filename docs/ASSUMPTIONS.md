# Assumptions and limitations

Every item below is a decision this pipeline makes on your behalf or a
limit it cannot overcome. The purpose of writing them down is that a
reader can decide whether a conclusion survives them. If you fork this
for a real study, edit this file rather than deleting it.

## 1. Study design

**Repeated samples per donor.** The default analysis assumes each donor
contributes a baseline and a post-treatment sample, so donor is modelled
as a random intercept. If your design is cross-sectional, pass
`--random-effect none` and be aware that the treatment estimate is then
confounded with between-person variation, which in gut data is large.
In the demonstration run donor accounted for a median 84 percent of CLR
variance per taxon against 0.7 percent for the treatment.

**Two groups only.** Both implementations here handle a two-level factor.
Multi-group and dose-response designs need the R `ANCOMBC::ancombc2`
global test, which is available through `R/03_differential_abundance.R`
but not through the Python package.

**Batch is not confounded with the arm.** If every treated sample was
extracted in one batch and every control in another, no model separates
them. `scripts/03_ordination.py` reports batch variance so you can check
this before spending effort on testing.

## 2. Sequencing and processing

**Single sequencing run for the error model.** DADA2 learns a
run-specific error model. Samples from multiple runs must be processed
separately and their sequence tables merged. Pooling runs into one error
model degrades variant resolution.

**truncLen must leave overlap.** For V4 with 515F/806R the merged
amplicon is about 253 bp. Truncation lengths that sum to less than the
amplicon plus roughly 20 bp of overlap cause silent failure at the merge
step. Check `read_tracking.tsv`; a large drop at `merged` is a parameter
error, not biology.

**Primers are assumed removed.** If primers remain in the reads DADA2
will treat them as biological sequence. Remove them with cutadapt first.

**Chimera removal is de novo consensus.** This can over-remove in samples
with very high diversity and under-remove when a chimera is abundant in
most samples.

## 3. Taxonomy and reference databases

**Species assignment is unreliable from 16S.** The gene does not resolve
most species. Only exact-match species assignment is used, and even an
exact unique match does not identify a unique organism, because distinct
species share identical V4 sequences.

**SILVA and GTDB are not interchangeable.** SILVA is curated on rRNA
phylogeny. GTDB is built on whole-genome relatedness and has renamed and
split many genera. A genus label from one database does not denote the
same clade as the same label from the other. Do not merge results across
databases and do name the release in every figure caption.

**Databases are biased towards cultured and clinically studied
organisms.** Unassigned means unrepresented, not absent, and the
unassigned fraction varies systematically by body site and host species.
Report the unassigned percentage at each rank; the script prints it.

**Classification is better against a region-trimmed reference.** If the
primer pair is known, extract the amplified region from the reference
before training.

**Bootstrap values are not probabilities.** The `minBoot` threshold of 80
measures stability of the assignment to k-mer resampling, not the
probability that the organism belongs to that taxon.

**Do not aggregate before testing.** Collapsing ASVs to genus can merge
taxa with opposite responses and cancel real effects. Test at ASV level,
aggregate for presentation.

## 4. Compositionality

**Absolute abundance is not recoverable.** This is the hardest limit.
A uniform change in total bacterial load is mathematically invisible to
every method here. The simulator can create exactly this situation with
`--biomass-ratio`, and the pipeline output looks identical whether or not
the load changed. Resolving it requires external quantification: flow
cytometry, qPCR of total 16S copies, or a synthetic spike-in. If you do
not have one, phrase every result as relative.

**CLR uses a moving reference.** The geometric mean of the sample is the
reference, so when a large share of the community shifts one way,
unchanged taxa acquire opposite-sign artefacts. This was measured in the
demonstration run: the mean CLR coefficient across genuinely null taxa
was −0.164 when 30 of 300 taxa rose together. Use `--denom iqlr` in
`scripts/04_differential_abundance.py`, or rely on the ANCOM-BC2 bias
correction, which reduced the mean null effect from −0.337 to −0.037.

**16S copy number varies between one and more than fifteen per genome.**
Read counts are therefore not proportional to cell counts even within a
single sample. No copy-number correction is applied, because correction
factors are themselves model-dependent and database-dependent. Results
describe sequence abundance.

**Zero handling is a choice with consequences.** ALDEx2 avoids a point
replacement by sampling a Dirichlet posterior. ANCOM-BC2 uses a
pseudocount and then tests whether the answer depends on it. The
`passed_ss` column is the output of that test and features failing it
should not be reported. In the demonstration run 128 of 299 features
failed.

**Rarefying is not used and is not recommended.** Discarding reads to
equalise library sizes throws away real information and does not address
compositionality, which is a property of the data-generating process
rather than of unequal depth.

## 5. Contamination

**Contaminant identification requires negative controls or DNA
quantitation.** With neither, `scripts/02_decontam.py` exits and says so
rather than substituting a prevalence filter, which does not identify
contaminants.

**Contaminant profiles are batch-specific.** Kit lots differ, so scores
are computed within batch and combined when a batch column exists.

**Low-biomass studies are the dangerous case.** In skin, lung, placenta,
or a heavily antibiotic-treated gut, reagent DNA can dominate. Because
contaminant load is roughly fixed in absolute terms, its relative
abundance rises when real biomass falls, so contaminants appear enriched
by exactly the treatments that reduce biomass. Any perturbation study
that could change biomass must run this step.

**The threshold of 0.1 is permissive on purpose.** Retaining a
contaminant in a differential abundance test is more costly than losing a
genuine rare taxon. Report the threshold and the number removed.

## 6. Statistical inference

**Benjamini-Hochberg assumes positive regression dependence.** Amplicon
features are strongly dependent through the constant-sum constraint.
BH is the field convention and is used here, but the nominal level should
be read as approximate. The simulation study measures the realised false
discovery proportion rather than assuming the nominal one holds.

**ANCOM-BC2 assumes most taxa are not differentially abundant.** The bias
estimate is the location of the bulk of the per-taxon effect distribution.
If a perturbation moves most of the community one way, the correction
absorbs real signal and the method under-reports. The simulator's
`asymmetric` flag stresses this deliberately.

**The bias estimate carries uncertainty.** It is propagated into each
taxon's standard error here. Implementations that ignore it are
anti-conservative.

**A random effect does not correct compositional bias.** These are
separate problems. The demonstration run shows a mixed model on CLR values
returning a realised false discovery proportion of 0.383 while ANCOM-BC2
on the same data returns 0.

**What the random effect buys depends on the design.** In a paired or
crossover layout, donor is crossed with the arm and blocking on it
increases power; ignoring donor loses findings rather than inventing
them. In a nested layout, where each donor sits in one arm and
contributes several samples, ignoring donor is pseudoreplication and
inflates significance. In the nested simulation, ignoring donor gave a
realised false discovery proportion of 0.50. Do not cite the paired
argument to justify a nested design or the reverse.

**Mixed model convergence is not guaranteed.** 270 of 299 taxa converged
in the demonstration run. Non-convergence falls back to a fixed-effect fit
and is flagged in the `converged` column. Do not silently treat a fallback
fit as a mixed one.

**Realised FDP in one dataset is a noisy draw.** Single-dataset
performance numbers should not be read as calibration. Use
`scripts/06_simulation_study.py` with several replicates.

## 7. The simulator

The generator models overdispersion, sparsity, donor random effects,
batch effects, variable library size, reagent contamination, optional
biomass change, and the multinomial sampling that creates
compositionality. It does not model phylogenetic correlation between
taxa, non-linear dose response, zero-inflation beyond sampling zeros,
temporal autocorrelation, or body-site-specific community structure.
Performance on simulated data is a necessary condition for trusting a
method, not a sufficient one.

## 8. Reproducibility

All random seeds are explicit and exposed as arguments. The R taxonomy
script records the reference path, md5 and package versions.
`sessionInfo()` is written alongside the R results. Manual metadata
corrections must be logged in `data/metadata/curation_log.md`, because
public sample metadata is frequently incomplete or inconsistently
labelled and an uncorrected sample sheet is the most common silent error
in reanalysis of public data.
