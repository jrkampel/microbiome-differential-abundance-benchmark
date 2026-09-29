# mdab: differential abundance methods for 16S data, compared honestly

A compositional analysis pipeline for 16S amplicon data that runs ALDEx2
and ANCOM-BC2 on the same dataset, models donor and batch as random
effects, screens for reagent contaminants with decontam, and then does
the thing most microbiome pipelines do not: measures whether the methods
are right, using simulated data with known ground truth.

The claim this repository makes is narrow and defensible. It does not
claim to identify the correct method. It shows where two accepted methods
disagree, why they disagree, and what each disagreement costs in false
discoveries or lost power under conditions you can vary.

---

## Why compositionality is the whole problem

A sequencer returns a roughly fixed number of reads per sample. That
number is set by library loading and flow cell capacity, not by how much
bacterial DNA was in the specimen. A sample is therefore a fixed-size
random draw from an unknown pool, and the data carry only relative
information. If one taxon truly blooms tenfold, the read counts of every
other taxon fall even though nothing happened to them.

The centred log-ratio removes the dependence on the total:

```
clr(x_ij) = log(x_ij) - (1/D) * sum_k log(x_kj)
```

Differences of CLR values are log ratios, and a log ratio does not change
when the library size changes. Two costs come with this. CLR values sum
to zero within a sample, so the feature covariance is singular. And the
reference is the sample's own geometric mean, so when a large share of
the community moves in one direction the reference moves with it, and
unchanged taxa acquire apparent changes of the opposite sign. This
repository demonstrates that artefact rather than asserting it. See
"Result 4" below.

Zeros break the logarithm, so zero handling is a modelling decision. The
two methods compared here handle it in opposite ways, and that is most of
why they disagree.

### ALDEx2

ALDEx2 refuses to choose a zero replacement. It treats each sample's true
composition as unknown and draws Monte Carlo instances from the posterior
under a Jeffreys prior:

```
p_j | n_j ~ Dirichlet(n_j + 1/2)
```

Sparse features get wide posteriors, so uncertainty from low counts
propagates into the test instead of being fixed beforehand. Each instance
is CLR-transformed and tested, and the reported q value is the expectation
of the per-instance Benjamini-Hochberg value. Its effect size is the
median between-group difference divided by the larger within-group
dispersion, a standardised quantity rather than a fold change.

The consequence is that ALDEx2 is conservative on rare taxa, sometimes
severely so.

### ANCOM-BC2

ANCOM-BC2 models the observed count as the true absolute abundance times
an unknown per-sample sampling fraction:

```
log O_ij = log A_ij + d_j + e_ij
```

The term `d_j` is not identifiable from one sample, but its group-level
difference is estimable if most taxa are not differentially abundant. Fit
a per-taxon linear model, take the bulk of the distribution of estimated
group effects, and that location is the bias. Subtract it and run ordinary
inference.

That assumption is load-bearing and it is stated here because it can
fail. If the perturbation genuinely moves most of the community one way,
the bias estimate absorbs part of the real signal. An antibiotic arm that
collapses a phylum is exactly such a case.

Because the correction leaves an ordinary linear model behind, ANCOM-BC2
accepts covariates and a random intercept, which is why it can handle
repeated donor sampling directly. It also runs a pseudocount sensitivity
analysis, refitting across plausible pseudocounts and flagging taxa whose
significance depends on that arbitrary choice. In the demonstration run
below, 128 of 299 features failed that check.

---

## What this repository contains

```
src/mdab/            readable Python implementations of the statistics
  compositions.py    closure, zero replacement, CLR, ALR, IQLR, Aitchison
  aldex2.py          Dirichlet Monte Carlo CLR with Welch and Wilcoxon
  ancombc.py         bias-corrected LM/LMM, structural zeros, sensitivity
  decontam.py        prevalence, frequency and combined contaminant scores
  mixed.py           per-taxon mixed models on CLR, ICC, variance partition
  simulate.py        generator with donors, batch, contaminants, truth
  compare.py         concordance, agreement, scoring against truth
R/                   the reference implementations, authoritative on real data
  01_dada2.R         reads to ASVs
  02_assign_taxonomy.R  SILVA or GTDB, with provenance recording
  03_differential_abundance.R  decontam, ALDEx2, ANCOM-BC2, concordance
scripts/             the runnable pipeline
  00_fetch_public_data.py   ENA portal download and metadata triage
  01_simulate_dataset.py    ground-truth dataset when ENA is unreachable
  02_decontam.py            contaminant screen
  03_ordination.py          Aitchison PCA, blocked PERMANOVA, variance partition
  04_differential_abundance.py  the method comparison
  05_random_effects.py      mixed model against pseudoreplication
  06_simulation_study.py    repeated simulations, FDR calibration and power
tests/               invariance and behaviour tests, 15 of them
docs/ASSUMPTIONS.md  every assumption and limitation, in one place
```

Two implementations of the same statistics exist on purpose. The R
scripts call the packages the field uses and are authoritative for real
results. The Python package exists so the logic can be read, so the
pipeline runs without Bioconductor, and so the simulation study can be
driven from one process.

---

## Quick start

```bash
pip install -r requirements.txt

# ground-truth dataset, no network needed
python scripts/01_simulate_dataset.py

python scripts/02_decontam.py
python scripts/03_ordination.py
python scripts/04_differential_abundance.py
python scripts/05_random_effects.py
python scripts/06_simulation_study.py --replicates 4 --donors 30

PYTHONPATH=src python -m pytest tests -q
```

On real data:

```bash
python scripts/00_fetch_public_data.py --accession PRJNAxxxxxx --metadata-only
# curate data/metadata/*_samples.tsv by hand, then
python scripts/00_fetch_public_data.py --accession PRJNAxxxxxx --max-samples 40
Rscript R/01_dada2.R data/raw data/processed 240 200 4
Rscript R/02_assign_taxonomy.R data/processed path/to/silva_trainset.fa.gz
Rscript R/03_differential_abundance.R data/processed results/tables arm donor
```

`make all` runs the simulated route end to end.

---

## Results from the demonstration run

Simulated study: 30 donors sampled at baseline and after treatment, 300
taxa, 30 of them truly changed by a log fold change of 1.5 in one
direction, plus 12 reagent contaminants, 6 negative controls, 3 batches,
and a genuine 40 percent drop in absolute bacterial load under treatment.
Library sizes span 326 to 97,530 reads. 22.1 percent of the count matrix
is zero.

**Result 1. Contaminants are recoverable when controls exist.** The
combined prevalence and frequency method recovered all 12 planted
contaminants with no false flags. They carried only 0.59 percent of reads
here, but their relative abundance rises automatically when real biomass
falls, which is exactly the mechanism by which a contaminant impersonates
a treatment effect in a drug arm.

**Result 2. Donor dominates, so the design dictates the model.** Median
share of CLR variance per feature:

| factor | median variance explained |
|---|---|
| donor | 0.843 |
| batch | 0.098 |
| treatment arm | 0.007 |

PERMANOVA on Aitchison distance, permuted within donor, gives F = 2.03,
p = 0.001 for the arm. The treatment signal is real and it is small
relative to between-person variation. Median intraclass correlation for
donor across taxa is 0.741.

**Result 3. The methods agree on direction and disagree on threshold.**

| | ALDEx2 | ANCOM-BC2 |
|---|---|---|
| features called at q < 0.05 | 23 | 28 |
| true positives | 23 | 28 |
| false positives | 0 | 0 |
| sensitivity | 0.767 | 0.933 |
| realised FDP | 0.000 | 0.000 |

Jaccard overlap 0.759, Spearman correlation between the ALDEx2 effect and
the ANCOM-BC2 log fold change 0.902. Of the 7 disagreements, 6 were
ANCOM-BC2 calling something ALDEx2 declined to call and 1 the reverse.
Their median prevalence was 0.87 against 0.93 for all features, so the
disagreements sit at the sparser end, as expected. The estimated log
sampling-fraction difference was −0.300.

**Result 4. What a random effect does depends on the design, and it does
not fix compositional bias either way.**

Two things get conflated in this area, so both are measured separately.

*The paired design, which is the demonstration dataset.* Each donor gives
a baseline and a treated sample, so donor is crossed with the arm.
Blocking on donor removes between-donor variance from the residual and
therefore **increases** power. Fitting `clr ~ arm + (1 | donor)` called 47
features against 22 for ordinary least squares, with 26 calls unique to
the mixed model and only 1 unique to OLS. Nothing here is inflated by
ignoring donor; findings are simply lost.

*The nested design, which is genuine pseudoreplication.* Each donor
belongs to one arm and contributes three samples, so the replicates are
correlated within arm. Here ignoring donor is a real error:

| | called | TP | FP | realised FDP |
|---|---|---|---|---|
| donor modelled | 2 | 2 | 0 | 0.00 |
| donor ignored | 24 | 12 | 12 | 0.50 |

Half of everything reported is false. `scripts/05_random_effects.py`
runs both, via `simulate_study(design=...)`.

*And separately, the random effect does not correct compositional bias.*
In the paired run the mixed model called 29 true features plus 18 false
ones, a realised FDP of 0.383. All 18 false positives had negative
coefficients, and the mean coefficient across all genuinely null taxa was
−0.164. That is the CLR reference shift: 30 taxa rising together dragged
the sample geometric mean up, so unchanged taxa appear to fall.
ANCOM-BC2's bias correction pulls the mean null effect from −0.337 to
−0.037 and the false positives to zero on the same data.

Pseudoreplication and compositional bias are different problems. A random
effect addresses the first and does nothing about the second.

**Result 5. Calibration and power over repeated simulations.** Four
replicates per condition, 30 donors, 200 taxa, 20 truly changed:

| true log fold change | ALDEx2 sensitivity | ANCOM-BC2 sensitivity | ALDEx2 FDP | ANCOM-BC2 FDP |
|---|---|---|---|---|
| 0.5 | 0.000 | 0.000 | 0.000 | see note |
| 1.0 | 0.013 | 0.312 | 0.000 | 0.000 |
| 2.0 | 0.850 | 1.000 | 0.000 | 0.012 |

Note on the 0.5 row: ANCOM-BC2 averaged 0.5 features called per
replicate there, and one replicate called a single false feature, giving
a mean FDP of 0.50 in `simulation_summary.tsv`. A false discovery
proportion computed on one discovery is not a calibration estimate and
should not be read as one. With this few calls the quantity is
undefined in practice; raise the replicate count before quoting it.

Elsewhere both methods control the false discovery rate at or below
nominal. The power gap at moderate effect sizes is large. If you use ALDEx2 alone on a
study with effects near a log fold change of 1, you will find almost
nothing, and the absence will not be evidence of absence.

Figures are written to `results/figures/`: `decontam.png`,
`ordination.png`, `volcano.png`, `effect_concordance.png`,
`random_effects.png`, `simulation_calibration.png`.

---

## How to read a disagreement

`results/tables/disagreements.tsv` lists only the features where the two
methods differ, sorted by prevalence. The practical reading:

- Called by both, passes pseudocount sensitivity, absolute ALDEx2 effect
  above 1: report it.
- Called by ANCOM-BC2 only, at low prevalence: treat as a lead, not a
  finding. This is where ALDEx2's wide Dirichlet posterior is telling you
  the counts are too sparse to support the claim.
- Called by ALDEx2 only: rare, and usually means the ANCOM-BC2 bias
  estimate absorbed the effect. Check whether the community moved
  one-sidedly.
- Fails pseudocount sensitivity: do not report it regardless of q value.
  A finding that depends on whether you added 0.5 or 2 to the counts is
  not a finding.

Agreement between the two is not evidence of correctness. Both rest on
log ratios and both share the same blind spot, which is the next section.

---

## Assumptions and limits

The full list is in [`docs/ASSUMPTIONS.md`](docs/ASSUMPTIONS.md). The
ones that should never be omitted from a write-up:

**Absolute abundance is unrecoverable.** Every method here reports
relative change. In the demonstration run the treated arm genuinely had
40 percent less bacterial load, and nothing in the output reflects that,
because no compositional method can see it. A taxon reported as
"increased" may have decreased in absolute terms while decreasing less
than its neighbours. Without flow cytometry, qPCR of total 16S copies, or
a spike-in, this ambiguity cannot be resolved and must be stated in the
results, not buried in the methods.

**16S does not resolve species.** Genus-level calls from a 250 bp region
are usually reliable, species-level calls usually are not. The pipeline
uses exact-match species assignment only, and even a unique exact match
does not identify a unique organism.

**Databases disagree with each other.** SILVA is curated on rRNA
phylogeny; GTDB is built on whole-genome relatedness and renames and
splits many genera. A genus name in one is not the same entity in the
other. Never mix the two in one analysis, and name the exact release in
every caption. `R/02_assign_taxonomy.R` writes the reference path, its
md5 and the package version to `taxonomy_provenance.txt` for this reason.

**Reference databases are biased towards cultured organisms.**
Unassigned does not mean absent, it means unrepresented, and the
unassigned fraction is not random with respect to body site or host.

**16S copy number varies between organisms**, from one to more than
fifteen. Read counts are therefore not proportional to cell counts even
within a sample. Copy-number correction exists and is itself
model-dependent, so this pipeline does not apply it and instead reports
sequence abundance, not cell abundance.

**Filtering changes results.** The prevalence filter changes which
features enter the CLR geometric mean and therefore changes every CLR
value. The threshold is a parameter, it is reported, and conclusions
should be checked across a range of it.

**The simulation is a model, not reality.** It captures overdispersion,
sparsity, donor effects, batch, contamination and compositionality. It
does not capture phylogenetic correlation between taxa, non-linear
dose-response, or the peculiarities of any specific body site. Good
performance here is necessary, not sufficient.

---

## Installation

```
python >= 3.10
numpy, scipy, pandas, statsmodels, matplotlib, pytest
```

For the R route: R >= 4.3 with Bioconductor packages `dada2`, `phyloseq`,
`decontam`, `ALDEx2`, `ANCOMBC`, `TreeSummarizedExperiment`. See
`environment.yml`.

## Licence

MIT. See `LICENSE`.

## References

Aitchison J (1986). *The Statistical Analysis of Compositional Data.*
Gloor GB et al. (2017). Microbiome datasets are compositional. *Front Microbiol.*
Fernandes AD et al. (2014). Unifying the analysis of high-throughput sequencing datasets. *Microbiome.*
Lin H, Peddada SD (2020). Analysis of compositions of microbiomes with bias correction. *Nat Commun.*
Lin H, Peddada SD (2024). Multigroup analysis of compositions of microbiomes with covariate adjustments and repeated measures. *Nat Methods.*
Davis NM et al. (2018). Simple statistical identification and removal of contaminant sequences. *Microbiome.*
Callahan BJ et al. (2016). DADA2: high-resolution sample inference. *Nat Methods.*
Nearing JT et al. (2022). Microbiome differential abundance methods produce different results across 38 datasets. *Nat Commun.*
