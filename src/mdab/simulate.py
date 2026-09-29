"""Generate 16S-like count data with a known differential abundance truth.

Purpose
-------
Every comparison of differential abundance methods on real data suffers
from the same problem: there is no ground truth, so agreement between
methods is reported as if it were accuracy. Simulation fixes that. The
generator below produces data with the properties that actually break
these methods, so false discovery rate and power can be measured rather
than assumed.

What is modelled
----------------
1. **Absolute abundance, then sampling.** Each specimen has a latent
   absolute community. Reads are drawn multinomially from it with a
   library size that varies by an order of magnitude. This is the
   mechanism that makes the data compositional. Nothing downstream ever
   sees the latent values, but the evaluation scripts do.

2. **Donor random effects.** Taxon log-abundances carry a per-donor
   offset with standard deviation ``donor_sd``. Inter-individual variation
   in gut microbiomes is large relative to most treatment effects, which
   is why paired or repeated-measures designs dominate the field and why
   ignoring donor is not a minor sin.

3. **Overdispersion.** Counts are far more variable than a multinomial
   allows. Modelled as lognormal noise on the latent abundance before
   sampling, giving a mean-variance relationship close to what is seen in
   real 16S data.

4. **Sparsity.** A large fraction of taxa are rare, so zeros arise both
   because a taxon is absent and because it was present but not sampled.
   Distinguishing the two is impossible from the data alone, which is the
   whole reason zero handling is contested.

5. **A batch effect** applied multiplicatively to a subset of taxa, to
   test whether modelling batch as a random effect recovers the truth.

6. **Reagent contaminants** at fixed absolute load, so their relative
   abundance rises automatically when real biomass falls. Negative
   control samples contain contaminants only.

7. **Optional biomass change under treatment.** When ``biomass_ratio`` is
   not 1, the treatment arm has genuinely less total bacterial load. No
   compositional method can detect this, and it is the single most
   important caveat to put in a README.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

__all__ = ["simulate_study", "SimulatedStudy"]


@dataclass
class SimulatedStudy:
    counts: pd.DataFrame
    metadata: pd.DataFrame
    truth: pd.DataFrame
    latent: pd.DataFrame

    @property
    def biological_samples(self) -> pd.Index:
        return self.metadata.index[self.metadata["sample_type"] == "sample"]


def simulate_study(
    n_donors: int = 30,
    n_taxa: int = 300,
    n_da: int = 30,
    effect_lfc: float = 1.5,
    donor_sd: float = 1.2,
    overdispersion_sd: float = 0.6,
    n_batches: int = 3,
    batch_sd: float = 0.5,
    n_controls: int = 6,
    n_contaminants: int = 12,
    contaminant_load: float = 400.0,
    biomass_ratio: float = 1.0,
    depth_log10_range: tuple[float, float] = (4.0, 5.0),
    asymmetric: bool = True,
    design: str = "paired",
    n_reps: int = 2,
    seed: int = 7,
) -> SimulatedStudy:
    """Simulate a two-arm study with repeated donors.

    Parameters
    ----------
    design : "paired" or "nested". This choice decides what a random
        effect buys you, and the two cases are routinely conflated.

        ``paired``  each donor contributes one baseline and one treated
        sample, so donor is crossed with the arm. Blocking on donor
        removes between-donor variance from the residual and therefore
        *increases* power. Ignoring donor here loses findings.

        ``nested``  each donor belongs to one arm and contributes
        ``n_reps`` samples. Donor is nested within arm. Ignoring donor
        here treats correlated replicates as independent, shrinks the
        standard error and *inflates* significance. This is
        pseudoreplication in the strict sense.

        Both designs occur in real microbiome studies. A crossover drug
        trial is paired; a cohort with several stool samples per
        participant per arm is nested.
    n_reps : samples per donor in the nested design.
    n_donors : number of donors.
    n_da : number of truly differentially abundant taxa.
    effect_lfc : log fold change applied to those taxa, in natural log
        units on absolute abundance.
    asymmetric : if True all true effects go in the same direction. This
        is the adversarial case for ANCOM-BC2, because the assumption
        that the bulk of taxa are null is still satisfied but the bias
        estimate has to separate a real one-sided signal from a sampling
        fraction shift.
    biomass_ratio : total absolute bacterial load in the treated arm
        relative to baseline. Values below 1 mimic an antibiotic.
    """
    rng = np.random.default_rng(seed)

    taxa = [f"ASV{i:04d}" for i in range(n_taxa)]
    contam_taxa = [f"CONTAM{i:03d}" for i in range(n_contaminants)]
    all_features = taxa + contam_taxa

    # Baseline community: log-normal abundances give the long tail of rare
    # taxa that characterises amplicon data.
    base_log = rng.normal(loc=4.0, scale=2.0, size=n_taxa)
    base_abundance = np.exp(base_log)

    da_idx = rng.choice(n_taxa, size=n_da, replace=False)
    if asymmetric:
        signs = np.ones(n_da)
    else:
        signs = rng.choice([-1.0, 1.0], size=n_da)
    true_lfc = np.zeros(n_taxa)
    true_lfc[da_idx] = signs * effect_lfc

    donor_offsets = rng.normal(0.0, donor_sd, size=(n_donors, n_taxa))

    batch_of_donor = rng.integers(0, n_batches, size=n_donors)
    batch_effect = rng.normal(0.0, batch_sd, size=(n_batches, n_taxa))

    contaminant_profile = rng.dirichlet(np.ones(n_contaminants) * 0.7)

    rows, meta_rows, latent_rows = [], [], []

    if design not in {"paired", "nested"}:
        raise ValueError("design must be 'paired' or 'nested'")
    donor_arm = rng.permutation(
        np.array(["baseline"] * (n_donors // 2) + ["treated"] * (n_donors - n_donors // 2))
    )

    for d in range(n_donors):
        if design == "paired":
            plan = [("baseline", 0), ("treated", 0)]
        else:
            plan = [(donor_arm[d], r) for r in range(n_reps)]

        for arm, rep in plan:
            log_abund = (
                np.log(base_abundance)
                + donor_offsets[d]
                + batch_effect[batch_of_donor[d]]
                + rng.normal(0.0, overdispersion_sd, size=n_taxa)
            )
            if arm == "treated":
                log_abund = log_abund + true_lfc
            abund = np.exp(log_abund)
            if arm == "treated":
                abund = abund * biomass_ratio

            # Reagent contamination is a fixed absolute amount, so its
            # share rises when real biomass falls. This is the mechanism
            # by which contaminants masquerade as treatment effects.
            contam_abs = contaminant_load * contaminant_profile
            full_abs = np.concatenate([abund, contam_abs])
            probs = full_abs / full_abs.sum()

            depth = int(10 ** rng.uniform(*depth_log10_range))
            counts = rng.multinomial(depth, probs)

            sample_id = (f"D{d:03d}_{arm[:3]}" if design == "paired"
                         else f"D{d:03d}_{arm[:3]}_r{rep}")
            rows.append(pd.Series(counts, index=all_features, name=sample_id))
            latent_rows.append(pd.Series(full_abs, index=all_features, name=sample_id))
            meta_rows.append(
                {
                    "sample_id": sample_id,
                    "donor": f"D{d:03d}",
                    "arm": arm,
                    "batch": f"B{batch_of_donor[d]}",
                    "sample_type": "sample",
                    "is_control": False,
                    "dna_ng_ul": float(np.exp(rng.normal(2.0, 0.4))
                                       * (biomass_ratio if arm == "treated" else 1.0)),
                    "library_size": int(counts.sum()),
                }
            )

    # Negative controls: reagent DNA only, shallow and noisy.
    for c in range(n_controls):
        contam_abs = contaminant_load * rng.dirichlet(contaminant_profile * 40 + 0.5)
        probs = np.concatenate([np.full(n_taxa, 1e-9), contam_abs])
        probs = probs / probs.sum()
        depth = int(10 ** rng.uniform(2.5, 3.5))
        counts = rng.multinomial(depth, probs)
        sample_id = f"NEG{c:02d}"
        rows.append(pd.Series(counts, index=all_features, name=sample_id))
        latent_rows.append(pd.Series(np.concatenate([np.zeros(n_taxa), contam_abs]),
                                     index=all_features, name=sample_id))
        meta_rows.append(
            {
                "sample_id": sample_id,
                "donor": "NA",
                "arm": "control",
                "batch": f"B{rng.integers(0, n_batches)}",
                "sample_type": "negative_control",
                "is_control": True,
                "dna_ng_ul": float(np.exp(rng.normal(-1.5, 0.4))),
                "library_size": int(counts.sum()),
            }
        )

    counts_df = pd.DataFrame(rows).astype("int64")
    meta_df = pd.DataFrame(meta_rows).set_index("sample_id")
    meta_df = meta_df.loc[counts_df.index]

    truth = pd.DataFrame(
        {
            "true_lfc": np.concatenate([true_lfc, np.zeros(n_contaminants)]),
            "is_da": np.concatenate([true_lfc != 0, np.zeros(n_contaminants, bool)]),
            "is_contaminant": np.concatenate(
                [np.zeros(n_taxa, bool), np.ones(n_contaminants, bool)]
            ),
            "baseline_abundance": np.concatenate(
                [base_abundance, contaminant_load * contaminant_profile]
            ),
        },
        index=all_features,
    )

    latent_df = pd.DataFrame(latent_rows)

    return SimulatedStudy(counts=counts_df, metadata=meta_df, truth=truth, latent=latent_df)
