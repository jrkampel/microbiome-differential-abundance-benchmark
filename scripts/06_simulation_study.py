#!/usr/bin/env python3
"""Repeat the simulation many times and measure calibration and power.

A single dataset shows where two methods disagree. It cannot show which
one is better calibrated, because the realised false discovery proportion
in one dataset is a noisy draw. This script repeats the whole pipeline
across replicate simulations and across a grid of conditions, then
reports the mean false discovery proportion against the nominal level and
the sensitivity at that level.

Conditions varied
-----------------
``effect``           size of the true log fold change
``n_donors``         sample size
``asymmetric``       whether all true effects share a direction, which
                     stresses the ANCOM-BC2 bias assumption
``biomass_ratio``    whether total load really changes, which no
                     compositional method can detect and which therefore
                     tests nothing about the methods but everything about
                     how results should be worded

This is the figure that makes the repository more than a tutorial. It
turns "the methods disagree" into a statement about which disagreement
costs you.
"""

from __future__ import annotations

import argparse
import itertools
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mdab import (  # noqa: E402
    aldex2,
    ancombc2,
    concordance_table,
    evaluate_against_truth,
    prevalence_filter,
    simulate_study,
)


def _sem(values) -> float:
    """Standard error of the mean, zero when there is a single replicate."""
    n = len(values)
    return float(np.std(values, ddof=1) / np.sqrt(n)) if n > 1 else 0.0


def one_replicate(seed: int, effect: float, n_donors: int, asymmetric: bool,
                  biomass_ratio: float, alpha: float, mc_samples: int,
                  use_random_effect: bool) -> pd.DataFrame:
    study = simulate_study(
        n_donors=n_donors,
        n_taxa=200,
        n_da=20,
        effect_lfc=effect,
        asymmetric=asymmetric,
        biomass_ratio=biomass_ratio,
        n_controls=0,
        n_contaminants=0,
        seed=seed,
    )
    bio = study.biological_samples
    counts = prevalence_filter(study.counts.loc[bio], 0.10)
    meta = study.metadata.loc[bio]

    a = aldex2(counts, meta["arm"], mc_samples=mc_samples, seed=seed)
    b = ancombc2(
        counts,
        meta["arm"],
        random_effect=meta["donor"] if use_random_effect else None,
        sensitivity_pseudocounts=(0.5, 1.0, 2.0),
        prv_cut=0.0,
    )
    conc = concordance_table(a.table, b.table, alpha=alpha)
    perf = evaluate_against_truth(conc, study.truth["is_da"].astype(bool))
    perf = perf.reset_index()
    perf["seed"] = seed
    perf["effect"] = effect
    perf["n_donors"] = n_donors
    perf["asymmetric"] = asymmetric
    perf["biomass_ratio"] = biomass_ratio
    perf["delta_hat"] = b.delta
    return perf


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--replicates", type=int, default=10)
    ap.add_argument("--effects", type=float, nargs="+", default=[0.5, 1.0, 2.0])
    ap.add_argument("--donors", type=int, nargs="+", default=[10, 30])
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--mc-samples", type=int, default=64)
    ap.add_argument("--asymmetric", type=int, nargs="+", default=[1],
                    help="1 for one-sided true effects, 0 for balanced, both to compare")
    ap.add_argument("--no-random-effect", action="store_true")
    ap.add_argument("--results", type=Path, default=Path("results"))
    args = ap.parse_args()

    grid = list(itertools.product(args.effects, args.donors, args.asymmetric,
                                  range(args.replicates)))
    print(f"{len(grid)} fits, this takes a few minutes")

    rows = []
    t0 = time.time()
    for k, (effect, n_donors, asym, rep) in enumerate(grid, 1):
        rows.append(
            one_replicate(
                seed=1000 * rep + k,
                effect=effect,
                n_donors=n_donors,
                asymmetric=bool(asym),
                biomass_ratio=1.0,
                alpha=args.alpha,
                mc_samples=args.mc_samples,
                use_random_effect=not args.no_random_effect,
            )
        )
        if k % 5 == 0 or k == len(grid):
            print(f"  {k}/{len(grid)}  {time.time() - t0:.0f}s")

    res = pd.concat(rows, ignore_index=True)
    tdir = args.results / "tables"
    tdir.mkdir(parents=True, exist_ok=True)
    res.to_csv(tdir / "simulation_raw.tsv", sep="\t", index=False)

    summary = (
        res.groupby(["method", "effect", "n_donors", "asymmetric"])
        .agg(
            mean_FDP=("FDP", "mean"),
            se_FDP=("FDP", _sem),
            mean_sensitivity=("sensitivity", "mean"),
            mean_called=("n_called", "mean"),
            mean_delta_hat=("delta_hat", "mean"),
        )
        .reset_index()
    )
    summary.to_csv(tdir / "simulation_summary.tsv", sep="\t", index=False)
    print("\n" + summary.round(3).to_string(index=False))

    fdir = args.results / "figures"
    fdir.mkdir(parents=True, exist_ok=True)
    methods = ["ALDEx2", "ANCOM-BC2", "intersection", "union"]
    colours = dict(zip(methods, ["#C44E52", "#4C72B0", "#8172B2", "#55A868"]))

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    for method in methods:
        sub = summary[summary["method"] == method].groupby("effect").mean(numeric_only=True)
        axes[0].errorbar(sub.index, sub["mean_FDP"], yerr=sub["se_FDP"],
                         marker="o", capsize=3, color=colours[method], label=method)
        axes[1].plot(sub.index, sub["mean_sensitivity"], marker="o",
                     color=colours[method], label=method)
    axes[0].axhline(args.alpha, ls="--", color="grey", label=f"nominal {args.alpha}")
    axes[0].set_xlabel("true log fold change")
    axes[0].set_ylabel("mean false discovery proportion")
    axes[0].set_title("Calibration")
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].set_xlabel("true log fold change")
    axes[1].set_ylabel("mean sensitivity")
    axes[1].set_title("Power")
    axes[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(fdir / "simulation_calibration.png", dpi=160)
    print(f"\nfigure written to {fdir / 'simulation_calibration.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
