#!/usr/bin/env python3
"""Run ALDEx2 and ANCOM-BC2 on the same data and characterise disagreement.

Outputs
-------
tables/aldex2_results.tsv
tables/ancombc2_results.tsv
tables/concordance.tsv           one row per feature with both calls
tables/agreement_summary.tsv
tables/disagreements.tsv         only the features where the calls differ
tables/method_performance.tsv    only when ground truth is available
figures/volcano.png
figures/effect_concordance.png

How to read the disagreement table
----------------------------------
Sort it by prevalence. In almost every dataset the disagreements
concentrate at low prevalence, because that is where the Dirichlet
posterior is wide and ALDEx2 declines to call, while the linear model
sees a mean shift and does. Whether that is ALDEx2 being appropriately
cautious or ANCOM-BC2 being appropriately powerful depends entirely on
whether the effects are real, which is why this repository includes a
simulator.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mdab import (  # noqa: E402
    agreement_summary,
    aldex2,
    ancombc2,
    concordance_table,
    evaluate_against_truth,
    prevalence_filter,
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", type=Path, default=Path("data/processed"))
    ap.add_argument("--counts", default="counts_decontaminated.tsv")
    ap.add_argument("--group-column", default="arm")
    ap.add_argument("--random-effect", default="donor",
                    help="metadata column for the random intercept; 'none' disables it")
    ap.add_argument("--mc-samples", type=int, default=128)
    ap.add_argument("--denom", default="all", choices=["all", "iqlr"])
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--min-effect", type=float, default=0.0,
                    help="optional ALDEx2 effect-size screen, 1.0 is the usual value")
    ap.add_argument("--min-prevalence", type=float, default=0.10)
    ap.add_argument("--results", type=Path, default=Path("results"))
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    counts = pd.read_csv(args.indir / args.counts, sep="\t", index_col=0)
    meta = pd.read_csv(args.indir / "metadata.tsv", sep="\t", index_col=0).loc[counts.index]
    counts = prevalence_filter(counts, args.min_prevalence)

    groups = meta[args.group_column]
    random = None if args.random_effect == "none" else meta[args.random_effect]

    print(f"{counts.shape[0]} samples, {counts.shape[1]} features after filtering")
    print(f"groups: {groups.value_counts().to_dict()}")

    print("running ALDEx2 ...")
    a = aldex2(counts, groups, mc_samples=args.mc_samples, denom=args.denom, seed=args.seed)

    print("running ANCOM-BC2 ...")
    b = ancombc2(counts, groups, random_effect=random, prv_cut=0.0)
    print(f"estimated log sampling-fraction difference delta = {b.delta:.4f} "
          f"(se {b.delta_se:.4f})")
    print(f"structural zeros: {int(b.table['structural_zero'].sum())}")
    print(f"failed pseudocount sensitivity: {int((~b.table['passed_ss']).sum())}")

    tdir = args.results / "tables"
    fdir = args.results / "figures"
    tdir.mkdir(parents=True, exist_ok=True)
    fdir.mkdir(parents=True, exist_ok=True)

    a.table.to_csv(tdir / "aldex2_results.tsv", sep="\t")
    b.table.to_csv(tdir / "ancombc2_results.tsv", sep="\t")
    b.sensitivity.to_csv(tdir / "ancombc2_pseudocount_sensitivity.tsv", sep="\t")

    conc = concordance_table(a.table, b.table, alpha=args.alpha, min_effect=args.min_effect)
    conc.to_csv(tdir / "concordance.tsv", sep="\t")

    summary = agreement_summary(conc)
    summary.to_frame("value").to_csv(tdir / "agreement_summary.tsv", sep="\t")
    print("\nagreement:")
    print(summary.round(3).to_string())

    disagree = conc[conc["agreement"].isin(["aldex2_only", "ancombc2_only"])]
    disagree.sort_values("prevalence").to_csv(tdir / "disagreements.tsv", sep="\t")
    if len(disagree):
        print(f"\n{len(disagree)} disagreements; median prevalence "
              f"{disagree['prevalence'].median():.2f} against "
              f"{conc['prevalence'].median():.2f} overall")

    truth_path = args.indir / "truth.tsv"
    if truth_path.exists():
        truth = pd.read_csv(truth_path, sep="\t", index_col=0)["is_da"].astype(bool)
        perf = evaluate_against_truth(conc, truth)
        perf.to_csv(tdir / "method_performance.tsv", sep="\t")
        print("\nperformance against simulated truth:")
        print(perf.round(3).to_string())

    # Volcano panels.
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    ax = axes[0]
    sig = conc["aldex_call"]
    ax.scatter(conc.loc[~sig, "aldex_effect"], -np.log10(conc.loc[~sig, "aldex_q"].clip(1e-12)),
               s=12, color="#BBBBBB", label="not called")
    ax.scatter(conc.loc[sig, "aldex_effect"], -np.log10(conc.loc[sig, "aldex_q"].clip(1e-12)),
               s=18, color="#C44E52", label="called")
    ax.axhline(-np.log10(args.alpha), ls="--", lw=0.8, color="grey")
    ax.set_xlabel("ALDEx2 standardised effect")
    ax.set_ylabel("-log10 expected BH q")
    ax.set_title("ALDEx2")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    sig = conc["ancom_call"]
    ax.scatter(conc.loc[~sig, "ancom_lfc"], -np.log10(conc.loc[~sig, "ancom_q"].clip(1e-12)),
               s=12, color="#BBBBBB", label="not called")
    ax.scatter(conc.loc[sig, "ancom_lfc"], -np.log10(conc.loc[sig, "ancom_q"].clip(1e-12)),
               s=18, color="#4C72B0", label="called")
    ax.axhline(-np.log10(args.alpha), ls="--", lw=0.8, color="grey")
    ax.set_xlabel("ANCOM-BC2 bias-corrected log fold change")
    ax.set_ylabel("-log10 BH q")
    ax.set_title("ANCOM-BC2")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(fdir / "volcano.png", dpi=160)

    # Concordance panel.
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    colours = {"both": "#8172B2", "aldex2_only": "#C44E52",
               "ancombc2_only": "#4C72B0", "neither": "#CCCCCC"}
    ax = axes[0]
    for lvl, colour in colours.items():
        sel = conc["agreement"] == lvl
        if sel.any():
            ax.scatter(conc.loc[sel, "ancom_lfc"], conc.loc[sel, "aldex_effect"],
                       s=18, color=colour, label=f"{lvl} ({int(sel.sum())})", alpha=0.85)
    ax.axhline(0, lw=0.6, color="grey")
    ax.axvline(0, lw=0.6, color="grey")
    ax.set_xlabel("ANCOM-BC2 log fold change")
    ax.set_ylabel("ALDEx2 standardised effect")
    ax.set_title(f"Effect concordance (Spearman {summary['spearman_effect_vs_lfc']:.2f})")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    for lvl, colour in colours.items():
        sel = conc["agreement"] == lvl
        if sel.any():
            ax.scatter(conc.loc[sel, "prevalence"], conc.loc[sel, "aldex_effect"].abs(),
                       s=18, color=colour, alpha=0.85, label=lvl)
    ax.set_xlabel("prevalence across samples")
    ax.set_ylabel("|ALDEx2 effect|")
    ax.set_title("Disagreement concentrates at low prevalence")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(fdir / "effect_concordance.png", dpi=160)

    print(f"\ntables in {tdir}, figures in {fdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
