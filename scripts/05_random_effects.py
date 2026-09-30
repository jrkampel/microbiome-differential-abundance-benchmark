#!/usr/bin/env python3
"""Quantify what modelling donor as a random effect changes.

The comparison is deliberately blunt. The same CLR values are tested
twice, once with a random intercept for donor and once with ordinary
least squares that pretends the samples are independent. If the OLS fit
calls many more features, that difference is the number of findings the
study would have reported purely by treating repeated measures as
replicates.

The intraclass correlation is reported per taxon. It is the fraction of
CLR variance attributable to donor identity, and the argument for the
mixed model is quantitative, not stylistic: with an ICC of 0.6 and two
samples per donor the effective sample size is close to the number of
donors, not the number of samples.
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

from mdab import mixed_model_scan, prevalence_filter  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", type=Path, default=Path("data/processed"))
    ap.add_argument("--counts", default="counts_decontaminated.tsv")
    ap.add_argument("--group-column", default="arm")
    ap.add_argument("--random-effect", default="donor")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--min-prevalence", type=float, default=0.10)
    ap.add_argument("--nested-demo", action="store_true", default=True,
                    help="also simulate a nested design, where ignoring donor "
                         "is genuine pseudoreplication")
    ap.add_argument("--results", type=Path, default=Path("results"))
    args = ap.parse_args()

    counts = pd.read_csv(args.indir / args.counts, sep="\t", index_col=0)
    meta = pd.read_csv(args.indir / "metadata.tsv", sep="\t", index_col=0).loc[counts.index]
    counts = prevalence_filter(counts, args.min_prevalence)

    print(f"fitting {counts.shape[1]} mixed models ...")
    res = mixed_model_scan(counts, meta[args.group_column], meta[args.random_effect])

    tdir = args.results / "tables"
    fdir = args.results / "figures"
    tdir.mkdir(parents=True, exist_ok=True)
    fdir.mkdir(parents=True, exist_ok=True)
    res.to_csv(tdir / "mixed_model_results.tsv", sep="\t")

    n_mixed = int((res["q_mixed"] < args.alpha).sum())
    n_ols = int((res["q_ols"] < args.alpha).sum())
    only_ols = int(((res["q_ols"] < args.alpha) & (res["q_mixed"] >= args.alpha)).sum())
    only_mixed = int(((res["q_mixed"] < args.alpha) & (res["q_ols"] >= args.alpha)).sum())
    icc = res["icc_donor"].dropna()

    print(f"converged: {int(res['converged'].sum())}/{len(res)}")
    print(f"median ICC for donor: {icc.median():.3f} "
          f"(quartiles {icc.quantile(0.25):.3f}, {icc.quantile(0.75):.3f})")
    print(f"features called by mixed model: {n_mixed}")
    print(f"features called ignoring donor: {n_ols}")
    print(f"called only when donor is ignored: {only_ols}")
    print(f"called only when donor is modelled: {only_mixed}")

    truth_path = args.indir / "truth.tsv"
    if truth_path.exists():
        truth = pd.read_csv(truth_path, sep="\t", index_col=0)["is_da"].astype(bool)
        t = truth.reindex(res.index).fillna(False)
        for label, col in (("mixed", "q_mixed"), ("ols", "q_ols")):
            call = res[col] < args.alpha
            tp = int((call & t).sum())
            fp = int((call & ~t).sum())
            fdp = fp / max(tp + fp, 1)
            print(f"  {label}: TP {tp}, FP {fp}, realised FDP {fdp:.3f}")

    pd.Series(
        {
            "n_features": len(res),
            "median_icc_donor": float(icc.median()),
            "n_called_mixed": n_mixed,
            "n_called_ols": n_ols,
            "n_only_ols": only_ols,
            "n_only_mixed": only_mixed,
        }
    ).to_frame("value").to_csv(tdir / "random_effect_summary.tsv", sep="\t")

    # The nested contrast. In the paired design above, blocking on donor
    # removes between-donor variance from the residual and gains power.
    # That is not pseudoreplication. Pseudoreplication is the nested case,
    # where each donor sits in one arm and contributes several samples, so
    # ignoring donor treats correlated replicates as independent. The two
    # are routinely conflated, so both are shown.
    if args.nested_demo:
        from mdab import simulate_study
        print("\nnested design contrast (donors nested within arm):")
        nested = simulate_study(design="nested", n_reps=3, n_donors=20,
                                n_taxa=150, n_da=15, n_controls=0,
                                n_contaminants=0, seed=99)
        nbio = nested.biological_samples
        ncounts = prevalence_filter(nested.counts.loc[nbio], 0.1)
        nres = mixed_model_scan(ncounts, nested.metadata.loc[nbio, "arm"],
                                nested.metadata.loc[nbio, "donor"])
        ntruth = nested.truth["is_da"].reindex(nres.index).fillna(False)
        for label, col in (("donor modelled", "q_mixed"), ("donor ignored", "q_ols")):
            call = nres[col] < args.alpha
            tp = int((call & ntruth).sum())
            fp = int((call & ~ntruth).sum())
            print(f"  {label:15s} called {int(call.sum()):3d}  TP {tp:3d}  FP {fp:3d}  "
                  f"realised FDP {fp / max(tp + fp, 1):.3f}")
        nres.to_csv(tdir / "mixed_model_nested_results.tsv", sep="\t")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    ax = axes[0]
    ax.hist(icc, bins=30, color="#4C72B0")
    ax.axvline(icc.median(), color="crimson", ls="--",
               label=f"median {icc.median():.2f}")
    ax.set_xlabel("ICC for donor")
    ax.set_ylabel("features")
    ax.set_title("Variance attributable to donor")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    lo = 1e-12
    ax.scatter(-np.log10(res["p_ols"].clip(lo)), -np.log10(res["p_mixed"].clip(lo)),
               s=14, alpha=0.7, color="#55A868")
    lim = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([0, lim], [0, lim], color="grey", lw=0.8)
    ax.axhline(-np.log10(args.alpha), ls=":", lw=0.8, color="grey")
    ax.axvline(-np.log10(args.alpha), ls=":", lw=0.8, color="grey")
    ax.set_xlabel("-log10 p, donor ignored")
    ax.set_ylabel("-log10 p, donor modelled")
    ax.set_title("Paired design: blocking on donor gains power")

    ax = axes[2]
    ax.scatter(res["se_ols"], res["se_mixed"], s=14, alpha=0.7, color="#C44E52")
    lim = max(res["se_ols"].max(), res["se_mixed"].max())
    ax.plot([0, lim], [0, lim], color="grey", lw=0.8)
    ax.set_xlabel("standard error, donor ignored")
    ax.set_ylabel("standard error, donor modelled")
    ax.set_title("Standard errors")
    fig.tight_layout()
    fig.savefig(fdir / "random_effects.png", dpi=160)
    print(f"figure written to {fdir / 'random_effects.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
