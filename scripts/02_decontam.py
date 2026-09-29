#!/usr/bin/env python3
"""Identify and remove reagent contaminants using negative controls.

Run this before any differential abundance testing, never after. A
contaminant that survives into the test will often be significant,
because its relative abundance tracks biomass and biomass tracks the
treatment.

If the study has no negative controls, this script exits with a clear
message rather than inventing a substitute. There is no sound way to
identify contaminants without either controls or DNA quantitation, and
claiming otherwise in a README is worse than admitting the gap. Record
the absence in ``docs/ASSUMPTIONS.md``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mdab import is_contaminant  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", type=Path, default=Path("data/processed"))
    ap.add_argument("--threshold", type=float, default=0.1,
                    help="decontam score cutoff; the package default of 0.1 is "
                         "permissive on purpose")
    ap.add_argument("--method", default="combined",
                    choices=["prevalence", "frequency", "combined", "either"])
    ap.add_argument("--results", type=Path, default=Path("results"))
    args = ap.parse_args()

    counts = pd.read_csv(args.indir / "counts.tsv", sep="\t", index_col=0)
    meta = pd.read_csv(args.indir / "metadata.tsv", sep="\t", index_col=0)

    if not meta["is_control"].any():
        print("no negative controls in this dataset; skipping decontam.")
        print("record this limitation in docs/ASSUMPTIONS.md and do not")
        print("substitute a prevalence filter for contaminant removal.")
        counts.to_csv(args.indir / "counts_decontaminated.tsv", sep="\t")
        return 0

    res = is_contaminant(
        counts,
        is_control=meta["is_control"].to_numpy(bool),
        concentration=meta["dna_ng_ul"].to_numpy(float) if "dna_ng_ul" in meta else None,
        method=args.method,
        threshold=args.threshold,
    )

    (args.results / "tables").mkdir(parents=True, exist_ok=True)
    (args.results / "figures").mkdir(parents=True, exist_ok=True)
    res.table.sort_values("score").to_csv(args.results / "tables" / "decontam_scores.tsv",
                                          sep="\t")

    clean = res.clean(counts)
    clean = clean.loc[~meta["is_control"].to_numpy(bool)]
    clean = clean.loc[:, clean.sum(axis=0) > 0]
    clean.to_csv(args.indir / "counts_decontaminated.tsv", sep="\t")

    n_removed = len(res.contaminants)
    reads_removed = counts[res.contaminants].to_numpy().sum() / counts.to_numpy().sum()
    print(f"method {res.method}, threshold {args.threshold}")
    print(f"flagged {n_removed} of {counts.shape[1]} features as contaminants")
    print(f"they carry {reads_removed:.2%} of all reads")
    print(f"remaining after removing controls: {clean.shape[0]} samples, "
          f"{clean.shape[1]} features")

    if (args.indir / "truth.tsv").exists():
        truth = pd.read_csv(args.indir / "truth.tsv", sep="\t", index_col=0)
        t = truth["is_contaminant"].reindex(res.table.index).fillna(False)
        called = res.table["contaminant"]
        tp = int((called & t).sum())
        fp = int((called & ~t).sum())
        fn = int((~called & t).sum())
        print(f"against simulated truth: recovered {tp}/{tp + fn} contaminants, "
              f"{fp} false flags")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].hist(res.table["score"].dropna(), bins=40, color="#4C72B0")
    ax[0].axvline(args.threshold, color="crimson", ls="--",
                  label=f"threshold {args.threshold}")
    ax[0].set_xlabel("decontam score (low means contaminant)")
    ax[0].set_ylabel("features")
    ax[0].set_title("Score distribution")
    ax[0].legend(frameon=False)

    prev_ctrl = (counts.loc[meta["is_control"].to_numpy(bool)] > 0).mean(axis=0)
    prev_samp = (counts.loc[~meta["is_control"].to_numpy(bool)] > 0).mean(axis=0)
    flagged = res.table["contaminant"].reindex(prev_ctrl.index).fillna(False)
    ax[1].scatter(prev_samp[~flagged], prev_ctrl[~flagged], s=12, alpha=0.5,
                  color="#55A868", label="retained")
    ax[1].scatter(prev_samp[flagged], prev_ctrl[flagged], s=18, alpha=0.9,
                  color="crimson", label="flagged")
    ax[1].plot([0, 1], [0, 1], color="grey", lw=0.8)
    ax[1].set_xlabel("prevalence in true samples")
    ax[1].set_ylabel("prevalence in negative controls")
    ax[1].set_title("Prevalence contrast")
    ax[1].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(args.results / "figures" / "decontam.png", dpi=160)
    print(f"figure written to {args.results / 'figures' / 'decontam.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
