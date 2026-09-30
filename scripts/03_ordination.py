#!/usr/bin/env python3
"""CLR transform, Aitchison ordination, PERMANOVA and variance partition.

This is the step that tells you whether the study is analysable at all.
If donor or batch dominates the ordination and the treatment arm is
invisible, no differential abundance method will rescue it, and the
correct thing to write in the README is that the design is confounded.

PERMANOVA here is run two ways. The unrestricted version permutes sample
labels freely and is wrong for a paired design, because it ignores that
donors are blocks. The restricted version permutes only within donor,
which is the valid null for a repeated-measures layout. Reporting both
makes the point concrete.
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

from mdab import aitchison_distance, clr, prevalence_filter, variance_partition  # noqa: E402


def permanova(dist: np.ndarray, labels: np.ndarray, blocks: np.ndarray | None,
              n_perm: int = 999, seed: int = 0) -> tuple[float, float]:
    """PERMANOVA pseudo-F with optional within-block permutation."""
    rng = np.random.default_rng(seed)

    def pseudo_f(lab):
        n = len(lab)
        total_ss = (dist ** 2).sum() / (2 * n)
        within = 0.0
        groups = np.unique(lab)
        for g in groups:
            idx = np.flatnonzero(lab == g)
            if len(idx) < 2:
                continue
            sub = dist[np.ix_(idx, idx)]
            within += (sub ** 2).sum() / (2 * len(idx))
        between = total_ss - within
        a, b = len(groups) - 1, n - len(groups)
        return (between / a) / (within / b) if within > 0 and a > 0 and b > 0 else np.nan

    observed = pseudo_f(labels)
    null = np.empty(n_perm)
    for i in range(n_perm):
        if blocks is None:
            perm = rng.permutation(labels)
        else:
            perm = labels.copy()
            for blk in np.unique(blocks):
                idx = np.flatnonzero(blocks == blk)
                perm[idx] = rng.permutation(labels[idx])
        null[i] = pseudo_f(perm)
    p = (np.sum(null >= observed) + 1) / (n_perm + 1)
    return float(observed), float(p)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", type=Path, default=Path("data/processed"))
    ap.add_argument("--counts", default="counts_decontaminated.tsv")
    ap.add_argument("--min-prevalence", type=float, default=0.10)
    ap.add_argument("--results", type=Path, default=Path("results"))
    args = ap.parse_args()

    counts = pd.read_csv(args.indir / args.counts, sep="\t", index_col=0)
    meta = pd.read_csv(args.indir / "metadata.tsv", sep="\t", index_col=0).loc[counts.index]

    before = counts.shape[1]
    counts = prevalence_filter(counts, args.min_prevalence)
    print(f"prevalence filter at {args.min_prevalence:.0%}: {before} -> {counts.shape[1]} features")

    z = clr(counts)
    dist = aitchison_distance(counts)

    centred = z - z.mean(axis=0)
    u, s, vt = np.linalg.svd(centred.to_numpy(), full_matrices=False)
    scores = u * s
    var_explained = s ** 2 / (s ** 2).sum()

    arms = meta["arm"].to_numpy()
    donors = meta["donor"].to_numpy()
    batches = meta["batch"].to_numpy()

    f_free, p_free = permanova(dist.to_numpy(), arms, None)
    f_block, p_block = permanova(dist.to_numpy(), arms, donors)
    f_batch, p_batch = permanova(dist.to_numpy(), batches, None)

    print(f"PERMANOVA arm, free permutation      F={f_free:.2f} p={p_free:.3f}  (invalid here)")
    print(f"PERMANOVA arm, permuted within donor F={f_block:.2f} p={p_block:.3f}  (valid)")
    print(f"PERMANOVA batch, free permutation    F={f_batch:.2f} p={p_batch:.3f}")

    vp = variance_partition(counts, meta[["arm", "donor", "batch"]])
    (args.results / "tables").mkdir(parents=True, exist_ok=True)
    vp.to_csv(args.results / "tables" / "variance_partition.tsv", sep="\t")
    print("\nmedian share of CLR variance explained per feature:")
    print(vp.median().round(3).to_string())

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    for ax, colour_by in zip(axes[:2], ["arm", "batch"]):
        for lvl in pd.unique(meta[colour_by]):
            sel = (meta[colour_by] == lvl).to_numpy()
            ax.scatter(scores[sel, 0], scores[sel, 1], s=28, alpha=0.85, label=str(lvl))
        ax.set_xlabel(f"PC1 ({var_explained[0]:.1%})")
        ax.set_ylabel(f"PC2 ({var_explained[1]:.1%})")
        ax.set_title(f"Aitchison PCA by {colour_by}")
        ax.legend(frameon=False, fontsize=8)

    # Paired lines make the within-donor shift visible even when the
    # between-donor spread swamps it in the unpaired view.
    ax = axes[2]
    for donor in pd.unique(donors):
        sel = np.flatnonzero(donors == donor)
        if len(sel) == 2:
            order = np.argsort(meta["arm"].to_numpy()[sel])
            pts = scores[sel][order]
            ax.plot(pts[:, 0], pts[:, 1], color="grey", lw=0.7, alpha=0.7, zorder=1)
    for lvl, colour in zip(["baseline", "treated"], ["#4C72B0", "#C44E52"]):
        sel = (meta["arm"] == lvl).to_numpy()
        ax.scatter(scores[sel, 0], scores[sel, 1], s=28, color=colour, label=lvl, zorder=2)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("Paired within-donor shift")
    ax.legend(frameon=False, fontsize=8)

    fig.tight_layout()
    (args.results / "figures").mkdir(parents=True, exist_ok=True)
    fig.savefig(args.results / "figures" / "ordination.png", dpi=160)

    pd.DataFrame(
        {
            "test": ["arm_free", "arm_within_donor", "batch_free"],
            "pseudo_F": [f_free, f_block, f_batch],
            "p_value": [p_free, p_block, p_batch],
            "valid_for_design": [False, True, True],
        }
    ).to_csv(args.results / "tables" / "permanova.tsv", sep="\t", index=False)

    z.to_csv(args.indir / "clr.tsv", sep="\t")
    print(f"\nfigure written to {args.results / 'figures' / 'ordination.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
