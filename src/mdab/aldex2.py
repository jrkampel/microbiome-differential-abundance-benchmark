"""A transparent reimplementation of the ALDEx2 two-group test.

Why reimplement rather than only call the R package
---------------------------------------------------
The R package remains the reference implementation and ``R/05_aldex2.R``
in this repository calls it. This module exists so that the statistical
logic is inspectable, so the pipeline can run in environments without
Bioconductor, and so that the method comparison in ``scripts/`` can be
validated against simulated data with known ground truth.

The algorithm
-------------
1. The observed count vector of sample ``j`` is treated as one multinomial
   draw from an unknown composition. The posterior of that composition
   under a Jeffreys prior is

   .. math:: p_j \\mid n_j \\sim \\mathrm{Dirichlet}(n_j + 1/2)

   so a feature observed zero times has a posterior concentrated near but
   not at zero, with width set by the library size. This is the key
   design choice: zero handling and count uncertainty are the same
   operation, and both propagate into the test rather than being fixed
   before it.

2. Each Monte Carlo instance is CLR-transformed, using either the full
   geometric mean (``denom="all"``) or the interquartile variance subset
   (``denom="iqlr"``).

3. Each instance is tested independently. The reported p value is the
   expectation over instances, and the reported q value is the
   expectation of the per-instance Benjamini-Hochberg value. Averaging
   adjusted values instance-wise is what the package does and it is
   deliberately conservative.

4. The effect size is the median between-group difference divided by the
   larger within-group dispersion. It is a standardised quantity, not a
   fold change, and the package authors argue it is more reproducible
   across studies than a p value. An absolute effect above 1 is the
   conventional screening threshold.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

__all__ = ["aldex2", "AldexResult"]


@dataclass
class AldexResult:
    table: pd.DataFrame
    mc_samples: int
    denom: str

    def significant(self, alpha: float = 0.05, min_effect: float = 0.0,
                    column: str = "we.eBH") -> pd.Index:
        sel = (self.table[column] < alpha) & (self.table["effect"].abs() >= min_effect)
        return self.table.index[sel]


def _bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p values, monotone enforced."""
    p = np.asarray(p, float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def _effect_size(z: np.ndarray, a_idx: np.ndarray, b_idx: np.ndarray,
                 rng: np.random.Generator) -> np.ndarray:
    """Standardised effect for one Monte Carlo instance.

    ``diff.btw`` is the median of the between-group differences taken over
    random pairings of samples, ``diff.win`` is the larger of the two
    within-group median absolute differences. Dividing one by the other
    gives a dimensionless quantity in CLR units that does not depend on
    sample size the way a t statistic does.
    """
    za, zb = z[a_idx], z[b_idx]
    na, nb = len(a_idx), len(b_idx)
    n_pairs = max(na, nb)

    pick_a = rng.integers(0, na, n_pairs)
    pick_b = rng.integers(0, nb, n_pairs)
    diff_btw = np.median(zb[pick_b] - za[pick_a], axis=0)

    def within(block, n):
        i = rng.integers(0, n, n_pairs)
        j = rng.integers(0, n, n_pairs)
        return np.median(np.abs(block[i] - block[j]), axis=0)

    diff_win = np.maximum(within(za, na), within(zb, nb))
    diff_win = np.where(diff_win < 1e-12, 1e-12, diff_win)
    return diff_btw / diff_win


def aldex2(counts: pd.DataFrame, groups, mc_samples: int = 128,
           denom: str = "all", seed: int = 0,
           test: str = "both") -> AldexResult:
    """Run the ALDEx2 two-group test.

    Parameters
    ----------
    counts : DataFrame, samples as rows, features as columns, raw integer
        counts. Do not pass rarefied or relative data. Rarefying discards
        real information and does not solve compositionality; the Dirichlet
        step already models the sampling.
    groups : sequence of two distinct labels, one per sample. The second
        label encountered is treated as the numerator, so a positive
        effect means higher in that group.
    mc_samples : number of Dirichlet instances. 128 is the package
        default. Below about 64 the expected p values are visibly noisy.
    denom : "all" or "iqlr".

    Returns
    -------
    AldexResult with columns rab.all, rab.win.<g>, diff.btw, diff.win,
    effect, overlap, we.ep, we.eBH, wi.ep, wi.eBH.
    """
    rng = np.random.default_rng(seed)
    groups = pd.Series(list(groups), index=counts.index)
    levels = list(pd.unique(groups))
    if len(levels) != 2:
        raise ValueError(f"expected exactly two groups, found {levels}")
    ref, alt = levels
    a_idx = np.flatnonzero((groups == ref).to_numpy())
    b_idx = np.flatnonzero((groups == alt).to_numpy())

    n_feat = counts.shape[1]
    values = counts.to_numpy(dtype=float)

    we_p = np.zeros((mc_samples, n_feat))
    we_bh = np.zeros((mc_samples, n_feat))
    wi_p = np.zeros((mc_samples, n_feat))
    wi_bh = np.zeros((mc_samples, n_feat))
    effects = np.zeros((mc_samples, n_feat))
    btw = np.zeros((mc_samples, n_feat))
    win = np.zeros((mc_samples, n_feat))
    rab_all = np.zeros((mc_samples, n_feat))
    rab_a = np.zeros((mc_samples, n_feat))
    rab_b = np.zeros((mc_samples, n_feat))
    overlap = np.zeros((mc_samples, n_feat))

    for m in range(mc_samples):
        # Dirichlet posterior draw, one per sample, Jeffreys prior.
        inst = np.empty_like(values)
        for j in range(values.shape[0]):
            inst[j] = rng.gamma(values[j] + 0.5, 1.0)
            inst[j] /= inst[j].sum()

        z = _clr_from_props(inst, denom)

        t_stat, p_t = stats.ttest_ind(z[b_idx], z[a_idx], equal_var=False, axis=0)
        we_p[m] = np.nan_to_num(p_t, nan=1.0)
        we_bh[m] = _bh(we_p[m])

        u_stat, p_u = stats.mannwhitneyu(z[b_idx], z[a_idx], alternative="two-sided", axis=0)
        wi_p[m] = np.nan_to_num(p_u, nan=1.0)
        wi_bh[m] = _bh(wi_p[m])

        effects[m] = _effect_size(z, a_idx, b_idx, rng)
        btw[m] = np.median(z[b_idx], axis=0) - np.median(z[a_idx], axis=0)
        win[m] = np.maximum(_mad(z[a_idx]), _mad(z[b_idx]))
        rab_all[m] = np.median(z, axis=0)
        rab_a[m] = np.median(z[a_idx], axis=0)
        rab_b[m] = np.median(z[b_idx], axis=0)
        overlap[m] = _overlap(z, a_idx, b_idx)

    table = pd.DataFrame(
        {
            "rab.all": rab_all.mean(axis=0),
            f"rab.win.{ref}": rab_a.mean(axis=0),
            f"rab.win.{alt}": rab_b.mean(axis=0),
            "diff.btw": np.median(btw, axis=0),
            "diff.win": np.median(win, axis=0),
            "effect": np.median(effects, axis=0),
            "overlap": overlap.mean(axis=0),
            "we.ep": we_p.mean(axis=0),
            "we.eBH": we_bh.mean(axis=0),
            "wi.ep": wi_p.mean(axis=0),
            "wi.eBH": wi_bh.mean(axis=0),
        },
        index=counts.columns,
    )
    if test == "welch":
        table = table.drop(columns=["wi.ep", "wi.eBH"])
    elif test == "wilcoxon":
        table = table.drop(columns=["we.ep", "we.eBH"])
    return AldexResult(table=table, mc_samples=mc_samples, denom=denom)


def _clr_from_props(props: np.ndarray, denom: str) -> np.ndarray:
    log_p = np.log(props)
    if denom == "all":
        return log_p - log_p.mean(axis=1, keepdims=True)
    if denom == "iqlr":
        base = log_p - log_p.mean(axis=1, keepdims=True)
        var = base.var(axis=0, ddof=1)
        lo, hi = np.percentile(var, [25, 75])
        keep = (var >= lo) & (var <= hi)
        if keep.sum() < 2:
            keep = np.ones_like(var, dtype=bool)
        return log_p - log_p[:, keep].mean(axis=1, keepdims=True)
    raise ValueError("denom must be 'all' or 'iqlr'")


def _mad(block: np.ndarray) -> np.ndarray:
    med = np.median(block, axis=0)
    return np.median(np.abs(block - med), axis=0)


def _overlap(z: np.ndarray, a_idx: np.ndarray, b_idx: np.ndarray) -> np.ndarray:
    """Proportion of the two CLR distributions that overlap.

    A non-parametric companion to the effect size. Values near 0.5 mean
    the groups are indistinguishable for that feature regardless of what
    the p value says.
    """
    za, zb = z[a_idx], z[b_idx]
    mid = 0.5 * (np.median(za, axis=0) + np.median(zb, axis=0))
    below_b = (zb < mid).mean(axis=0)
    above_a = (za > mid).mean(axis=0)
    return np.minimum(below_b + above_a, 1.0)
