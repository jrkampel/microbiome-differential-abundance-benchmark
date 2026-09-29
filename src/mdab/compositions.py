"""Compositional transforms for amplicon count tables.

Conventions used throughout this package
----------------------------------------
Count tables are stored as pandas DataFrames with samples as rows and
features (ASVs or collapsed taxa) as columns. This matches the
scikit-learn / statsmodels orientation and is the transpose of the
phyloseq otu_table default, so conversion scripts transpose explicitly.

The central problem these functions address is that sequencing returns a
fixed number of reads per sample. The observed vector is therefore a
composition, carrying only relative information, and lives on the simplex
rather than in real space. Log-ratio transforms map the simplex into real
space so that ordinary linear methods become admissible.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = [
    "closure",
    "multiplicative_replacement",
    "czm_replacement",
    "clr",
    "alr",
    "iqlr_clr",
    "aitchison_distance",
    "prevalence_filter",
]


def closure(counts: np.ndarray) -> np.ndarray:
    """Rescale each row to sum to one.

    Parameters
    ----------
    counts : array of shape (n_samples, n_features), non-negative.

    Notes
    -----
    Closure is the operation that makes the constant-sum constraint
    explicit. It is idempotent and it is the reason correlations between
    raw relative abundances are spurious: forcing rows to sum to one
    induces negative dependence between features even when the underlying
    absolute abundances are independent.
    """
    counts = np.asarray(counts, dtype=float)
    totals = counts.sum(axis=1, keepdims=True)
    if np.any(totals <= 0):
        raise ValueError("every sample must have at least one read")
    return counts / totals


def multiplicative_replacement(counts: np.ndarray, delta: float | None = None) -> np.ndarray:
    """Replace zeros by a small value and rescale the non-zero parts.

    This is the simplest zero replacement that preserves ratios among the
    observed parts. Each zero becomes ``delta`` and the non-zero parts are
    multiplied by ``(1 - n_zeros * delta)`` so the row still sums to one.

    ``delta`` defaults to ``0.65 / N`` where ``N`` is the row total, the
    value recommended by Martin-Fernandez and colleagues. The choice is a
    modelling assumption, not a neutral default, which is exactly why
    ANCOM-BC2 runs a pseudocount sensitivity analysis and why ALDEx2
    avoids a point replacement altogether.
    """
    counts = np.asarray(counts, dtype=float)
    out = np.empty_like(counts)
    for i, row in enumerate(counts):
        total = row.sum()
        d = delta if delta is not None else 0.65 / total
        prop = row / total
        zeros = prop == 0
        n_zero = int(zeros.sum())
        if n_zero == 0:
            out[i] = prop
            continue
        prop = prop * (1.0 - n_zero * d)
        prop[zeros] = d
        out[i] = prop
    return out


def czm_replacement(counts: np.ndarray, frac: float = 0.65) -> np.ndarray:
    """Count zero multiplicative replacement, the zCompositions default.

    Identical in form to :func:`multiplicative_replacement` but the
    replacement value scales with the per-sample detection limit
    ``frac / N_j``, so deeply sequenced samples receive a smaller
    replacement. This matters because a fixed pseudocount of 1 is a much
    larger perturbation in a shallow sample than in a deep one and
    therefore couples zero handling to library size.
    """
    counts = np.asarray(counts, dtype=float)
    out = np.empty_like(counts)
    for i, row in enumerate(counts):
        total = row.sum()
        out[i] = multiplicative_replacement(row[None, :], delta=frac / total)[0]
    return out


def clr(counts, pseudo: str | float = "czm", return_frame: bool = True):
    r"""Centred log-ratio transform.

    .. math::
        \mathrm{clr}(x_{ij}) = \log x_{ij} - \frac{1}{D}\sum_{k=1}^{D}\log x_{kj}

    The reference is the geometric mean of the sample, so CLR values are
    log ratios against an internal reference and are invariant to library
    size. Two consequences are worth stating in any write-up:

    1. CLR values sum to zero within a sample by construction, so the
       features are linearly dependent and the covariance matrix is
       singular. Methods requiring full rank need ILR or a subcomposition.
    2. The reference moves if the community composition changes globally.
       Under a strong, one-directional perturbation the geometric mean
       shifts and unchanged taxa acquire apparent changes of the opposite
       sign. The IQLR variant below reduces this.

    Parameters
    ----------
    pseudo : "czm", "mult", or a float added to every count.
    """
    is_frame = isinstance(counts, pd.DataFrame)
    values = counts.to_numpy(dtype=float) if is_frame else np.asarray(counts, float)

    if pseudo == "czm":
        props = czm_replacement(values)
    elif pseudo == "mult":
        props = multiplicative_replacement(values)
    else:
        props = closure(values + float(pseudo))

    log_p = np.log(props)
    out = log_p - log_p.mean(axis=1, keepdims=True)

    if is_frame and return_frame:
        return pd.DataFrame(out, index=counts.index, columns=counts.columns)
    return out


def alr(counts, reference: int | str, pseudo: str | float = "czm"):
    """Additive log-ratio against one chosen feature.

    Useful when a genuine internal reference exists, for example a spike-in
    or a taxon known to be unaffected by the perturbation. Unlike CLR the
    result is full rank, but every value is conditional on the reference
    being stable, which is a strong claim that should be justified.
    """
    is_frame = isinstance(counts, pd.DataFrame)
    values = counts.to_numpy(dtype=float) if is_frame else np.asarray(counts, float)
    props = czm_replacement(values) if pseudo == "czm" else closure(values + float(pseudo))
    if isinstance(reference, str):
        if not is_frame:
            raise TypeError("string reference requires a DataFrame")
        ref_idx = list(counts.columns).index(reference)
    else:
        ref_idx = reference
    log_p = np.log(props)
    out = log_p - log_p[:, [ref_idx]]
    out = np.delete(out, ref_idx, axis=1)
    if is_frame:
        cols = [c for k, c in enumerate(counts.columns) if k != ref_idx]
        return pd.DataFrame(out, index=counts.index, columns=cols)
    return out


def iqlr_clr(counts, pseudo: str | float = "czm"):
    """CLR using only interquartile-variance features as the reference set.

    ALDEx2 offers this as ``denom="iqlr"``. The denominator is restricted
    to features whose CLR variance falls in the interquartile range, that
    is, features that are neither invariant nor wildly variable. When a
    large fraction of the community is asymmetrically perturbed, the
    ordinary geometric mean reference drifts and produces sign-flipped
    artefacts; restricting the reference to moderately variable features
    stabilises it.
    """
    is_frame = isinstance(counts, pd.DataFrame)
    values = counts.to_numpy(dtype=float) if is_frame else np.asarray(counts, float)
    base = clr(values, pseudo=pseudo, return_frame=False)
    var = base.var(axis=0, ddof=1)
    lo, hi = np.percentile(var, [25, 75])
    keep = (var >= lo) & (var <= hi)
    if keep.sum() < 2:
        keep = np.ones_like(var, dtype=bool)

    props = czm_replacement(values) if pseudo == "czm" else closure(values + float(pseudo))
    log_p = np.log(props)
    ref = log_p[:, keep].mean(axis=1, keepdims=True)
    out = log_p - ref
    if is_frame:
        return pd.DataFrame(out, index=counts.index, columns=counts.columns)
    return out


def aitchison_distance(counts, pseudo: str | float = "czm") -> pd.DataFrame:
    """Euclidean distance in CLR space.

    This is the only distance on compositions that is simultaneously
    scale invariant, subcompositionally dominant and permutation
    invariant. Bray-Curtis and unweighted UniFrac satisfy none of the
    three, which is why ordinations built on them can shift when a single
    dominant taxon is added or removed.
    """
    z = clr(counts, pseudo=pseudo, return_frame=False)
    diff = z[:, None, :] - z[None, :, :]
    d = np.sqrt((diff ** 2).sum(axis=-1))
    idx = counts.index if isinstance(counts, pd.DataFrame) else None
    return pd.DataFrame(d, index=idx, columns=idx)


def prevalence_filter(counts: pd.DataFrame, min_prevalence: float = 0.10,
                      min_count: int = 1) -> pd.DataFrame:
    """Drop features seen in fewer than ``min_prevalence`` of samples.

    Filtering is unavoidable because the vast majority of ASVs appear in a
    handful of samples and contribute nothing but multiple-testing burden.
    It is also a decision that changes results, because it changes the CLR
    reference. Report the threshold, and check that conclusions survive a
    reasonable range of it.
    """
    present = (counts >= min_count).mean(axis=0)
    keep = present[present >= min_prevalence].index
    return counts.loc[:, keep]
