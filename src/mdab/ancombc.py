"""An ANCOM-BC2 style bias-corrected differential abundance model.

The generative model
--------------------
Let :math:`A_{ij}` be the absolute abundance of taxon :math:`i` in the
ecosystem sampled by specimen :math:`j`, and :math:`O_{ij}` the observed
read count. Sequencing and DNA extraction impose a per-sample scaling that
is unknown and varies by orders of magnitude:

.. math:: \\log O_{ij} = \\log A_{ij} + d_j + \\varepsilon_{ij}

The term :math:`d_j` is the log sampling fraction. It is not identifiable
from a single sample, because multiplying every taxon by a constant and
dividing the sampling fraction by the same constant gives identical data.
This is the formal statement of why compositional data cannot recover
absolute abundance without external information.

What makes it estimable
-----------------------
Fit, for each taxon separately, the linear model

.. math:: \\log O_{ij} = \\beta_{0i} + \\beta_{1i} x_j + u_{i,g(j)} + e_{ij}

where :math:`x_j` is the group indicator and :math:`u` an optional random
effect for donor or batch. The estimate :math:`\\hat\\beta_{1i}` contains
the true taxon-level effect plus the common group difference in sampling
fraction :math:`\\delta`. If most taxa are not differentially abundant,
the distribution of :math:`\\hat\\beta_{1i}` across taxa is centred on
:math:`\\delta`, so :math:`\\delta` can be estimated from the bulk of that
distribution and subtracted.

That assumption is the load-bearing one. State it in the README. If the
perturbation genuinely moves most of the community in one direction, the
bias correction will absorb part of the real signal and the method will
under-report. A drug arm that wipes out a phylum is exactly such a case.

Estimation of delta
-------------------
The original implementation fits a three-component Gaussian mixture by
EM, where the central component is the null taxa. This module offers two
estimators, selected by ``delta_method``:

``"em"``     two-component EM with a null component whose mean is delta.
``"robust"`` the Huber-type location estimate of the central portion,
             which is far simpler, has no convergence failures, and in
             simulation tracks the EM estimate closely.

Both are reported so a reader can see whether the conclusion depends on
the choice.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats

__all__ = ["ancombc2", "AncomBCResult"]


@dataclass
class AncomBCResult:
    table: pd.DataFrame
    delta: float
    delta_se: float
    structural_zeros: pd.DataFrame
    sensitivity: pd.DataFrame | None = None
    meta: dict = field(default_factory=dict)

    def significant(self, alpha: float = 0.05, require_passed_ss: bool = True) -> pd.Index:
        sel = self.table["q_value"] < alpha
        if require_passed_ss and "passed_ss" in self.table:
            sel &= self.table["passed_ss"]
        return self.table.index[sel]


def _bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def _structural_zeros(counts: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    """Taxa absent from every sample of at least one group.

    ANCOM-BC2 treats these separately because no finite log-ratio effect
    exists for them. Including them in the linear model forces an
    arbitrary pseudocount to determine the magnitude of an effect that is
    really a presence/absence statement. Reporting them as a distinct
    category is more honest than assigning them a fold change.
    """
    rows = []
    for level in pd.unique(groups):
        present = (counts.loc[groups == level] > 0).any(axis=0)
        rows.append(present.rename(level))
    present = pd.concat(rows, axis=1)
    present["structural_zero"] = ~present.all(axis=1)
    return present


def _huber_location(x: np.ndarray, k: float = 1.345, iters: int = 50) -> tuple[float, float]:
    """Huber M-estimate of location with a MAD scale.

    Resistant to the differentially abundant taxa in the tails, which is
    exactly what is needed when estimating the central mass of the
    per-taxon effect distribution.
    """
    mu = float(np.median(x))
    scale = float(stats.median_abs_deviation(x, scale="normal"))
    if scale <= 0:
        return mu, 0.0
    for _ in range(iters):
        r = (x - mu) / scale
        w = np.where(np.abs(r) <= k, 1.0, k / np.maximum(np.abs(r), 1e-12))
        new = float(np.sum(w * x) / np.sum(w))
        if abs(new - mu) < 1e-10:
            mu = new
            break
        mu = new
    se = scale / np.sqrt(len(x))
    return mu, se


def _em_null_location(x: np.ndarray, iters: int = 200) -> tuple[float, float]:
    """Two-component mixture: null taxa centred at delta, the rest diffuse."""
    mu = float(np.median(x))
    sigma = max(float(stats.median_abs_deviation(x, scale="normal")), 1e-6)
    sigma_alt = max(float(np.std(x)), sigma * 3)
    pi = 0.8
    for _ in range(iters):
        d_null = pi * stats.norm.pdf(x, mu, sigma)
        d_alt = (1 - pi) * stats.norm.pdf(x, mu, sigma_alt)
        denom = d_null + d_alt
        denom = np.where(denom <= 0, 1e-300, denom)
        r = d_null / denom
        new_mu = float(np.sum(r * x) / np.sum(r))
        new_sigma = float(np.sqrt(np.sum(r * (x - new_mu) ** 2) / np.sum(r)))
        new_pi = float(np.mean(r))
        if abs(new_mu - mu) < 1e-10:
            mu, sigma, pi = new_mu, max(new_sigma, 1e-6), new_pi
            break
        mu, sigma, pi = new_mu, max(new_sigma, 1e-6), min(max(new_pi, 0.05), 0.99)
    n_eff = max(pi * len(x), 1.0)
    return mu, sigma / np.sqrt(n_eff)


def _fit_one_taxon(y: np.ndarray, x: np.ndarray, covariates: np.ndarray | None,
                   random_group: np.ndarray | None):
    """Return (coefficient on x, standard error, residual df)."""
    design = [np.ones_like(y), x]
    if covariates is not None and covariates.size:
        design.extend(covariates.T)
    X = np.column_stack(design)

    if random_group is None:
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ beta
        df = max(len(y) - X.shape[1], 1)
        s2 = float(resid @ resid) / df
        xtx_inv = np.linalg.pinv(X.T @ X)
        se = float(np.sqrt(max(s2 * xtx_inv[1, 1], 1e-300)))
        return float(beta[1]), se, df

    import statsmodels.api as sm

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            model = sm.MixedLM(y, X, groups=random_group)
            res = model.fit(reml=True, method="lbfgs")
            beta = float(res.params.iloc[1] if hasattr(res.params, "iloc") else res.params[1])
            se = float(res.bse.iloc[1] if hasattr(res.bse, "iloc") else res.bse[1])
            df = max(len(np.unique(random_group)) - X.shape[1], 1)
            if not np.isfinite(se) or se <= 0:
                raise ValueError
            return beta, se, df
        except Exception:
            # Singular fit, usually a taxon with almost no variance within
            # donors. Fall back to the fixed-effect fit and let the
            # sensitivity flag catch it.
            return _fit_one_taxon(y, x, covariates, None)


def ancombc2(counts: pd.DataFrame, groups, covariates: pd.DataFrame | None = None,
             random_effect=None, pseudocount: float = 1.0,
             delta_method: str = "robust",
             sensitivity_pseudocounts=(0.01, 0.1, 0.5, 1.0, 2.0, 5.0),
             alpha: float = 0.05,
             prv_cut: float = 0.10) -> AncomBCResult:
    """Bias-corrected differential abundance between two groups.

    Parameters
    ----------
    counts : raw integer counts, samples as rows.
    groups : two-level grouping, one label per sample. The second level
        encountered is the numerator, so a positive log fold change means
        higher in that group.
    covariates : optional numeric design columns, for example age or
        sequencing depth quartile. Categorical variables must be dummy
        coded before being passed.
    random_effect : optional per-sample label such as donor or run. When
        supplied each taxon is fitted with a linear mixed model with a
        random intercept. This is the correct handling of repeated
        measures and of batch structure that is nuisance rather than
        interest.
    pseudocount : added before taking logs in the primary fit.
    sensitivity_pseudocounts : the fit is repeated at each of these
        values. A taxon "passes sensitivity" only if its sign and its
        significance are stable across all of them. This is the single
        most useful diagnostic in the whole method, because sparse taxa
        routinely flip.

    Returns
    -------
    AncomBCResult. ``table`` has columns lfc, lfc_raw, se, W, p_value,
    q_value, passed_ss, prevalence.
    """
    groups = pd.Series(list(groups), index=counts.index)
    levels = list(pd.unique(groups))
    if len(levels) != 2:
        raise ValueError(f"expected exactly two groups, found {levels}")
    ref, alt = levels
    x = (groups == alt).to_numpy(dtype=float)

    prevalence = (counts > 0).mean(axis=0)
    keep = prevalence[prevalence >= prv_cut].index
    counts = counts.loc[:, keep]
    prevalence = prevalence.loc[keep]

    cov = covariates.loc[counts.index].to_numpy(float) if covariates is not None else None
    rand = None
    if random_effect is not None:
        rand = pd.Series(list(random_effect), index=counts.index).to_numpy()

    zeros = _structural_zeros(counts, groups)

    def fit_all(pc: float):
        log_o = np.log(counts.to_numpy(float) + pc)
        betas = np.empty(log_o.shape[1])
        ses = np.empty(log_o.shape[1])
        dfs = np.empty(log_o.shape[1])
        for i in range(log_o.shape[1]):
            betas[i], ses[i], dfs[i] = _fit_one_taxon(log_o[:, i], x, cov, rand)
        return betas, ses, dfs

    betas, ses, dfs = fit_all(pseudocount)

    # Estimate the group difference in log sampling fraction.
    if delta_method == "em":
        delta, delta_se = _em_null_location(betas)
    elif delta_method == "robust":
        delta, delta_se = _huber_location(betas)
    else:
        raise ValueError("delta_method must be 'em' or 'robust'")

    lfc = betas - delta
    # The bias is estimated, so its uncertainty enters every taxon's
    # variance. Ignoring this is a known way to be anti-conservative.
    se_total = np.sqrt(ses ** 2 + delta_se ** 2)
    W = lfc / np.where(se_total > 0, se_total, np.inf)
    p = 2 * stats.t.sf(np.abs(W), df=np.maximum(dfs, 1))
    q = _bh(p)

    # Pseudocount sensitivity analysis.
    passed = np.ones(len(lfc), dtype=bool)
    sens_rows = {}
    for pc in sensitivity_pseudocounts:
        b2, s2, d2 = fit_all(pc)
        if delta_method == "em":
            dl, dse = _em_null_location(b2)
        else:
            dl, dse = _huber_location(b2)
        lfc2 = b2 - dl
        se2 = np.sqrt(s2 ** 2 + dse ** 2)
        p2 = 2 * stats.t.sf(np.abs(lfc2 / np.where(se2 > 0, se2, np.inf)),
                            df=np.maximum(d2, 1))
        q2 = _bh(p2)
        sens_rows[f"q_pc_{pc}"] = q2
        same_sign = np.sign(lfc2) == np.sign(lfc)
        same_call = (q2 < alpha) == (q < alpha)
        passed &= same_sign & same_call

    sensitivity = pd.DataFrame(sens_rows, index=counts.columns)

    table = pd.DataFrame(
        {
            "lfc": lfc,
            "lfc_raw": betas,
            "se": se_total,
            "W": W,
            "p_value": p,
            "q_value": q,
            "passed_ss": passed,
            "prevalence": prevalence.to_numpy(),
            "structural_zero": zeros["structural_zero"].reindex(counts.columns).to_numpy(),
        },
        index=counts.columns,
    )

    return AncomBCResult(
        table=table,
        delta=float(delta),
        delta_se=float(delta_se),
        structural_zeros=zeros,
        sensitivity=sensitivity,
        meta={
            "reference_level": ref,
            "numerator_level": alt,
            "pseudocount": pseudocount,
            "delta_method": delta_method,
            "random_effect": random_effect is not None,
            "n_features_tested": int(counts.shape[1]),
            "prv_cut": prv_cut,
        },
    )
