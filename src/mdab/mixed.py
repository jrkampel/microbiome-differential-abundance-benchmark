"""Per-taxon linear mixed models on CLR values.

The statistical point
---------------------
Two samples from the same donor are more alike than two samples from
different donors. Treating them as independent replicates inflates the
effective sample size, shrinks the standard error and produces p values
that are too small. In a paired design the inflation can be severe,
because between-donor variance in gut composition routinely exceeds the
treatment effect by a factor of several.

The model fitted for each taxon is

.. math::
    z_{ij} = \\beta_0 + \\beta_1 \\mathrm{treat}_j + u_{d(j)} + e_{ij},
    \\quad u_d \\sim N(0, \\sigma_u^2), \\ e_{ij} \\sim N(0, \\sigma^2)

where :math:`z` is the CLR value. The intraclass correlation

.. math:: \\mathrm{ICC} = \\sigma_u^2 / (\\sigma_u^2 + \\sigma^2)

is reported per taxon. It is the fraction of variance attributable to
donor identity and is the quantity to cite when arguing that donor had to
be modelled. Values above roughly 0.5 are common for gut taxa.

Why CLR and not raw counts here
-------------------------------
A Gaussian mixed model needs an approximately symmetric, unbounded
response with stable variance. CLR values satisfy this far better than
counts, proportions or log proportions, and the log-ratio is what the
compositional argument says is the estimable quantity in the first place.
The alternative, a negative binomial GLMM on counts with an offset, is
defensible but couples the mean-variance model to the library size
normalisation in a way that is harder to defend when the offset is
itself compositional.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy import stats

from .compositions import clr

__all__ = ["mixed_model_scan", "variance_partition"]


def _bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def mixed_model_scan(counts: pd.DataFrame, groups, random_effect,
                     covariates: pd.DataFrame | None = None,
                     pseudo: str | float = "czm") -> pd.DataFrame:
    """Fit ``clr ~ group + (1 | random_effect)`` for every feature.

    Returns a table with the fixed-effect estimate, its standard error,
    the p value, the BH q value, the ICC, and the naive fixed-effect-only
    p value for comparison. The gap between the two p value columns is
    the cost of ignoring the design.
    """
    import statsmodels.api as sm

    z = clr(counts, pseudo=pseudo)
    groups = pd.Series(list(groups), index=counts.index)
    levels = list(pd.unique(groups))
    if len(levels) != 2:
        raise ValueError("mixed_model_scan expects two groups")
    x = (groups == levels[1]).to_numpy(float)
    rand = pd.Series(list(random_effect), index=counts.index).to_numpy()

    design = [np.ones(len(x)), x]
    if covariates is not None:
        design.extend(covariates.loc[counts.index].to_numpy(float).T)
    X = np.column_stack(design)

    rows = []
    for feature in z.columns:
        y = z[feature].to_numpy(float)

        beta_ols, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ beta_ols
        df_ols = max(len(y) - X.shape[1], 1)
        s2 = float(resid @ resid) / df_ols
        se_ols = float(np.sqrt(max(s2 * np.linalg.pinv(X.T @ X)[1, 1], 1e-300)))
        p_ols = 2 * stats.t.sf(abs(beta_ols[1] / se_ols), df_ols)

        beta, se, icc, p_mixed, converged = np.nan, np.nan, np.nan, np.nan, False
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                res = sm.MixedLM(y, X, groups=rand).fit(reml=True, method="lbfgs")
                beta = float(np.asarray(res.params)[1])
                se = float(np.asarray(res.bse)[1])
                sigma_u = float(np.asarray(res.cov_re)[0, 0])
                sigma_e = float(res.scale)
                icc = sigma_u / (sigma_u + sigma_e) if (sigma_u + sigma_e) > 0 else np.nan
                df_mixed = max(len(np.unique(rand)) - X.shape[1], 1)
                p_mixed = 2 * stats.t.sf(abs(beta / se), df_mixed) if se > 0 else 1.0
                converged = bool(res.converged)
            except Exception:
                pass

        rows.append(
            {
                "feature": feature,
                "beta_mixed": beta,
                "se_mixed": se,
                "p_mixed": p_mixed,
                "icc_donor": icc,
                "converged": converged,
                "beta_ols": float(beta_ols[1]),
                "se_ols": se_ols,
                "p_ols": float(p_ols),
            }
        )

    out = pd.DataFrame(rows).set_index("feature")
    out["q_mixed"] = _bh(out["p_mixed"].fillna(1.0).to_numpy())
    out["q_ols"] = _bh(out["p_ols"].to_numpy())
    return out


def variance_partition(counts: pd.DataFrame, factors: pd.DataFrame,
                       pseudo: str | float = "czm") -> pd.DataFrame:
    """Share of CLR variance explained by each categorical factor.

    A blunt but useful diagnostic to put in a README figure: if batch
    explains more variance than the treatment arm, no amount of clever
    testing rescues the study, and the honest conclusion is that the
    design is confounded.
    """
    z = clr(counts, pseudo=pseudo)
    out = {}
    for name in factors.columns:
        labels = factors.loc[z.index, name]
        grand = z.mean(axis=0)
        total = ((z - grand) ** 2).sum(axis=0)
        between = 0.0
        for lvl, block in z.groupby(labels.to_numpy()):
            between = between + len(block) * (block.mean(axis=0) - grand) ** 2
        out[name] = (between / total.replace(0, np.nan)).to_numpy()
    return pd.DataFrame(out, index=z.columns)
