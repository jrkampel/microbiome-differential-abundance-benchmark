"""Compare differential abundance callers and, when truth is known, score them.

The comparison is the point of the repository. Two defensible methods run
on the same data will not return the same list, and the interesting
question is not which one is right but where and why they differ. Three
summaries are produced:

``concordance_table``
    One row per feature with each method's effect estimate, q value and
    call, plus an agreement label. This is the artefact to inspect by
    hand.

``agreement_summary``
    Counts of the four agreement categories and the rank correlation of
    the effect estimates. Rank correlation is the right statistic here
    because the two methods report effects on different scales, a
    standardised CLR difference against a bias-corrected log fold change.

``evaluate_against_truth``
    Realised false discovery proportion, sensitivity and specificity when
    a ground truth column is available, which in practice means simulated
    data or a spike-in study.

A caution to keep in the README: agreement between two methods is not
evidence of correctness. Both use log ratios and both inherit the same
blind spot, namely that a uniform change in absolute biomass is invisible.
Two methods can agree perfectly and both be wrong about the biology.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

__all__ = ["concordance_table", "agreement_summary", "evaluate_against_truth"]


def concordance_table(aldex_table: pd.DataFrame, ancom_table: pd.DataFrame,
                      alpha: float = 0.05, min_effect: float = 0.0,
                      require_ss: bool = True) -> pd.DataFrame:
    """Align the two result tables on shared features and label agreement."""
    shared = aldex_table.index.intersection(ancom_table.index)
    a = aldex_table.loc[shared]
    b = ancom_table.loc[shared]

    aldex_call = (a["we.eBH"] < alpha) & (a["effect"].abs() >= min_effect)
    ancom_call = b["q_value"] < alpha
    if require_ss and "passed_ss" in b:
        ancom_call = ancom_call & b["passed_ss"]

    label = np.where(
        aldex_call & ancom_call, "both",
        np.where(aldex_call, "aldex2_only",
                 np.where(ancom_call, "ancombc2_only", "neither")),
    )

    out = pd.DataFrame(
        {
            "aldex_effect": a["effect"],
            "aldex_diff_btw": a["diff.btw"],
            "aldex_q": a["we.eBH"],
            "aldex_call": aldex_call,
            "ancom_lfc": b["lfc"],
            "ancom_q": b["q_value"],
            "ancom_passed_ss": b["passed_ss"] if "passed_ss" in b else True,
            "ancom_call": ancom_call,
            "prevalence": b["prevalence"] if "prevalence" in b else np.nan,
            "agreement": label,
        },
        index=shared,
    )
    return out.sort_values("ancom_q")


def agreement_summary(conc: pd.DataFrame) -> pd.Series:
    counts = conc["agreement"].value_counts()
    n_a = int(conc["aldex_call"].sum())
    n_b = int(conc["ancom_call"].sum())
    both = int((conc["agreement"] == "both").sum())
    union = n_a + n_b - both
    rho, p_rho = stats.spearmanr(conc["aldex_effect"], conc["ancom_lfc"])
    return pd.Series(
        {
            "n_features": len(conc),
            "n_aldex2": n_a,
            "n_ancombc2": n_b,
            "n_both": both,
            "n_aldex2_only": int(counts.get("aldex2_only", 0)),
            "n_ancombc2_only": int(counts.get("ancombc2_only", 0)),
            "jaccard": both / union if union else np.nan,
            "spearman_effect_vs_lfc": rho,
            "spearman_p": p_rho,
        }
    )


def evaluate_against_truth(conc: pd.DataFrame, truth: pd.Series) -> pd.DataFrame:
    """Realised FDP, sensitivity and specificity per method.

    ``truth`` is a boolean Series indexed like ``conc``. The false
    discovery proportion reported here is the realised proportion in this
    one dataset, not an expectation, so it will fluctuate around the
    nominal level. Averaging over repeated simulations is the correct way
    to assess calibration and ``scripts/06_simulation_study.py`` does that.
    """
    t = truth.reindex(conc.index).fillna(False).to_numpy(bool)
    rows = []
    for name, call in (("ALDEx2", conc["aldex_call"]), ("ANCOM-BC2", conc["ancom_call"]),
                       ("intersection", conc["agreement"] == "both"),
                       ("union", conc["aldex_call"] | conc["ancom_call"])):
        c = call.to_numpy(bool)
        tp = int((c & t).sum())
        fp = int((c & ~t).sum())
        fn = int((~c & t).sum())
        tn = int((~c & ~t).sum())
        rows.append(
            {
                "method": name,
                "TP": tp, "FP": fp, "FN": fn, "TN": tn,
                "n_called": tp + fp,
                "FDP": fp / (tp + fp) if (tp + fp) else 0.0,
                "sensitivity": tp / (tp + fn) if (tp + fn) else np.nan,
                "specificity": tn / (tn + fp) if (tn + fp) else np.nan,
            }
        )
    return pd.DataFrame(rows).set_index("method")
