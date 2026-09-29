"""Contaminant identification following the decontam statistics.

Why this step exists
--------------------
DNA extraction kits, PCR reagents and water are not sterile. They contain
a reproducible community, usually Burkholderia, Ralstonia, Bradyrhizobium,
Pseudomonas, Cutibacterium and similar. In a stool study this is a small
fraction of the reads and often ignorable. In a low-biomass study such as
skin, lung lavage, placenta or a heavily antibiotic-treated gut, reagent
DNA can dominate and produce entirely spurious "findings" that correlate
with batch, because the contaminant load correlates with the extraction
kit lot rather than with biology.

A drug arm is the dangerous case. If treatment reduces bacterial biomass,
treated samples contain proportionally more reagent DNA, and contaminants
will appear "enriched by treatment" with very convincing statistics. Any
differential abundance analysis of a perturbation that plausibly changes
biomass should run this step and say so.

Two statistics
--------------
Prevalence method. Build the 2x2 table of presence and absence across
true samples and negative controls and test independence. A taxon that is
present in a large share of controls relative to its presence in true
samples is scored as a contaminant. Implemented here as a chi-square test
with a Fisher exact fallback for small tables, matching the behaviour of
``isContaminant(method="prevalence")``.

Frequency method. If a DNA concentration was measured per sample, a
contaminant's relative abundance should be inversely proportional to
total DNA, since it is a roughly fixed absolute amount diluted by a
variable amount of real template. The method fits, for each taxon, the
log-log regression of relative frequency on concentration and compares
the fit of a slope of -1, the contaminant expectation, against a slope of
0, the non-contaminant expectation.

Interpretation of the score
---------------------------
The score is a p value under the non-contaminant null, so small means
contaminant. The default threshold of 0.1 is deliberately permissive
compared with the usual 0.05 because the cost of retaining a contaminant
in a differential abundance analysis is higher than the cost of losing a
genuine low-prevalence taxon. Report the threshold and the number removed.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

__all__ = ["is_contaminant", "DecontamResult"]


@dataclass
class DecontamResult:
    table: pd.DataFrame
    threshold: float
    method: str

    @property
    def contaminants(self) -> pd.Index:
        return self.table.index[self.table["contaminant"]]

    def clean(self, counts: pd.DataFrame) -> pd.DataFrame:
        return counts.drop(columns=self.contaminants, errors="ignore")


def _prevalence_score(counts: pd.DataFrame, is_control: np.ndarray) -> np.ndarray:
    presence = (counts.to_numpy() > 0)
    n_ctrl = int(is_control.sum())
    n_samp = int((~is_control).sum())
    if n_ctrl == 0:
        raise ValueError("no negative controls supplied")

    scores = np.ones(presence.shape[1])
    for i in range(presence.shape[1]):
        a = int(presence[is_control, i].sum())          # present in control
        b = n_ctrl - a                                   # absent in control
        c = int(presence[~is_control, i].sum())          # present in sample
        d = n_samp - c
        table = np.array([[a, b], [c, d]])
        if table.min() < 0 or (a + c) == 0:
            scores[i] = 1.0
            continue
        if table.min() < 5:
            _, p = stats.fisher_exact(table, alternative="greater")
        else:
            chi2, p_two, _, _ = stats.chi2_contingency(table, correction=True)
            # one sided towards enrichment in controls
            enriched = (a / max(n_ctrl, 1)) > (c / max(n_samp, 1))
            p = p_two / 2 if enriched else 1 - p_two / 2
        scores[i] = float(np.clip(p, 0, 1))
    return scores


def _frequency_score(counts: pd.DataFrame, concentration: np.ndarray) -> np.ndarray:
    freq = counts.to_numpy(float)
    freq = freq / freq.sum(axis=1, keepdims=True)
    conc = np.asarray(concentration, float)
    if np.any(conc <= 0):
        raise ValueError("DNA concentrations must be positive")
    log_c = np.log(conc)

    scores = np.ones(freq.shape[1])
    for i in range(freq.shape[1]):
        obs = freq[:, i]
        keep = obs > 0
        if keep.sum() < 4:
            scores[i] = 1.0
            continue
        y = np.log(obs[keep])
        xc = log_c[keep]
        # Contaminant model: slope fixed at -1. Non contaminant: slope 0.
        resid_cont = y - (-1.0) * xc
        ss_cont = float(np.var(resid_cont) * keep.sum())
        ss_non = float(np.var(y) * keep.sum())
        df = max(keep.sum() - 1, 1)
        if ss_cont <= 0:
            scores[i] = 0.0
            continue
        f = ss_non / ss_cont
        scores[i] = float(stats.f.sf(f, df, df)) if f > 1 else 1.0
    return scores


def is_contaminant(counts: pd.DataFrame, is_control=None, concentration=None,
                   method: str = "prevalence", threshold: float = 0.1,
                   batch=None) -> DecontamResult:
    """Score every feature for being a reagent contaminant.

    Parameters
    ----------
    counts : samples as rows, including the negative control samples.
    is_control : boolean per sample, True for extraction or PCR blanks.
    concentration : per sample DNA quantitation, required for the
        frequency and combined methods.
    method : "prevalence", "frequency", "combined" or "either".
    batch : optional per sample batch label. Scores are computed within
        batch and combined by taking the minimum, which is the correct
        behaviour when each extraction round has its own kit lot and
        therefore its own contaminant profile.
    """
    if batch is not None:
        batch = pd.Series(list(batch), index=counts.index)
        per_batch = []
        for lvl in pd.unique(batch):
            sel = (batch == lvl).to_numpy()
            sub = counts.loc[sel]
            sub_ctrl = np.asarray(is_control)[sel] if is_control is not None else None
            sub_conc = np.asarray(concentration)[sel] if concentration is not None else None
            if sub_ctrl is not None and sub_ctrl.sum() == 0:
                continue
            res = is_contaminant(sub, sub_ctrl, sub_conc, method, threshold)
            per_batch.append(res.table["score"].rename(lvl))
        if not per_batch:
            raise ValueError("no batch contained a negative control")
        scores = pd.concat(per_batch, axis=1).min(axis=1).reindex(counts.columns).fillna(1.0)
        table = pd.DataFrame({"score": scores})
        table["contaminant"] = table["score"] < threshold
        return DecontamResult(table=table, threshold=threshold, method=method + "+batch")

    if method == "prevalence":
        score = _prevalence_score(counts, np.asarray(is_control, bool))
        parts = {"p_prev": score}
    elif method == "frequency":
        score = _frequency_score(counts, concentration)
        parts = {"p_freq": score}
    elif method in {"combined", "either"}:
        p_prev = _prevalence_score(counts, np.asarray(is_control, bool))
        p_freq = _frequency_score(counts, concentration)
        parts = {"p_prev": p_prev, "p_freq": p_freq}
        if method == "either":
            score = np.minimum(p_prev, p_freq)
        else:
            # Fisher's method on the two independent statistics.
            chi = -2 * (np.log(np.clip(p_prev, 1e-300, 1)) + np.log(np.clip(p_freq, 1e-300, 1)))
            score = stats.chi2.sf(chi, df=4)
    else:
        raise ValueError("unknown method")

    table = pd.DataFrame({**parts, "score": score}, index=counts.columns)
    table["prev_in_controls"] = (
        (counts.loc[np.asarray(is_control, bool)] > 0).mean(axis=0)
        if is_control is not None else np.nan
    )
    table["contaminant"] = table["score"] < threshold
    return DecontamResult(table=table, threshold=threshold, method=method)
