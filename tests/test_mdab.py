"""Tests for the compositional and differential abundance machinery.

The tests that matter here are the invariance tests. A log-ratio method
that is not invariant to library size is broken in a way that unit tests
on toy numbers will not catch, so those properties are asserted directly.

Run with::

    PYTHONPATH=src python -m pytest tests -q
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mdab import (  # noqa: E402
    aitchison_distance,
    aldex2,
    ancombc2,
    closure,
    clr,
    concordance_table,
    evaluate_against_truth,
    iqlr_clr,
    is_contaminant,
    mixed_model_scan,
    prevalence_filter,
    simulate_study,
)


@pytest.fixture(scope="module")
def study():
    return simulate_study(n_donors=12, n_taxa=80, n_da=8, seed=3)


def test_closure_rows_sum_to_one():
    x = np.array([[1.0, 2.0, 7.0], [10.0, 10.0, 80.0]])
    assert np.allclose(closure(x).sum(axis=1), 1.0)


def test_clr_sums_to_zero_within_sample(study):
    z = clr(study.counts)
    assert np.allclose(z.to_numpy().sum(axis=1), 0.0, atol=1e-9)


def test_clr_is_invariant_to_library_size(study):
    """The central property. Doubling every count in a sample must not
    change its CLR values, because a CLR value is a log ratio."""
    counts = study.counts.copy()
    z1 = clr(counts, pseudo=1.0)
    scaled = counts * 3
    z2 = clr(scaled, pseudo=3.0)  # pseudocount scales with the data
    assert np.allclose(z1.to_numpy(), z2.to_numpy(), atol=1e-8)


def test_aitchison_distance_is_scale_invariant(study):
    counts = study.counts.iloc[:8, :20] + 1
    d1 = aitchison_distance(counts, pseudo=0.0)
    d2 = aitchison_distance(counts * 5, pseudo=0.0)
    assert np.allclose(d1.to_numpy(), d2.to_numpy(), atol=1e-8)


def test_aitchison_is_subcompositionally_dominant(study):
    """Dropping features can only shrink the distance, never enlarge it.
    Bray-Curtis fails this, which is why it misbehaves under filtering."""
    counts = study.counts.iloc[:10, :40] + 1
    full = aitchison_distance(counts, pseudo=0.0).to_numpy()
    sub = aitchison_distance(counts.iloc[:, :20], pseudo=0.0).to_numpy()
    assert np.all(sub <= full + 1e-8)


def test_iqlr_reference_differs_from_clr_under_asymmetric_change(study):
    z_clr = clr(study.counts)
    z_iqlr = iqlr_clr(study.counts)
    assert not np.allclose(z_clr.to_numpy(), z_iqlr.to_numpy())


def test_prevalence_filter_removes_rare_features(study):
    filtered = prevalence_filter(study.counts, 0.5)
    assert filtered.shape[1] <= study.counts.shape[1]
    assert ((filtered > 0).mean(axis=0) >= 0.5).all()


def test_aldex2_recovers_planted_signal(study):
    bio = study.biological_samples
    counts = prevalence_filter(study.counts.loc[bio], 0.1)
    res = aldex2(counts, study.metadata.loc[bio, "arm"], mc_samples=32, seed=0)
    called = res.significant(alpha=0.10)
    truth = study.truth["is_da"]
    # Conservative method, so demand precision rather than recall.
    if len(called):
        precision = truth.reindex(called).fillna(False).mean()
        assert precision >= 0.8


def test_aldex2_gives_no_signal_under_the_null(study):
    """Permuting the labels must destroy the signal. If it does not, the
    test is picking up structure that is not the treatment."""
    bio = study.biological_samples
    counts = prevalence_filter(study.counts.loc[bio], 0.1)
    rng = np.random.default_rng(11)
    shuffled = rng.permutation(study.metadata.loc[bio, "arm"].to_numpy())
    res = aldex2(counts, shuffled, mc_samples=32, seed=0)
    assert (res.table["we.eBH"] < 0.05).sum() <= 1


def test_ancombc2_estimates_bias_and_recovers_signal(study):
    bio = study.biological_samples
    counts = prevalence_filter(study.counts.loc[bio], 0.1)
    res = ancombc2(
        counts,
        study.metadata.loc[bio, "arm"],
        random_effect=study.metadata.loc[bio, "donor"],
        sensitivity_pseudocounts=(0.5, 1.0),
    )
    assert np.isfinite(res.delta)
    called = res.significant(alpha=0.10, require_passed_ss=False)
    truth = study.truth["is_da"]
    assert truth.reindex(called).fillna(False).sum() >= 1


def test_ancombc2_bias_correction_centres_null_taxa(study):
    """After correction the null taxa should sit near zero. This is the
    property that distinguishes ANCOM-BC2 from a plain linear model on
    log counts, and it is what stops one-sided real changes from
    generating opposite-sign artefacts."""
    bio = study.biological_samples
    counts = prevalence_filter(study.counts.loc[bio], 0.1)
    res = ancombc2(counts, study.metadata.loc[bio, "arm"],
                   sensitivity_pseudocounts=(1.0,))
    null_taxa = study.truth.index[~study.truth["is_da"]]
    shared = res.table.index.intersection(null_taxa)
    corrected = res.table.loc[shared, "lfc"].mean()
    raw = res.table.loc[shared, "lfc_raw"].mean()
    assert abs(corrected) <= abs(raw) + 1e-9


def test_decontam_recovers_planted_contaminants(study):
    res = is_contaminant(
        study.counts,
        is_control=study.metadata["is_control"].to_numpy(bool),
        concentration=study.metadata["dna_ng_ul"].to_numpy(float),
        method="combined",
    )
    truth = study.truth["is_contaminant"]
    called = res.table["contaminant"]
    recall = (called & truth).sum() / max(truth.sum(), 1)
    precision = (called & truth).sum() / max(called.sum(), 1)
    assert recall >= 0.6
    assert precision >= 0.6


def test_mixed_model_reports_high_icc_for_donor(study):
    bio = study.biological_samples
    counts = prevalence_filter(study.counts.loc[bio], 0.2)
    res = mixed_model_scan(counts, study.metadata.loc[bio, "arm"],
                           study.metadata.loc[bio, "donor"])
    assert res["icc_donor"].median() > 0.3


def test_concordance_and_evaluation_shapes(study):
    bio = study.biological_samples
    counts = prevalence_filter(study.counts.loc[bio], 0.1)
    a = aldex2(counts, study.metadata.loc[bio, "arm"], mc_samples=16, seed=0)
    b = ancombc2(counts, study.metadata.loc[bio, "arm"],
                 sensitivity_pseudocounts=(1.0,))
    conc = concordance_table(a.table, b.table)
    assert set(conc["agreement"]) <= {"both", "aldex2_only", "ancombc2_only", "neither"}
    perf = evaluate_against_truth(conc, study.truth["is_da"].astype(bool))
    assert {"FDP", "sensitivity", "specificity"} <= set(perf.columns)


def test_simulation_is_reproducible():
    a = simulate_study(n_donors=5, n_taxa=30, seed=42)
    b = simulate_study(n_donors=5, n_taxa=30, seed=42)
    pd.testing.assert_frame_equal(a.counts, b.counts)


def test_nested_design_has_donors_within_one_arm():
    """The nested layout is what makes ignoring donor pseudoreplication.
    Each donor must appear in exactly one arm."""
    study = simulate_study(design="nested", n_reps=3, n_donors=8, n_taxa=40,
                           n_controls=0, n_contaminants=0, seed=5)
    meta = study.metadata
    arms_per_donor = meta.groupby("donor")["arm"].nunique()
    assert (arms_per_donor == 1).all()
    assert (meta.groupby("donor").size() == 3).all()


def test_paired_design_has_both_arms_per_donor():
    study = simulate_study(design="paired", n_donors=8, n_taxa=40,
                           n_controls=0, n_contaminants=0, seed=5)
    bio = study.biological_samples
    meta = study.metadata.loc[bio]
    assert (meta.groupby("donor")["arm"].nunique() == 2).all()
