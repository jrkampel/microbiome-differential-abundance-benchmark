"""mdab: microbiome differential abundance benchmarking.

A small, readable implementation of the compositional workflow used for
16S amplicon data, together with a simulator that supplies ground truth
so the methods can be scored rather than merely compared.

Scope and honesty statement
---------------------------
This package does not replace DADA2, ALDEx2, ANCOM-BC2 or decontam. Those
are the reference implementations and the R scripts in ``R/`` call them.
What lives here is a second, independent implementation of the same
statistics, written to be read. Agreement between the two implementations
on the same data is checked in ``tests/``.
"""

from .aldex2 import AldexResult, aldex2
from .ancombc import AncomBCResult, ancombc2
from .compare import agreement_summary, concordance_table, evaluate_against_truth
from .compositions import (
    aitchison_distance,
    alr,
    closure,
    clr,
    czm_replacement,
    iqlr_clr,
    multiplicative_replacement,
    prevalence_filter,
)
from .decontam import DecontamResult, is_contaminant
from .mixed import mixed_model_scan, variance_partition
from .simulate import SimulatedStudy, simulate_study

__version__ = "0.1.0"

__all__ = [
    "clr", "alr", "iqlr_clr", "closure", "czm_replacement",
    "multiplicative_replacement", "aitchison_distance", "prevalence_filter",
    "aldex2", "AldexResult",
    "ancombc2", "AncomBCResult",
    "is_contaminant", "DecontamResult",
    "mixed_model_scan", "variance_partition",
    "simulate_study", "SimulatedStudy",
    "concordance_table", "agreement_summary", "evaluate_against_truth",
]
