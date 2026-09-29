#!/usr/bin/env python3
"""Generate the benchmark dataset with known ground truth.

Run this when you cannot reach ENA, or when you want to check that the
downstream analysis recovers signal you planted. The output files have
exactly the same shape as the ones produced by the DADA2 route, so every
later script is agnostic about which one made them.

Outputs
-------
``data/processed/counts.tsv``   features as columns, samples as rows
``data/processed/metadata.tsv`` one row per sample
``data/processed/truth.tsv``    per-feature ground truth, simulation only
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mdab import simulate_study  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", type=int, default=30)
    ap.add_argument("--taxa", type=int, default=300)
    ap.add_argument("--da", type=int, default=30, help="truly changed taxa")
    ap.add_argument("--effect", type=float, default=1.5, help="log fold change")
    ap.add_argument("--biomass-ratio", type=float, default=0.6,
                    help="absolute load in the treated arm; below 1 mimics a drug "
                         "that reduces biomass, which no compositional method can see")
    ap.add_argument("--symmetric", action="store_true",
                    help="split effect directions instead of the harder one-sided case")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--outdir", type=Path, default=Path("data/processed"))
    args = ap.parse_args()

    study = simulate_study(
        n_donors=args.donors,
        n_taxa=args.taxa,
        n_da=args.da,
        effect_lfc=args.effect,
        biomass_ratio=args.biomass_ratio,
        asymmetric=not args.symmetric,
        seed=args.seed,
    )

    args.outdir.mkdir(parents=True, exist_ok=True)
    study.counts.to_csv(args.outdir / "counts.tsv", sep="\t")
    study.metadata.to_csv(args.outdir / "metadata.tsv", sep="\t")
    study.truth.to_csv(args.outdir / "truth.tsv", sep="\t")

    n_samp = int((study.metadata["sample_type"] == "sample").sum())
    n_ctrl = int((study.metadata["sample_type"] == "negative_control").sum())
    print(f"samples {n_samp}, negative controls {n_ctrl}, features {study.counts.shape[1]}")
    print(f"library size median {int(study.metadata['library_size'].median())}, "
          f"range {int(study.metadata['library_size'].min())}"
          f"-{int(study.metadata['library_size'].max())}")
    print(f"zero fraction {float((study.counts == 0).to_numpy().mean()):.2%}")
    print(f"truly differential {int(study.truth['is_da'].sum())}, "
          f"contaminants {int(study.truth['is_contaminant'].sum())}")
    print(f"written to {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
