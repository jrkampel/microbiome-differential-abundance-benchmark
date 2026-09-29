# Contributing

## Running things locally

```bash
pip install -e ".[dev]"
pytest
make all
```

An editable install means `import mdab` works everywhere. The scripts
also insert `src/` on the path themselves, so they run without it.

## What a change needs

**A statistical change needs a test that would fail without it.** The
tests that matter here are the invariance tests in `tests/test_mdab.py`:
CLR invariance to library size, Aitchison scale invariance and
subcompositional dominance, and the label-permutation null check. A
method change that breaks one of those is wrong regardless of how the
results look.

**A change to a reported number means updating the README.** The results
section quotes specific counts from the demonstration run. If the
simulator or a method changes, rerun `make all` and update them, or
delete them. Stale numbers in a README are worse than no numbers.

**Assumptions go in `docs/ASSUMPTIONS.md`, not in a code comment.** If a
change introduces a new assumption or removes one, that file is the
place it gets recorded.

## What this repository deliberately does not do

It does not pick a winner between ALDEx2 and ANCOM-BC2, and pull requests
framed that way will be asked to reframe. The point is to characterise
where they differ and what each difference costs.

It does not claim to recover absolute abundance. Any change implying it
does needs external quantification in the design, not a cleverer
transform.
