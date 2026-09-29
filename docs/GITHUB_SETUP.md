# Pushing this to GitHub

The repository is already initialised with one commit on `main`. Nothing
below rewrites history; it only adds a remote and pushes.

## 1. Create the empty remote

On github.com, create a new repository. Do **not** tick "add a README",
"add .gitignore" or "choose a licence": all three already exist here and
an initialised remote will cause a rejected push on the first attempt.

Suggested name: `microbiome-da-benchmark`

Suggested description, which is short enough for the sidebar:

> Compositional differential abundance for 16S data: ALDEx2 and ANCOM-BC2
> compared against simulated ground truth, with donor random effects and
> decontam.

Suggested topics: `microbiome`, `16s-rrna`, `bioinformatics`,
`compositional-data`, `differential-abundance`, `aldex2`, `ancombc`,
`dada2`, `reproducible-research`

## 2. Push

```bash
cd microbiome-da-benchmark

# set your own identity on the existing commit
git config user.name  "Your Name"
git config user.email "your@email"
git commit --amend --reset-author --no-edit

git remote add origin git@github.com:USERNAME/microbiome-da-benchmark.git
git push -u origin main
```

Use the HTTPS URL instead if you have not set up SSH keys:

```bash
git remote add origin https://github.com/USERNAME/microbiome-da-benchmark.git
```

HTTPS will ask for a personal access token, not your password. Generate
one under Settings, Developer settings, Personal access tokens, with
`repo` scope.

## 3. Check CI

`.github/workflows/tests.yml` runs on the first push. It does two jobs:

- **test** installs the package on Python 3.10, 3.11 and 3.12, runs
  `ruff check`, then `pytest`.
- **pipeline** runs the simulated route end to end and fails if any stage
  errors or stops producing its expected outputs, then uploads
  `results/` as a build artefact.

Both pass locally as of the initial commit. If the pipeline job times out
on the free runner, drop `scripts/06_simulation_study.py` from it; it is
already excluded for that reason.

Add the badge to the top of the README once the first run is green:

```markdown
[![tests](https://github.com/USERNAME/microbiome-da-benchmark/actions/workflows/tests.yml/badge.svg)](https://github.com/USERNAME/microbiome-da-benchmark/actions/workflows/tests.yml)
```

## 4. Before anyone else reads it

Two things are currently overstated by omission rather than by claim, and
both are cheap to fix:

1. **`scripts/00_fetch_public_data.py` has never been run against a live
   accession.** It is written and documented but untested. Either run one
   real accession through it and commit the resulting sample sheet, or add
   a line to the README saying the simulation route is validated and the
   real-data route is implemented but untested. The second is honest and
   takes a minute; the first is better.

2. **Every number in the README results section is from simulated data.**
   The section heading says so. Keep it explicit anywhere the numbers get
   quoted away from that heading, including a CV or cover letter.

## 5. Things worth doing next, roughly in order of value

- Run one public dataset end to end through the R route and add its
  results as a second section. This is the single largest credibility
  gain available.
- Add a `--denom iqlr` comparison to the results, since the repository
  already argues that the moving CLR reference matters and currently
  demonstrates the problem without demonstrating that particular fix.
- Cross-check the Python implementations against the R packages on the
  same input and commit the comparison. `tests/` currently checks
  internal consistency and statistical behaviour, not agreement with
  ALDEx2 and ANCOMBC themselves.
- Pin package versions in `environment.yml` once you have a working R
  environment, so the Bioconductor route is reproducible.
