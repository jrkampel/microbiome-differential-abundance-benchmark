# Choosing a public 16S dataset with a perturbation arm

## What to look for

A dataset is usable here if it has all four of these. Check before
downloading anything, because the download is the slow part and the
metadata is where studies fail.

1. **A real perturbation arm.** A drug, an antibiotic course, a dietary
   intervention, a faecal transplant, a probiotic. Observational
   case-control studies work too but confound the treatment with disease
   status.
2. **Repeated sampling per subject.** Before and after at minimum. This
   is what makes the random effect meaningful and it substantially
   increases power, because between-person variation is removed from the
   comparison.
3. **Machine-readable sample metadata** that actually names the arm and
   the subject. This is the usual failure point. ENA sample attributes
   are submitter-entered and frequently omit the treatment variable, in
   which case you need the paper's supplementary tables.
4. **A stated primer pair and region.** Without it you cannot set DADA2
   truncation lengths correctly or trim the reference database.

Desirable but rare: negative controls deposited alongside the samples,
and any absolute quantification such as qPCR of total 16S copies or flow
cytometry counts. If a study has either, prefer it, because it lets you
say something about absolute abundance that no compositional method can.

## Where to look

**ENA Portal API.** The best route. `scripts/00_fetch_public_data.py`
uses `result=read_run` for FASTQ paths and `result=sample` for
attributes. Returns plain TSV with direct FTP links, no toolkit needed.

**SRA.** Mirrors ENA content. Needs `fasterq-dump` from the SRA toolkit
and is slower. Use when a study is not mirrored to ENA. SRA Run Selector
often has better-curated metadata columns than ENA sample attributes.

**Qiita.** Studies come pre-processed with standardised sample metadata
under the Qiimp templates, which is a large advantage. Some studies need
a login, and the deposited processing may not match what you want, so
download raw and reprocess.

**MGnify.** Useful for finding studies by biome, then follow the
accession back to ENA for raw reads.

**curatedMetagenomicData.** Shotgun rather than amplicon, but the
metadata curation is the standard to aim for, and worth reading as an
example of how a sample sheet should look.

## Triage workflow

```bash
# 1. Metadata only, fast
python scripts/00_fetch_public_data.py --accession PRJNAxxxxxx --metadata-only

# 2. Inspect. Does a treatment column exist? A subject column?
#    How many runs per subject? Single or paired end? Which instrument?
column -t -s $'\t' data/metadata/PRJNAxxxxxx_runs.tsv | less -S

# 3. Curate by hand into data/metadata/sample_sheet.tsv with at least:
#    sample_id, donor, arm, batch, is_control
#    Log every manual correction in data/metadata/curation_log.md

# 4. Small pilot download before committing to the whole study
python scripts/00_fetch_public_data.py --accession PRJNAxxxxxx --max-samples 8
```

## Warnings

**Check `library_strategy`.** A study can contain both amplicon and
shotgun runs. Filter to `AMPLICON` before downloading.

**Check `library_layout`.** The DADA2 script assumes paired-end with
`_1` and `_2` suffixes. Single-end data needs the forward-only DADA2
route, which is a small edit but must be a deliberate one.

**Check for multiple sequencing runs.** If `instrument_model` or run
dates vary, learn error models separately per run and merge the sequence
tables afterwards. Treat run as a batch variable.

**Check whether reads are already demultiplexed and primer-trimmed.**
Some depositors upload raw multiplexed reads with barcodes still
attached. Running DADA2 on those produces nonsense.

**Do not trust a column named `treatment`.** Read the paper and confirm
what the levels mean, including which one is the control and whether a
"baseline" sample was taken before or after any washout.
