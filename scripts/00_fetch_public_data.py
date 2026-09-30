#!/usr/bin/env python3
"""Fetch a public 16S study with a drug or perturbation arm from ENA.

This script talks to the ENA Portal API, which is the most reliable route
to raw reads. SRA mirrors the same data but its toolkit is heavier and
Qiita requires a login for some studies. All three are listed in
``docs/DATASETS.md`` with the tradeoffs.

Why ENA rather than SRA
-----------------------
The ENA filereport endpoint returns a plain TSV with direct FTP paths to
the FASTQ files, so a study can be downloaded with curl and no SRA
toolkit. The same endpoint returns library layout and instrument, which
you need in order to set the DADA2 truncation parameters correctly.

Candidate studies with a perturbation arm
-----------------------------------------
See ``docs/DATASETS.md``. The default below is a metformin study because
it has a clean drug arm, repeated sampling per donor, and published
results to compare against. Any accession with a treatment column in its
sample metadata will work.

Important caveat about metadata
-------------------------------
ENA sample attributes are entered by submitters and are frequently
incomplete, inconsistently spelled, or missing the variable you actually
need. Expect to curate the sample sheet by hand and to cross-reference
the paper's supplementary tables. Do not skip this and do not assume a
column called "treatment" means what you think it means. Record every
manual correction in ``data/metadata/curation_log.md``.

Usage
-----
    python scripts/00_fetch_public_data.py --accession PRJNA000000 \\
        --outdir data/raw --max-samples 40

Network note
------------
This script needs outbound access to ``ftp.ebi.ac.uk`` and
``www.ebi.ac.uk``. It will not run inside a sandbox that restricts
egress, which is why the rest of the pipeline can also be driven from the
simulator in ``scripts/01_simulate_dataset.py``.
"""

from __future__ import annotations

import argparse
import csv
import io
import sys
import urllib.parse
import urllib.request
from pathlib import Path

PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"

FIELDS = [
    "run_accession", "sample_accession", "experiment_accession",
    "library_strategy", "library_layout", "instrument_model",
    "read_count", "fastq_ftp", "fastq_md5", "sample_title",
    "sample_alias", "description",
]


def fetch_filereport(accession: str) -> list[dict]:
    query = urllib.parse.urlencode(
        {
            "accession": accession,
            "result": "read_run",
            "fields": ",".join(FIELDS),
            "format": "tsv",
            "limit": "0",
        }
    )
    url = f"{PORTAL}?{query}"
    print(f"[ena] {url}", file=sys.stderr)
    with urllib.request.urlopen(url, timeout=120) as resp:
        text = resp.read().decode("utf-8")
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def fetch_sample_attributes(accession: str) -> list[dict]:
    """Sample level metadata, where the treatment arm usually hides."""
    query = urllib.parse.urlencode(
        {
            "accession": accession,
            "result": "sample",
            "fields": "sample_accession,sample_title,description,"
                      "collection_date,host_subject_id,sample_alias",
            "format": "tsv",
            "limit": "0",
        }
    )
    with urllib.request.urlopen(f"{PORTAL}?{query}", timeout=120) as resp:
        text = resp.read().decode("utf-8")
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def download(url: str, dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[skip] {dest.name}", file=sys.stderr)
        return
    if not url.startswith(("http://", "https://", "ftp://")):
        url = "https://" + url
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"[get ] {dest.name}", file=sys.stderr)
    with urllib.request.urlopen(url, timeout=600) as resp, open(tmp, "wb") as fh:
        while chunk := resp.read(1 << 20):
            fh.write(chunk)
    tmp.rename(dest)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accession", required=True, help="ENA or SRA project accession")
    ap.add_argument("--outdir", default="data/raw", type=Path)
    ap.add_argument("--max-samples", type=int, default=0,
                    help="0 downloads everything; set a small number first")
    ap.add_argument("--metadata-only", action="store_true",
                    help="write the sample sheet and stop, useful for triage")
    args = ap.parse_args()

    runs = fetch_filereport(args.accession)
    if not runs:
        print("no runs returned; check the accession", file=sys.stderr)
        return 1

    strategies = {r["library_strategy"] for r in runs}
    layouts = {r["library_layout"] for r in runs}
    instruments = {r["instrument_model"] for r in runs}
    print(f"[info] {len(runs)} runs, strategy={strategies}, "
          f"layout={layouts}, instrument={instruments}", file=sys.stderr)
    if "AMPLICON" not in strategies:
        print("[warn] library_strategy is not AMPLICON; this may be shotgun data",
              file=sys.stderr)

    meta_dir = Path("data/metadata")
    meta_dir.mkdir(parents=True, exist_ok=True)
    sheet = meta_dir / f"{args.accession}_runs.tsv"
    with open(sheet, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        for row in runs:
            writer.writerow({k: row.get(k, "") for k in FIELDS})
    print(f"[out ] {sheet}", file=sys.stderr)

    try:
        samples = fetch_sample_attributes(args.accession)
        if samples:
            s_sheet = meta_dir / f"{args.accession}_samples.tsv"
            with open(s_sheet, "w", newline="") as fh:
                writer = csv.DictWriter(fh, fieldnames=list(samples[0].keys()),
                                        delimiter="\t")
                writer.writeheader()
                writer.writerows(samples)
            print(f"[out ] {s_sheet}", file=sys.stderr)
    except Exception as exc:
        print(f"[warn] sample attributes unavailable: {exc}", file=sys.stderr)

    if args.metadata_only:
        return 0

    selected = runs[: args.max_samples] if args.max_samples else runs
    for row in selected:
        urls = [u for u in row["fastq_ftp"].split(";") if u]
        for url in urls:
            download(url, args.outdir / Path(url).name)

    print(f"[done] {len(selected)} runs into {args.outdir}", file=sys.stderr)
    print("Next: verify md5 sums, then run R/01_dada2.R", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
