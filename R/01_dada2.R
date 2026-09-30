#!/usr/bin/env Rscript
# DADA2 processing of paired-end 16S reads into an ASV table.
#
# Why DADA2 rather than OTU clustering
# ------------------------------------
# Clustering at 97 percent identity was a workaround for sequencing error.
# It merges genuinely distinct organisms and produces units that are not
# comparable between studies. DADA2 instead learns the run-specific error
# model and resolves exact amplicon sequence variants, which are single
# nucleotide resolution and are comparable across studies because the
# sequence itself is the label. Use ASVs.
#
# The parameters that actually matter
# -----------------------------------
# truncLen is the one people get wrong. Reads must be truncated where
# quality collapses, but the forward and reverse truncation points must
# still leave enough overlap to merge. For the V4 region amplified by
# 515F/806R the amplicon is about 253 bp, so truncLen values must sum to
# at least 253 plus roughly 20 bp of overlap plus the primer lengths.
# Truncating too aggressively silently destroys the merge step and you
# end up with a tiny ASV table and no error message that says why.
#
# Always inspect plotQualityProfile before setting these. Do not copy the
# values below into a different study.

suppressPackageStartupMessages({
  library(dada2)
  library(Biostrings)
})

args <- commandArgs(trailingOnly = TRUE)
raw_dir    <- if (length(args) >= 1) args[[1]] else "data/raw"
out_dir    <- if (length(args) >= 2) args[[2]] else "data/processed"
trunc_fwd  <- if (length(args) >= 3) as.integer(args[[3]]) else 240L
trunc_rev  <- if (length(args) >= 4) as.integer(args[[4]]) else 200L
n_threads  <- if (length(args) >= 5) as.integer(args[[5]]) else 4L

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "filtered"), recursive = TRUE, showWarnings = FALSE)

fnFs <- sort(list.files(raw_dir, pattern = "_1\\.fastq\\.gz$", full.names = TRUE))
fnRs <- sort(list.files(raw_dir, pattern = "_2\\.fastq\\.gz$", full.names = TRUE))
stopifnot(length(fnFs) > 0, length(fnFs) == length(fnRs))
sample_names <- sub("_1\\.fastq\\.gz$", "", basename(fnFs))
message(sprintf("[dada2] %d paired samples", length(sample_names)))

# Quality profiles. Look at these before trusting truncLen.
pdf(file.path(out_dir, "quality_profiles.pdf"), width = 9, height = 5)
print(plotQualityProfile(fnFs[seq_len(min(4, length(fnFs)))]))
print(plotQualityProfile(fnRs[seq_len(min(4, length(fnRs)))]))
dev.off()

filtFs <- file.path(out_dir, "filtered", paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(out_dir, "filtered", paste0(sample_names, "_R_filt.fastq.gz"))

# maxEE is an expected-errors filter and is strictly better than a mean
# quality cutoff, because it accumulates the per-base error probabilities
# rather than averaging a log scale.
out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  truncLen = c(trunc_fwd, trunc_rev),
  maxN = 0, maxEE = c(2, 2), truncQ = 2,
  rm.phix = TRUE, compress = TRUE, multithread = n_threads
)
print(head(out))

keep <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[keep]; filtRs <- filtRs[keep]
sample_names <- sample_names[keep]

# The error model is learned per run. If samples come from several
# sequencing runs, split them and learn errors separately, then merge the
# sequence tables. Pooling runs into one error model is a common and
# consequential mistake.
errF <- learnErrors(filtFs, multithread = n_threads, randomize = TRUE)
errR <- learnErrors(filtRs, multithread = n_threads, randomize = TRUE)
pdf(file.path(out_dir, "error_model.pdf"), width = 8, height = 6)
print(plotErrors(errF, nominalQ = TRUE))
print(plotErrors(errR, nominalQ = TRUE))
dev.off()

# pool = "pseudo" shares information across samples so that a variant
# seen once in each of many samples can be resolved. Full pooling is more
# sensitive but scales badly.
ddF <- dada(filtFs, err = errF, multithread = n_threads, pool = "pseudo")
ddR <- dada(filtRs, err = errR, multithread = n_threads, pool = "pseudo")

merged <- mergePairs(ddF, filtFs, ddR, filtRs, minOverlap = 12, verbose = TRUE)
seqtab <- makeSequenceTable(merged)
message(sprintf("[dada2] %d ASVs before chimera removal", ncol(seqtab)))
print(table(nchar(getSequences(seqtab))))

# Chimeras are PCR artefacts formed from two parent templates. They are
# abundant and they look like novel taxa, so removal is not optional.
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = n_threads, verbose = TRUE)
message(sprintf("[dada2] %d ASVs after chimera removal, %.1f%% of reads retained",
                ncol(seqtab_nochim), 100 * sum(seqtab_nochim) / sum(seqtab)))

# Read tracking. A large drop at any stage is a parameter problem, not
# biology. A drop at merging almost always means truncLen was too short.
getN <- function(x) sum(getUniques(x))
track <- cbind(out[keep, ], sapply(ddF, getN), sapply(ddR, getN),
               sapply(merged, getN), rowSums(seqtab_nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR",
                     "merged", "nonchim")
rownames(track) <- sample_names
write.table(track, file.path(out_dir, "read_tracking.tsv"),
            sep = "\t", quote = FALSE)
print(track)

# Store sequences separately and give ASVs short stable names, otherwise
# every downstream table has 253-character column names.
asv_seqs <- colnames(seqtab_nochim)
asv_ids <- sprintf("ASV%04d", seq_along(asv_seqs))
writeXStringSet(setNames(DNAStringSet(asv_seqs), asv_ids),
                file.path(out_dir, "asv_sequences.fasta"))
colnames(seqtab_nochim) <- asv_ids
rownames(seqtab_nochim) <- sample_names

write.table(seqtab_nochim, file.path(out_dir, "counts.tsv"),
            sep = "\t", quote = FALSE)
saveRDS(seqtab_nochim, file.path(out_dir, "seqtab_nochim.rds"))
message("[dada2] done")
