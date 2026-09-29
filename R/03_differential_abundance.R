#!/usr/bin/env Rscript
# Reference-implementation run: decontam, ALDEx2 and ANCOM-BC2.
#
# This script is the authority for results on real data. The Python
# package in src/mdab reimplements the same statistics so that the logic
# is readable and so the pipeline can run without Bioconductor, and
# tests/test_against_reference.py compares the two on shared inputs.
# Where they differ, this script wins.

suppressPackageStartupMessages({
  library(phyloseq)
  library(decontam)
  library(ALDEx2)
  library(ANCOMBC)
  library(TreeSummarizedExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
proc_dir   <- if (length(args) >= 1) args[[1]] else "data/processed"
out_dir    <- if (length(args) >= 2) args[[2]] else "results/tables"
group_col  <- if (length(args) >= 3) args[[3]] else "arm"
random_col <- if (length(args) >= 4) args[[4]] else "donor"
alpha      <- 0.05
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts <- as.matrix(read.delim(file.path(proc_dir, "counts.tsv"), row.names = 1))
meta   <- read.delim(file.path(proc_dir, "metadata.tsv"), row.names = 1)
meta   <- meta[rownames(counts), , drop = FALSE]

ps <- phyloseq(
  otu_table(counts, taxa_are_rows = FALSE),
  sample_data(meta)
)

## ---------------------------------------------------------------------
## 1. decontam
## ---------------------------------------------------------------------
# Run before anything else. The combined method uses both the prevalence
# contrast against negative controls and the inverse relationship between
# contaminant frequency and total DNA. Threshold 0.1 is the package
# default and is intentionally permissive.
if ("is_control" %in% colnames(meta) && any(as.logical(meta$is_control))) {
  contam <- isContaminant(
    ps,
    method = if ("dna_ng_ul" %in% colnames(meta)) "combined" else "prevalence",
    neg = as.logical(meta$is_control),
    conc = if ("dna_ng_ul" %in% colnames(meta)) meta$dna_ng_ul else NULL,
    threshold = 0.1,
    batch = if ("batch" %in% colnames(meta)) meta$batch else NULL
  )
  write.table(contam, file.path(out_dir, "decontam_scores_R.tsv"),
              sep = "\t", quote = FALSE)
  message(sprintf("[decontam] flagged %d of %d features",
                  sum(contam$contaminant, na.rm = TRUE), nrow(contam)))
  ps <- prune_taxa(!contam$contaminant, ps)
} else {
  message("[decontam] no negative controls; skipping. Record this in ASSUMPTIONS.md")
}

# Drop the controls themselves and any feature that is now empty.
ps <- prune_samples(!as.logical(sample_data(ps)$is_control), ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# Prevalence filter. Report the threshold, it changes the CLR reference.
prev <- apply(otu_table(ps) > 0, 2, mean)
ps <- prune_taxa(prev >= 0.10, ps)
message(sprintf("[filter] %d samples, %d features retained",
                nsamples(ps), ntaxa(ps)))

counts_f <- as(otu_table(ps), "matrix")
meta_f <- data.frame(sample_data(ps))

## ---------------------------------------------------------------------
## 2. ALDEx2
## ---------------------------------------------------------------------
# ALDEx2 expects features as rows. mc.samples = 128 is the default;
# denom = "all" is the plain CLR, "iqlr" restricts the reference set and
# is the right choice when a large share of the community moves one way.
aldex_in <- t(counts_f)
conds <- as.character(meta_f[[group_col]])

set.seed(1)
ax <- aldex(
  aldex_in, conds,
  mc.samples = 128,
  test = "t",
  effect = TRUE,
  denom = "all",
  verbose = TRUE
)
write.table(ax, file.path(out_dir, "aldex2_results_R.tsv"), sep = "\t", quote = FALSE)
message(sprintf("[aldex2] %d features with we.eBH < %.2f",
                sum(ax$we.eBH < alpha), alpha))
message(sprintf("[aldex2] %d also with |effect| > 1",
                sum(ax$we.eBH < alpha & abs(ax$effect) > 1)))

## ---------------------------------------------------------------------
## 3. ANCOM-BC2
## ---------------------------------------------------------------------
# rand_formula puts a random intercept on donor, which is the correct
# handling of repeated measures. pseudo_sens = TRUE runs the pseudocount
# sensitivity analysis; taxa failing it should not be reported as
# findings no matter how small their q value.
tse <- TreeSummarizedExperiment(
  assays = list(counts = t(counts_f)),
  colData = meta_f
)

anc <- ancombc2(
  data = tse,
  assay_name = "counts",
  fix_formula = group_col,
  rand_formula = if (!is.na(random_col) && random_col %in% colnames(meta_f))
    paste0("(1 | ", random_col, ")") else NULL,
  p_adj_method = "BH",
  prv_cut = 0.10,
  lib_cut = 1000,
  group = group_col,
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = alpha,
  pseudo_sens = TRUE,
  global = FALSE,
  n_cl = 4
)

res <- anc$res
write.table(res, file.path(out_dir, "ancombc2_results_R.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
if (!is.null(anc$zero_ind)) {
  write.table(anc$zero_ind, file.path(out_dir, "ancombc2_structural_zeros_R.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

diff_cols <- grep("^diff_", colnames(res), value = TRUE)
ss_cols <- grep("^passed_ss", colnames(res), value = TRUE)
target <- setdiff(diff_cols, "diff_(Intercept)")[1]
message(sprintf("[ancombc2] %d features with diff = TRUE", sum(res[[target]])))
if (length(ss_cols)) {
  target_ss <- setdiff(ss_cols, "passed_ss_(Intercept)")[1]
  message(sprintf("[ancombc2] %d of those also pass pseudocount sensitivity",
                  sum(res[[target]] & res[[target_ss]])))
}

## ---------------------------------------------------------------------
## 4. Concordance
## ---------------------------------------------------------------------
shared <- intersect(rownames(ax), res$taxon)
comparison <- data.frame(
  taxon = shared,
  aldex_effect = ax[shared, "effect"],
  aldex_q = ax[shared, "we.eBH"],
  aldex_call = ax[shared, "we.eBH"] < alpha,
  ancom_lfc = res[match(shared, res$taxon), grep("^lfc_", colnames(res))[2]],
  ancom_call = res[match(shared, res$taxon), target]
)
comparison$agreement <- with(comparison, ifelse(
  aldex_call & ancom_call, "both",
  ifelse(aldex_call, "aldex2_only",
         ifelse(ancom_call, "ancombc2_only", "neither"))))
write.table(comparison, file.path(out_dir, "concordance_R.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(table(comparison$agreement))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("[done]")
