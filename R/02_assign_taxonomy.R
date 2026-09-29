#!/usr/bin/env Rscript
# Taxonomic assignment of ASVs against SILVA or GTDB.
#
# What the classifier actually does
# ---------------------------------
# assignTaxonomy implements the RDP naive Bayes classifier. The query
# sequence is decomposed into 8-mers, the likelihood of that 8-mer profile
# is computed under each reference genus, and the assignment is
# bootstrapped by resampling 8-mers. The bootstrap value is a measure of
# how stable the assignment is to the choice of k-mers, not a posterior
# probability that the organism is that genus. A confidence of 80 is the
# convention and is already permissive at genus level.
#
# The limits you must state in the README
# ---------------------------------------
# 1. The 16S gene does not resolve most species. Species level calls from
#    a 250 bp region are unreliable and frequently wrong. addSpecies only
#    reports exact matches, which is the honest version, and even then a
#    unique exact match does not imply a unique organism.
# 2. Databases disagree. SILVA is curated on rRNA phylogeny, GTDB on whole
#    genome relatedness with substantial renaming and splitting of genera.
#    A "Clostridium" in one is not a "Clostridium" in the other. Never mix
#    assignments from two databases in one analysis, and name the exact
#    release in every figure caption.
# 3. Reference databases are biased towards cultured and clinically
#    relevant organisms. Unassigned does not mean absent, it means
#    unrepresented, and the unassigned fraction is not random.
# 4. Classifiers do better when the reference is trimmed to the amplified
#    region. If you have the primers, extract the region first.
# 5. Collapsing ASVs to genus discards real resolution and can merge taxa
#    with opposite responses to the treatment. Run differential abundance
#    at ASV level and aggregate for presentation, not the other way round.

suppressPackageStartupMessages({
  library(dada2)
  library(Biostrings)
})

args <- commandArgs(trailingOnly = TRUE)
proc_dir  <- if (length(args) >= 1) args[[1]] else "data/processed"
ref_fasta <- if (length(args) >= 2) args[[2]] else "data/reference/silva_nr99_v138.2_toSpecies_trainset.fa.gz"
species_fasta <- if (length(args) >= 3) args[[3]] else NA_character_
n_threads <- if (length(args) >= 4) as.integer(args[[4]]) else 4L
min_boot  <- 80

if (!file.exists(ref_fasta)) {
  stop(paste0(
    "reference not found: ", ref_fasta, "\n",
    "Download a DADA2-formatted training set before running this.\n",
    "SILVA releases are on Zenodo, GTDB training sets are distributed\n",
    "by the DADA2 maintainers. Record the exact release and its md5 in\n",
    "docs/ASSUMPTIONS.md, because taxonomy is not reproducible without it."
  ))
}

seqs <- readDNAStringSet(file.path(proc_dir, "asv_sequences.fasta"))
asv_ids <- names(seqs)
seq_chr <- as.character(seqs)

message(sprintf("[tax] classifying %d ASVs against %s", length(seq_chr), basename(ref_fasta)))
taxa <- assignTaxonomy(seq_chr, ref_fasta, multithread = n_threads,
                       minBoot = min_boot, outputBootstraps = TRUE)

tax_table <- as.data.frame(taxa$tax, stringsAsFactors = FALSE)
boot_table <- as.data.frame(taxa$boot)
rownames(tax_table) <- asv_ids
rownames(boot_table) <- asv_ids

if (!is.na(species_fasta) && file.exists(species_fasta)) {
  message("[tax] exact-match species assignment")
  sp <- addSpecies(as.matrix(tax_table), species_fasta, allowMultiple = TRUE)
  tax_table <- as.data.frame(sp, stringsAsFactors = FALSE)
  rownames(tax_table) <- asv_ids
} else {
  message("[tax] skipping species assignment; 16S rarely supports it anyway")
}

# Report the unassigned fraction at every rank. This is the number to put
# in the README, because it bounds how much of the analysis can be
# interpreted biologically.
for (rank in colnames(tax_table)) {
  frac <- mean(is.na(tax_table[[rank]]))
  message(sprintf("[tax] unassigned at %-8s %.1f%%", rank, 100 * frac))
}

write.table(tax_table, file.path(proc_dir, "taxonomy.tsv"), sep = "\t", quote = FALSE)
write.table(boot_table, file.path(proc_dir, "taxonomy_bootstraps.tsv"),
            sep = "\t", quote = FALSE)

writeLines(
  c(
    "# Taxonomy provenance",
    paste("reference:", normalizePath(ref_fasta)),
    paste("md5:", tools::md5sum(ref_fasta)),
    paste("minBoot:", min_boot),
    paste("dada2 version:", as.character(packageVersion("dada2"))),
    paste("date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  ),
  file.path(proc_dir, "taxonomy_provenance.txt")
)
message("[tax] done")
