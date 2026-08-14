# Complete DADA2 Workflow Script with Per-Step Stats Tracking
# Built for Illumina paired-end reads (e.g., 16S amplicon data)
# Run in R or RStudio; adjust paths and parameters as needed.

# ----- Section 0: Load Packages and Set Paths ----- #

# On CHPC zlib is needed
#
# conda install anaconda::zlib
# conda install bioconda::bioconductor-rsamtools
# conda install liblzma
# conda install xz

# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) == 0) {
# stop("Please provide the path to the FASTQ directory as a command-line argument.")
# }
# path <- args[1]

library("dada2")

# Update this path to your FASTQ directory

path <- getwd()

# List forward and reverse FASTQ files

fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE)) # Forward reads
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE)) # Reverse reads

# Extract sample names from tokens 2 and 3 of underscore-delimited filename
# e.g., run01_sampleA_rep1_R1_001.fastq -> sampleA_rep1

sample.names <- sapply(strsplit(basename(fnFs), "_"), function(x) paste(x[2:3], collapse = "_"))

# Helper function for stats tracking (sums unique sequences)

getN <- function(x) sum(getUniques(x))

# Initialize tracking dataframe for stats per step

track <- data.frame(row.names = sample.names)

# ----- Section 1: Inspect Read Quality ----- #
# Visualize quality profiles to inform trimming (no numerical stats here)
# plotQualityProfile(fnFs[1:2])  # First 2 forward samples (adjust for more)
# plotQualityProfile(fnRs[1:2])  # First 2 reverse samples
# Tip: Look for quality drops; adjust truncLen/trimLeft in next section accordingly.

# ----- Section 2: Filter and Trim Reads ----- #
# Define output paths for filtered files

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
dir.create(file.path(path, "filtered"), showWarnings = FALSE) # Create dir if needed
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter and trim

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
  # truncLen = c(240, 160),  # Adjust: e.g., c(240, 160) for 2x250 bp reads
  # trimLeft = c(0, 10),     # Trim primers/adapters if present
  maxN = 0, maxEE = c(5, 5), rm.phix = TRUE,
  compress = TRUE, multithread = TRUE
)

# Update tracking: input and filtered reads

track$input    <- out[, "reads.in"]  # Raw read pairs
track$filtered <- out[, "reads.out"] # Filtered read pairs

# Display stats for this step

print("Stats after Filtering:")
print(track[, c("input", "filtered")])
print(summary(track$filtered / track$input)) # Retention rate

# Guard: remove samples with zero reads after filtering to prevent downstream errors

keep <- out[, "reads.out"] > 0

if (any(!keep)) {
  warning(paste(
    "Removing samples with 0 reads after filtering:",
    paste(sample.names[!keep], collapse = ", ")
  ))
}

filtFs       <- filtFs[keep]
filtRs       <- filtRs[keep]
sample.names <- sample.names[keep]
track        <- track[keep, , drop = FALSE]

# ----- Section 3: Learn Error Rates ----- #
# Learn errors from filtered reads

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Plot to verify (optional; no per-sample numerical stats here)

cat("Generating error rate plots...\n")

tryCatch(
  {
    png("errF_plot.png", width = 800, height = 600)
    print(plotErrors(errF, nominalQ = TRUE))
    dev.off()
    cat("Forward error rate plot saved: errF_plot.png\n")
  },
  error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat("Error generating forward error plot:", e$message, "\n")
  }
)

tryCatch(
  {
    png("errR_plot.png", width = 800, height = 600)
    print(plotErrors(errR, nominalQ = TRUE))
    dev.off()
    cat("Reverse error rate plot saved: errR_plot.png\n")
  },
  error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat("Error generating reverse error plot:", e$message, "\n")
  }
)

# Tip: Lines should fit points well; if not, check data quality.

# ----- Section 4: Denoise Samples (Sample Inference) ----- #
# Denoise forward and reverse reads

dadaFs <- dada(filtFs, err = errF, multithread = TRUE, pool = FALSE) # pool=TRUE for more sensitivity, but slower
dadaRs <- dada(filtRs, err = errR, multithread = TRUE, pool = FALSE)

# Update tracking: denoised unique sequences

track$denoisedF <- sapply(dadaFs, getN) # Forward
track$denoisedR <- sapply(dadaRs, getN) # Reverse

# Display updated stats

print("Stats after Denoising:")
print(track[, c("input", "filtered", "denoisedF", "denoisedR")])
print(summary(track$denoisedF / track$filtered)) # Retention rate

# ----- Section 5: Merge Paired Reads ----- #
# Merge forward and reverse reads

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs,
  minOverlap = 50, maxMismatch = 5, # Adjust for your amplicon overlap
  verbose = TRUE
)

# Update tracking: merged read pairs

track$merged <- sapply(mergers, getN)

# Display updated stats

print("Stats after Merging:")
print(track[, c("input", "filtered", "denoisedF", "denoisedR", "merged")])
print(summary(track$merged / track$denoisedF)) # Retention rate

# ----- Section 6: Construct Sequence Table ----- #
# Create ASV table

seqtab <- makeSequenceTable(mergers)
print(dim(seqtab)) # Dimensions: samples x ASVs

# Optional: Filter by sequence length (adjust range to match your primer set)
# V4 (515F/806R)   ~253 bp
# V3-V4 (341F/806R) ~400-450 bp
# seqtab <- seqtab[, nchar(colnames(seqtab)) %in% 400:450]  # Uncomment and adjust range

print(table(nchar(getSequences(seqtab)))) # Length distribution (stats)

# ----- Section 7: Remove Chimeras ----- #
# Remove chimeric sequences

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)

# Update tracking: non-chimeric reads

track$nonchim <- rowSums(seqtab.nochim)

# Display final pipeline stats

print("Final Stats after Chimera Removal:")
print(track)
print(summary(track$nonchim / track$merged)) # Chimera retention rate
print(sum(track$nonchim) / sum(track$input))  # Overall retention

# ----- Section 8: Export ASV Table and FASTA ----- #
# Assign ASV IDs (e.g., ASV_Merged1, ASV_Merged2, ...)

asv_seqs <- colnames(seqtab.nochim)
asv_ids  <- paste0("ASV_Merged", seq(1, length(asv_seqs)))

# Set ASV IDs as column names in the table

colnames(seqtab.nochim) <- asv_ids

# Export ASV table as TSV (ASVs as rows, samples as columns)

asv_table_path <- file.path(path, "ASV_table.tsv")
write.table(t(seqtab.nochim), file = asv_table_path, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
print(paste("ASV table exported to:", asv_table_path))

# Export ASV sequences as FASTA (paste0 avoids invalid space after >)

asv_fasta_path <- file.path(path, "ASVs.fasta")
fasta_lines    <- c(rbind(paste0(">", asv_ids), asv_seqs))
writeLines(fasta_lines, asv_fasta_path)
print(paste("ASV FASTA exported to:", asv_fasta_path))

# Export tracking stats as CSV for easy viewing

track_path <- file.path(path, "dada2_tracking_stats.csv")
write.csv(track, file = track_path, row.names = TRUE)
print(paste("Tracking stats exported to:", track_path))