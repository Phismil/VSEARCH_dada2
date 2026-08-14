library(dada2); library(Biostrings)

path <- getwd()

fnFs         <- sort(list.files(path, pattern = "R1_001.fastq", full.names = TRUE))

sample.names <- sapply(base::strsplit(basename(fnFs), "_"), function(x) paste(x[2:2], collapse = ""))

getN         <- function(x) sum(getUniques(x))

# Quality profile

#png("quality_profile_F.png", width = 1200, height = 800)
#3print(plotQualityProfile(fnFs[1:min(4, length(fnFs))]))
#dev.off()

# Filter — adjust trimLeft, maxEE, minLen as needed

filtFs <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq.gz"))

names(filtFs) <- sample.names

dir.create(file.path(path, "filtered"), showWarnings = FALSE)

out <- filterAndTrim(fnFs, filtFs, trimLeft = 0, maxN = 0, maxEE = 5,
                     minLen = 100, rm.phix = TRUE, compress = TRUE, multithread = TRUE)

keep         <- file.exists(filtFs)

filtFs       <- filtFs[keep]

sample.names <- sample.names[keep]

out          <- out[keep, , drop = FALSE]

track <- data.frame(input = out[,1], filtered = out[,2], row.names = sample.names)

cat(sprintf("Filtered: %d / %d (%.1f%%)\n",
    sum(track$filtered), sum(track$input),
    100 * sum(track$filtered) / sum(track$input)))

# Error rates

errF <- learnErrors(filtFs, multithread = TRUE)

png("error_rates_F.png", width = 1000, height = 700, res = 120)

print(plotErrors(errF, nominalQ = TRUE))

dev.off()

# Denoise

dadaFs <- dada(filtFs, err = errF, multithread = TRUE, pool = FALSE)

track$denoised <- if (is(dadaFs, "dada")) getN(dadaFs) else sapply(dadaFs, getN)

cat("Denoised:", sum(track$denoised), "\n")

# Sequence table

seqtab       <- makeSequenceTable(dadaFs)

track$tabled <- rowSums(seqtab)

cat("ASV table:", nrow(seqtab), "samples,", ncol(seqtab), "ASVs\n")

print(table(nchar(getSequences(seqtab))))

# Chimera removal

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = TRUE, verbose = TRUE)

track$nonchim <- rowSums(seqtab.nochim)

cat(sprintf("\nFinal: %d ASVs, %d reads\n", ncol(seqtab.nochim), sum(track$nonchim)))

# Summary

print(track)

print(summary(track$nonchim / track$tabled))          # chimera retention

cat("Overall retention:", sum(track$nonchim) / sum(track$input), "\n")

# Export

asv_seqs <- colnames(seqtab.nochim)

asv_ids  <- paste0("ASV_F_", seq_along(asv_seqs))

colnames(seqtab.nochim) <- asv_ids

write.table(t(seqtab.nochim), file.path(path, "ASV_table.tsv"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

writeLines(c(rbind(paste0(">",asv_ids), asv_seqs)), file.path(path, "ASVs.fasta"))

write.csv(track, file.path(path, "dada2_tracking_stats.csv"), row.names = TRUE)

cat("Exported: ASV_table.tsv, ASVs.fasta, dada2_tracking_stats.csv\n")
