











# VSEARCH_dada2

```bash
# VSEARCH
## COI

#!/bin/bash

# ==============================================================================
# VSEARCH PIPELINE FOR AMPLICON SEQUENCING
# ==============================================================================
# Description: Automates VSEARCH processing for Single-End and Paired-End reads.
# Dependencies: VSEARCH, Cutadapt, Trimmomatic, Java, Perl, BBMap
# ==============================================================================

# --- CONFIGURATION VARIABLES (EDIT THESE) ---
# #!/bin/bash

# Delete primer-trimmed reads

#
#
rm -rf merged unmergedf unmergedr otu

# Primers
PRIMER_FWD="GGWACWGGWTGAACWGTWTAYCCYCC"
PRIMER_REV="TANACYTCNGGRTGNCCRAARAAYCA"
RC_FWD="GGRGGRTAWACWGTTCAWCCWGTWCC"   # reverse complement of PRIMER_FWD
RC_REV="TGRTTYTTYGGNCAYCCNGARGTNTA"   # reverse complement of PRIMER_REV
LINKED_R1="^${PRIMER_FWD}...${RC_REV}"   # R1: FWD at 5' ... RC of REV at 3'
LINKED_R2="^${PRIMER_REV}...${RC_FWD}"   # R2: REV at 5' ... RC of FWD at 3'


#PRIMER_REV_NO_ANCHOR="TANACYTCNGGRTGNCCRAARAAYCA"



MARKER_LENGTH="317"
LABEL="COI"
MAX_EE="5"
MIN_LENGTH_MERGED="300"
MIN_OVERLAP="50"
MAX_DIF="5"
MAX_EE_FORWARD="5"
MAX_LENGTH_FORWARD="300"
THREADS="24"
MIN_SIZE="1"

# --- PATHS TO EXTERNAL SCRIPTS/DBS (EDIT THESE) ---
# Ensure these point to the correct locations on your cluster

PATH_MAP_SCRIPT="/mnt/lustre/users/aemami-khoyi/map.pl"

#PATH_TRIMMOMATIC="/apps/chpc/bio/trimmomatic/0.36/bin/trimmomatic.jar"
#PATH_BBMAP="/apps/chpc/bio/bbmap-38.95/readlength.sh"
#PATH_BLAST_DB=""
#BLAST_QUERY_FILE="merged.fasta" # Placeholder
# --- MODULE LOADING (CLUSTER SPECIFIC) ---

module purge
module add chpc/BIOMODULES
module add cutadapt/3.4
module add vsearch/2.18.0
# Add Java if needed for Trimmomatic explicitly, or rely on system java

# ==============================================================================
# FUNCTION: PAIRED END WORKFLOW
# ==============================================================================





for R1 in *_R1_001.fastq; do

    [ -e "$R1" ] || continue

    R2=${R1//R1_001.fastq/R2_001.fastq}
    #cutadapt -j 0 -g "$PRIMER_FWD" -G "$PRIMER_REV" --discard-untrimmed -o "ptrim_${R1}" -p "ptrim_${R2}" "$R1" "$R2"
    #cutadapt -j 0 -g "${LINKED_R1}" -G "${LINKED_R2}" --discard-untrimmed --pair-filter=any --minimum-length 200 -o "ptrim_${R1}" -p "ptrim_${R2}" "$R1" "$R2"
    cutadapt -j 0 -g "${PRIMER_FWD}" -a "${RC_REV}" -G "^${PRIMER_REV}" -A "${RC_FWD}" --error-rate 0.15 --discard-untrimmed --pair-filter=any -o "ptrim_${R1}" -p "ptrim_${R2}" "$R1" "$R2"
done


# 4. Merging Pairs

mkdir -p merged unmergedf unmergedr

    for f in ptrim_*_R1_001.fastq; do

    [ -e "$f" ] || continue

    r=$(sed -e "s/_R1_/_R2_/" <<< "$f")

    s=$(cut -d_ -f2,3 <<< "$f")

    echo "Merging Sample: $s"

    vsearch --fastq_mergepairs "$f" --reverse "$r" --fastq_minovlen "$MIN_OVERLAP" --fastq_maxdiffs "$MAX_DIF" --fastqout "${s}.merged.fastq" --fastq_eeout --fastqout_notmerged_fwd "${s}.notmergedforward.fastq" --fastqout_notmerged_rev "${s}.notmergedreverse.fastq"

    rm -rf "${s}.notmergedreverse.fastq"

    # Filtering

    vsearch --fastq_filter "${s}.merged.fastq" --fastq_maxee "$MAX_EE" --fastq_minlen "$MIN_LENGTH_MERGED" --fastq_maxlen "$MARKER_LENGTH" --fastq_maxns 0 --fastaout "${s}.filtered.fasta" --fasta_width 0


    vsearch --fastq_filter "${s}.notmergedforward.fastq" --fastq_maxee "$MAX_EE_FORWARD" --fastq_minlen 150 --fastq_maxlen "$MAX_LENGTH_FORWARD" --fastq_maxns 0 --fastaout "${s}.notmergedforward.filtered.fasta" --fasta_width 0

    # Dereplication & Move

    vsearch --derep_fulllength "${s}.filtered.fasta" --strand both --output "${s}.derep.merged.fasta" --sizeout --uc "${s}.derep.merged.uc" --relabel "${s}." --fasta_width 0


    mv *derep.merged* merged


    vsearch --derep_fulllength "${s}.notmergedforward.filtered.fasta" --strand both --output "${s}.derep.notmergedforward.fasta" --sizeout --uc "${s}.derep.notmergedforward.uc" --relabel "${s}." --fasta_width 0


    mv *derep.notmergedforward* unmergedf

done

# 5. Process Merged Reads

cd merged || exit

cat *derep*fasta > all.merged.fasta
vsearch --derep_fulllength all.merged.fasta --minuniquesize "$MIN_SIZE" --sizein --sizeout --fasta_width 0 --uc all.merged.derep.uc --output all.merged.derep.fasta
vsearch --cluster_size all.merged.derep.fasta --threads "$THREADS" --id 0.98 --strand both --sizein --sizeout --fasta_width 0 --uc all.merged.preclustered.uc --centroids all.merged.preclustered.fasta
vsearch --uchime_denovo all.merged.preclustered.fasta --sizein --sizeout --fasta_width 0 --nonchimeras all.merged.denovo.nonchimeras.fasta
perl "$PATH_MAP_SCRIPT" all.merged.derep.fasta all.merged.preclustered.uc all.merged.denovo.nonchimeras.fasta > all.merged.nonchimeras.derep.fasta
perl "$PATH_MAP_SCRIPT" all.merged.fasta all.merged.derep.uc all.merged.nonchimeras.derep.fasta > all.merged.nonchimeras.fasta
vsearch --cluster_size all.merged.nonchimeras.fasta --threads "$THREADS" --id 0.97 --strand both --sizein --sizeout --fasta_width 0 --uc all.merged.clustered.uc --relabel "OTU_M_${LABEL}_" --centroids all.otus.merged.fasta --otutabout all.tab_otus.merged.txt
cd ..

# 6. Process Unmerged Forward Reads

cd unmergedf || exit
cat *derep*.fasta > all.unmergedf.fasta
vsearch --derep_fulllength all.unmergedf.fasta --minuniquesize "$MIN_SIZE" --sizein --sizeout --fasta_width 0 --uc all.unmergedf.derep.uc --output all.unmergedf.derep.fasta
vsearch --cluster_size all.unmergedf.derep.fasta --threads "$THREADS" --id 0.97 --strand both --sizein --sizeout --fasta_width 0 --uc all.unmergedf.preclustered.uc --centroids all.unmergedf.preclustered.fasta
vsearch --uchime_denovo all.unmergedf.preclustered.fasta --sizein --sizeout --fasta_width 0 --nonchimeras all.unmergedf.denovo.nonchimeras.fasta
perl "$PATH_MAP_SCRIPT" all.unmergedf.derep.fasta all.unmergedf.preclustered.uc all.unmergedf.denovo.nonchimeras.fasta > all.unmergedf.nonchimeras.derep.fasta
perl "$PATH_MAP_SCRIPT" all.unmergedf.fasta all.unmergedf.derep.uc all.unmergedf.nonchimeras.derep.fasta > all.unmergedf.nonchimeras.fasta
vsearch --cluster_size all.unmergedf.nonchimeras.fasta --threads "$THREADS" --id 0.97 --strand both --sizein --sizeout --fasta_width 0 --uc all.unmergedf.clustered.uc --relabel "OTU_UnmergedF_${LABEL}_" --centroids all.otus.unmergedf.fasta --otutabout all.tab_otus.unmergedf.txt
cd ..

# 7. Final Output Organization & BLAST

mkdir -p otu

cp merged/all.otus.merged.fasta unmergedf/all.otus.unmergedf.fasta ./otu

cp merged/all.tab_otus.merged.txt unmergedf/all.tab_otus.unmergedf.txt ./otu

cd otu/ || exit


cat all.tab_otus.merged.txt > merged.tab.txt  && tail -n +2 all.tab_otus.unmergedf.txt >> merged_plus_reverse.tab.txt

cat all.otus.merged.fasta all.otus.unmergedf.fasta > merged_plus_forward.fasta

grep "#" all.tab_otus.merged.txt > CHECK.TXT && grep "#" all.tab_otus.unmergedf.txt >> CHECK.TXT

rm -f ptrim_*.fastq *.notmergedforward.fastq *.notmergedforward.filtered.fasta *.merged.fastq *.filtered.fasta && rmdir unmergedr 2>/dev/null && echo "Cleanup done"


```

```R
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

```





