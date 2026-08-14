#!/usr/bin/env bash

# ==============================================================================
# MMseqs2 ASV/OTU comparison + parallel BLAST preparation
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Run MMseqs2 search - top hit only
# ------------------------------------------------------------------------------

mmseqs easy-search ASVs.fasta all.otus.merged.fasta results.tsv tmp/ \
  --min-seq-id 0.99 -c 0.9 --cov-mode 0 --threads 8 --search-type 3 \
  --max-seqs 1 \
  --format-output "query,target,pident,qcov,tcov,alnlen,mismatch,evalue,bits,qseq,tseq"

# ------------------------------------------------------------------------------
# 2. Get matched/unmatched IDs
# ------------------------------------------------------------------------------

cut -f1 results.tsv | sort -u > hit_asvs.txt
cut -f2 results.tsv | sort -u > hit_otus.txt

seqkit seq -n ASVs.fasta            | awk '{print $1}' | sort -u > all_asvs.txt
seqkit seq -n all.otus.merged.fasta | awk '{print $1}' | sort -u > all_otus.txt

comm -23 all_asvs.txt hit_asvs.txt > specific_asvs.txt
comm -23 all_otus.txt hit_otus.txt > specific_otus.txt

# ------------------------------------------------------------------------------
# 3. Extract shared/specific sequences
# ------------------------------------------------------------------------------

seqkit grep -f hit_asvs.txt ASVs.fasta > shared_ASVs.fasta
seqkit grep -f hit_otus.txt all.otus.merged.fasta > shared_OTUs.fasta
seqkit grep -f specific_asvs.txt ASVs.fasta > specific_ASVs.fasta
seqkit grep -f specific_otus.txt all.otus.merged.fasta > specific_OTUs.fasta

# ------------------------------------------------------------------------------
# 4. Stats
# ------------------------------------------------------------------------------

TOTAL_ASVS=$(seqkit stats -T ASVs.fasta            | awk 'NR==2{print $4}')
TOTAL_OTUS=$(seqkit stats -T all.otus.merged.fasta | awk 'NR==2{print $4}')
SHARED_ASVS=$(seqkit stats -T shared_ASVs.fasta    | awk 'NR==2{print $4}')
SHARED_OTUS=$(seqkit stats -T shared_OTUs.fasta    | awk 'NR==2{print $4}')
SPEC_ASVS=$(seqkit stats -T specific_ASVs.fasta    | awk 'NR==2{print $4}')
SPEC_OTUS=$(seqkit stats -T specific_OTUs.fasta    | awk 'NR==2{print $4}')
TOTAL_HITS=$(wc -l < results.tsv)

{
echo "============================================"
echo "         SEARCH STATS (TOP HIT ONLY)        "
echo "============================================"
echo "Total ASVs (query):              $TOTAL_ASVS"
echo "Total OTUs (target):             $TOTAL_OTUS"
echo "Total hits in results.tsv:       $TOTAL_HITS"
echo "--------------------------------------------"
echo "Shared ASVs  (matched):          $SHARED_ASVS / $TOTAL_ASVS"
echo "Shared OTUs  (matched):          $SHARED_OTUS / $TOTAL_OTUS"
echo "Specific to ASVs (no hit):       $SPEC_ASVS   / $TOTAL_ASVS"
echo "Specific to OTUs (never hit):    $SPEC_OTUS   / $TOTAL_OTUS"
echo "--------------------------------------------"
printf "ASV match rate:                  %.1f%%\n" "$(echo "$SHARED_ASVS $TOTAL_ASVS" | awk '{printf "%.1f", $1/$2*100}')"
printf "OTU match rate:                  %.1f%%\n" "$(echo "$SHARED_OTUS $TOTAL_OTUS" | awk '{printf "%.1f", $1/$2*100}')"
echo "============================================"
} | tee 1_reciprocal.txt

# ------------------------------------------------------------------------------
# 5. Cleanup temporary ID files
# ------------------------------------------------------------------------------

rm hit_asvs.txt hit_otus.txt all_asvs.txt all_otus.txt

# ------------------------------------------------------------------------------
# 6. Split each FASTA category and run MPI BLAST
# ------------------------------------------------------------------------------

for f in specific_OTUs shared_OTUs specific_ASVs shared_ASVs; do

    mkdir -p "$f"

    cp "$f.fasta" "$f/"
    cp mpi_blast.sh "$f/"

    cd "$f" || exit

    perl /mnt/lustre/users/aemami-khoyi/fasta-splitter.pl --n-parts 10 *.fasta

    bash mpi_blast.sh

    cd ..

done