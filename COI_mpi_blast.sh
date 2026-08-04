#!/bin/bash

#for f in specific_OTUs shared_OTUs specific_ASVs shared_ASVs; do
#    mkdir -p "$f" && mv "$f.fasta" "$f/" && cp mpi_blast.sh "$f/"
#done




# ── Config ────────────────────────────────────────────────────────────────────
DB="/mnt/lustre/users/aemami-khoyi/ref_db/coi_2026.fasta"
PATTERN="*.part-*.fasta"
JOBS_PER_NODE=6
THREADS=4

# ── Get nodes from PBS ────────────────────────────────────────────────────────
NODE_LIST=($(cat $PBS_NODEFILE | sort -u))
NUM_NODES=${#NODE_LIST[@]}
echo "Found $NUM_NODES nodes: ${NODE_LIST[@]}"
WORKDIR=$(pwd)

# ── Collect input files ───────────────────────────────────────────────────────
FILES=($(ls "$WORKDIR"/$PATTERN 2>/dev/null))
TOTAL_FILES=${#FILES[@]}
(( TOTAL_FILES == 0 )) && { echo "No files found matching: $PATTERN"; exit 1; }
echo "Found $TOTAL_FILES files to process"

FILES_PER_NODE=$(( (TOTAL_FILES + NUM_NODES - 1) / NUM_NODES ))

# ── Local setup ───────────────────────────────────────────────────────────────
module add chpc/BIOMODULES blast gnu-parallel/20200322 2>/dev/null
mkdir -p ~/.parallel && touch ~/.parallel/will-cite

START_TIME=$(date +%s)

# ── Dispatch BLAST jobs to each node ─────────────────────────────────────────
for (( i = 0; i < NUM_NODES; i++ )); do
    NODE="${NODE_LIST[$i]}"
    START=$(( i * FILES_PER_NODE ))
    COUNT=$(( i == NUM_NODES - 1 ? TOTAL_FILES - START : FILES_PER_NODE ))
    BATCH=("${FILES[@]:$START:$COUNT}")
    echo "Sending ${#BATCH[@]} files to $NODE"

    ssh "$NODE" "
        . /etc/profile.d/modules.sh
        module add chpc/BIOMODULES blast gnu-parallel/20200322 2>/dev/null
        mkdir -p ~/.parallel && touch ~/.parallel/will-cite
        cd '$WORKDIR'
        parallel -j $JOBS_PER_NODE \
            '[ -f {} ] && blastn -query {} -db \"$DB\" \
             -out {}.blastn.out -num_threads $THREADS \
             -outfmt 6 -max_target_seqs 3 -evalue 0.00001 \
             -perc_identity 90 2>/dev/null' \
            ::: ${BATCH[*]}
    " &
done

wait

# ── Combine and report ────────────────────────────────────────────────────────
cat *.blastn.out > combined_blast_results.out 2>/dev/null
COMPLETED=$(ls *.blastn.out 2>/dev/null | wc -l)
DURATION=$(( ($(date +%s) - START_TIME) / 60 ))

echo "Done in ${DURATION} min | Completed: $COMPLETED / $TOTAL_FILES files"
