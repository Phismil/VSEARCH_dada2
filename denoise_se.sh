# Primers

rm -rf unmergedf

PATH_MAP_SCRIPT="/mnt/lustre/users/aemami-khoyi/map.pl"

PRIMER_FWD="GGWACWGGWTGAACWGTWTAYCCYCC"

MARKER_LENGTH="317"
LABEL="COI"
MAX_EE_FORWARD="1"
MAX_LENGTH_FORWARD="300"
THREADS="24"
MIN_SIZE="1"           # Retain all sequences including singletons
UNOISE_ALPHA="5.0"     # Higher than default (2.0) = more lenient = more ZOTUs retained

module purge
module add chpc/BIOMODULES
module add cutadapt/3.4
module add vsearch/2.18.0

mkdir -p unmergedf

for R1 in *_R1.fastq; do

    [ -e "$R1" ] || continue

    s=$(cut -d_ -f1,2 <<< "$R1")

    PTRIM="ptrim_${R1}"

    cutadapt -j 0 \
        -g "${PRIMER_FWD}" \
        --error-rate 0.15 \
        --discard-untrimmed \
        -o "${PTRIM}" \
        "$R1"

    vsearch --fastq_filter "${PTRIM}" \
            --fastq_maxee "$MAX_EE_FORWARD" \
            --fastq_minlen 100 \
            --fastq_maxlen "$MAX_LENGTH_FORWARD" \
            --fastq_maxns 0 \
            --fastaout "${s}.notmergedforward.filtered.fasta" \
            --fasta_width 0

    vsearch --derep_fulllength "${s}.notmergedforward.filtered.fasta" \
            --strand both \
            --output "${s}.derep.notmergedforward.fasta" \
            --sizeout \
            --uc "${s}.derep.notmergedforward.uc" \
            --relabel "${s}." \
            --fasta_width 0

    mv *derep.notmergedforward* unmergedf

done

cd unmergedf || exit

cat *derep*.fasta > all.unmergedf.fasta

cp all.unmergedf.fasta input.fasta

# Denoise with lenient settings

vsearch --cluster_unoise input.fasta \
        --threads "$THREADS" \
        --minsize "$MIN_SIZE" \
        --unoise_alpha "$UNOISE_ALPHA" \
        --sizein \
        --sizeout \
        --centroids denoised.fasta \
        --uc denoised.uc

# Chimera filtering

vsearch --uchime3_denovo denoised.fasta \
        --sizein \
        --sizeout \
        --nonchimeras unoise3_output.fasta

perl "$PATH_MAP_SCRIPT" input.fasta denoised.uc unoise3_output.fasta > ASVs_unoise3.fasta

vsearch --cluster_size ASVs_unoise3.fasta \
        --threads "$THREADS" \
        --id 1 \
        --strand both \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --uc denoise_clustered.uc \
        --relabel "ASV_${LABEL}_" \
        --centroids ASVs.fasta \
        --otutabout ASVs_tab.tsv