# VSEARCH
## COI

rm -rf merged unmergedf unmergedr otu

PRIMER_FWD="GGWACWGGWTGAACWGTWTAYCCYCC"
PRIMER_REV="TANACYTCNGGRTGNCCRAARAAYCA"

MARKER_LENGTH="317"
LABEL="COI"
MAX_EE="1"
MIN_LENGTH_MERGED="200"
MIN_OVERLAP="35"
MAX_DIF="35"
MAX_EE_FORWARD="1"
MAX_LENGTH_FORWARD="300"
THREADS="24"
MIN_SIZE="1"
UNOISE_ALPHA="5.0"

PATH_MAP_SCRIPT="/mnt/lustre/users/aemami-khoyi/map.pl"

module purge
module add chpc/BIOMODULES
module add cutadapt/3.4
module add vsearch/2.18.0

# ==============================================================================
# 1. Primer Trimming
# ==============================================================================

for R1 in *_R1.fastq; do

    [ -e "$R1" ] || continue

    R2=${R1//R1.fastq/R2.fastq}

    cutadapt -j 0 -g "$PRIMER_FWD" -G "$PRIMER_REV" --discard-untrimmed -o "ptrim_${R1}" -p "ptrim_${R2}" "$R1" "$R2"

done

# ==============================================================================
#  Merging & Filtering
# ==============================================================================

mkdir -p merged unmergedf unmergedr

for f in ptrim_*_R1.fastq; do

    [ -e "$f" ] || continue

    r=$(sed -e "s/_R1./_R2./" <<< "$f")

    s=$(cut -d_ -f2,3 <<< "$f")

    echo "Merging Sample: $s"

    vsearch --fastq_mergepairs "$f" --reverse "$r" --fastq_minovlen "$MIN_OVERLAP" --fastq_maxdiffs "$MAX_DIF" \
            --fastqout "${s}.merged.fastq" --fastq_eeout \
            --fastqout_notmerged_fwd "${s}.notmergedforward.fastq" \
            --fastqout_notmerged_rev "${s}.notmergedreverse.fastq"

    rm -rf "${s}.notmergedreverse.fastq"

    # Filtering

    vsearch --fastq_filter "${s}.merged.fastq" --fastq_maxee "$MAX_EE" --fastq_minlen "$MIN_LENGTH_MERGED" \
            --fastq_maxlen "$MARKER_LENGTH" --fastq_maxns 0 --fastaout "${s}.filtered.fasta" --fasta_width 0

    vsearch --fastq_filter "${s}.notmergedforward.fastq" --fastq_maxee "$MAX_EE_FORWARD" --fastq_minlen 150 \
            --fastq_maxlen "$MAX_LENGTH_FORWARD" --fastq_maxns 0 --fastaout "${s}.notmergedforward.filtered.fasta" --fasta_width 0

    # Dereplication & Move

    vsearch --derep_fulllength "${s}.filtered.fasta" --strand both --output "${s}.derep.merged.fasta" \
            --sizeout --uc "${s}.derep.merged.uc" --relabel "${s}." --fasta_width 0

    mv *derep.merged* merged

    vsearch --derep_fulllength "${s}.notmergedforward.filtered.fasta" --strand both --output "${s}.derep.notmergedforward.fasta" \
            --sizeout --uc "${s}.derep.notmergedforward.uc" --relabel "${s}." --fasta_width 0

    mv *derep.notmergedforward* unmergedf

done

rm -rf unmergedr

# ==============================================================================
#  Denoise Merged
# ==============================================================================

find . -maxdepth 1 -type f \( -name '*merge*' -o -name '*filtered*' -o -name 'ptrim_*' \) -delete

cd merged || exit

cat *derep*.fasta > all.merged.fasta

cp all.merged.fasta input.fasta

vsearch --cluster_unoise input.fasta \
        --threads "$THREADS" --minsize "$MIN_SIZE" --unoise_alpha "$UNOISE_ALPHA" \
        --sizein --sizeout --centroids denoised.fasta --uc denoised.uc

vsearch --uchime3_denovo denoised.fasta \
        --sizein --sizeout --nonchimeras unoise3_non_chimera_output.fasta

perl "$PATH_MAP_SCRIPT" input.fasta denoised.uc unoise3_non_chimera_output.fasta > ASVs_unoise3.fasta

vsearch --cluster_size ASVs_unoise3.fasta \
        --threads "$THREADS" --id 1 --strand both \
        --sizein --sizeout --fasta_width 0 \
        --uc denoise_clustered.uc \
        --relabel "ASV_merged_${LABEL}_" \
        --centroids ASVs_merged.fasta \
        --otutabout ASVs_tab_merged.tsv

cd ..

# ==============================================================================
# 4. Denoise Unmerged Forward
# ==============================================================================

cd unmergedf || exit

cat *derep*.fasta > all.unmergedf.fasta

cp all.unmergedf.fasta input.fasta

vsearch --cluster_unoise input.fasta \
        --threads "$THREADS" --minsize "$MIN_SIZE" --unoise_alpha "$UNOISE_ALPHA" \
        --sizein --sizeout --centroids denoised.fasta --uc denoised.uc

vsearch --uchime3_denovo denoised.fasta \
        --sizein --sizeout --nonchimeras unoise3_non_chimera_output.fasta

perl "$PATH_MAP_SCRIPT" input.fasta denoised.uc unoise3_non_chimera_output.fasta > ASVs_unoise3.fasta

vsearch --cluster_size ASVs_unoise3.fasta \
        --threads "$THREADS" --id 1 --strand both \
        --sizein --sizeout --fasta_width 0 \
        --uc denoise_clustered.uc \
        --relabel "ASV_forward_${LABEL}_" \
        --centroids ASVs_onlyf.fasta \
        --otutabout ASVs_tab_onlyf.tsv

cd ..