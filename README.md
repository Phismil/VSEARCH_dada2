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







