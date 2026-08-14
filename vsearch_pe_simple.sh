# VSEARCH
## COI

PRIMER_FWD="GGWACWGGWTGAACWGTWTAYCCYCC"
PRIMER_REV="TANACYTCNGGRTGNCCRAARAAYCA"

Round="c2"

MARKER_LENGTH="317"
LABEL="COI"
MAX_EE="1"
MIN_LENGTH_MERGED="200"
MIN_OVERLAP="35"
MAX_DIF="35"
MAX_EE_FORWARD="1"
MAX_LENGTH_FORWARD="300"
THREADS="24"
MIN_SIZE="2"

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

for R1 in *_R1.fastq; do

[ -e "$R1" ] || continue

R2=${R1//R1.fastq/R2.fastq}

cutadapt -j 0 -g "$PRIMER_FWD" -G "$PRIMER_REV" --discard-untrimmed -o "ptrim_${R1}" -p "ptrim_${R2}" "$R1" "$R2"

done

# 4. Merging Pairs

mkdir -p merged unmergedf unmergedr

for f in ptrim_*_R1.fastq; do

[ -e "$f" ] || continue

r=$(sed -e "s/_R1.fastq/_R2.fastq/" <<< "$f")

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

mkdir -p "$Round"_otu

cp merged/all.otus.merged.fasta unmergedf/all.otus.unmergedf.fasta ./"$Round"_otu

cp merged/all.tab_otus.merged.txt unmergedf/all.tab_otus.unmergedf.txt ./"$Round"_otu