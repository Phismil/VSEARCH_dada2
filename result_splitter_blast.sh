for lower in 100 151 201 251 281; do

    upper=$((lower + 50))

    seqkit seq -m $lower -M $upper all.asvs.fasta > bin_${lower}_${upper}.fasta

    mkdir -p bin_${lower}_${upper}

    mv bin_${lower}_${upper}.fasta bin_${lower}_${upper}/

    cp COI_mpi_blast.sh bin_${lower}_${upper}/

    cd bin_${lower}_${upper}

    perl /mnt/lustre/users/aemami-khoyi/fasta-splitter.pl --n-parts 100 bin_${lower}_${upper}.fasta

    bash COI_mpi_blast.sh

    cd ..

done

for dir in bin_*/; do

    if [[ -f "${dir}combined_blast_results.out" ]]; then

        cat "${dir}combined_blast_results.out" >> all_combined_blast_results.out

        echo "Added ${dir}combined_blast_results.out"

    else

        echo "WARNING: ${dir}combined_blast_results.out not found — skipping"

    fi

done