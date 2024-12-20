# Define a unique output file name using SLURM_ARRAY_TASK_ID
PRO_OUTPUT="protein_results_${SLURM_ARRAY_TASK_ID}.csv"
NT_OUTPUT="nt_results_${SLURM_ARRAY_TASK_ID}.csv"


# This is supposed to generate a .csv file for nucleotide search
blastn -query $1 -db /scratch/biol726310/BIOL7263_Genomics/array_project/databases/nucleotide_db/nt_viruses -out /scratch/biol726310/BIOL7263_Genomics/array_project/blast_results/$NT_OUTPUT -evalue 1e-5 -max_target_seqs 5 -num_threads 20 -outfmt "6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"


# This is supposed to generate a .csv file for protein search
diamond blastx -d /scratch/biol726310/BIOL7263_Genomics/array_project/databases/protein_db/virus_protein_db.dmnd -q $1 --outfmt 6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore --max-target-seqs 5 --threads 20 | sed 's/\t/,/g' > /scratch/biol726310/BIOL7263_Genomics/array_project/blast_results/$PRO_OUTPUT
