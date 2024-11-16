# Define a unique output file name using SLURM_ARRAY_TASK_ID
OUTPUT_FOLDER="spades_assembly_${SLURM_ARRAY_TASK_ID}"

# RNA SPAdes
spades.py --rna --threads 24 -o /scratch/biol726310/BIOL7263_Genomics/array_project/assembly_data/$OUTPUT_FOLDER -1 /scratch/biol726310/BIOL7263_Genomics/array_project/unmapped_reads/$1 -2 /scratch/biol726310/BIOL7263_Genomics/array_project/unmapped_reads/$2
