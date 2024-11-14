# Define a unique output file name using SLURM_ARRAY_TASK_ID
OUTPUT_FILE="output_${SLURM_ARRAY_TASK_ID}.sam"

# alignment of raw reads (trimmed) to reference
hisat2 -x /scratch/biol726310/BIOL7263_Genomics/array_project/reference_seq/indexed_genome -1 $1 -2 $2 -S /scratch/biol726310/BIOL7263_Genomics/array_project/alignment_data/$OUTPUT_FILE

