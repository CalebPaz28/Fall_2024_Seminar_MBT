# Define a unique output file name using SLURM_ARRAY_TASK_ID
OUTPUT_UNMAPPED_SORTED="unmapped_output_${SLURM_ARRAY_TASK_ID}.bam"
OUTPUT_FASTQ1="unmapped_output_${SLURM_ARRAY_TASK_ID}_1.fastq"
OUTPUT_FASTQ2="unmapped_output_${SLURM_ARRAY_TASK_ID}_2.fastq"


# Extract unmapped reads from the sorted BAM file
samtools view -b -f 12 $1 -o /scratch/biol726310/BIOL7263_Genomics/array_project/unmapped_reads/$OUTPUT_UNMAPPED_SORTED


# Convert the .BAM file to a format (using bedtools) that can be read by SPAdes such as fastq
bedtools bamtofastq -i $2 -fq /scratch/biol726310/BIOL7263_Genomics/array_project/unmapped_reads/$OUTPUT_FASTQ1 -fq2 /scratch/biol726310/BIOL7263_Genomics/array_project/unmapped_reads/$OUTPUT_FASTQ2


# Convert the .BAM file to a format (using samtools) that can be read by SPAdes such as fastq

samtools fastq $2 -1 /scratch/biol726310/BIOL7263_Genomics/array_project/unmapped_reads/$OUTPUT_FASTQ1 -2 /scratch/biol726310/BIOL7263_Genomics/array_project/unmapped_reads/$OUTPUT_FASTQ2
