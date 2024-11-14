# Define a unique output file name using SLURM_ARRAY_TASK_ID
OUTPUT_BAM="output_${SLURM_ARRAY_TASK_ID}.bam"
OUTPUT_SORT="sorted_output_${SLURM_ARRAY_TASK_ID}.bam"
OUTPUT_STATS="stats_output_${SLURM_ARRAY_TASK_ID}.txt"


# converting the SAM to BAM
samtools view -S -b /scratch/biol726310/BIOL7263_Genomics/array_project/alignment_data/$1 > /scratch/biol726310/BIOL7263_Genomics/array_project/alignment_data/$OUTPUT_BAM


# sorting the BAM file
samtools sort /scratch/biol726310/BIOL7263_Genomics/array_project/alignment_data/$2 -o /scratch/biol726310/BIOL7263_Genomics/array_project/alignment_data/$OUTPUT_SORT


# indexing the BAM (used for IGV) and should generate and .bai file
samtools index /scratch/biol726310/BIOL7263_Genomics/array_project/alignment_data/$1


# using samtools flagstat for checking the alignment
samtools flagstat /scratch/biol726310/BIOL7263_Genomics/array_project/alignment_data/$3 > /scratch/biol726310/BIOL7263_Genomics/array_project/alignment_data/$OUTPUT_STATS


# using qualimap to check the alignment (doesnt work as well with the array)
qualimap bamqc -outdir /scratch/biol726310/BIOL7263_Genomics/array_project/alignment_data/bamqc -bam /scratch/biol726310/BIOL7263_Genomics/array_project/alignment_data/$3 -gff /scratch/biol726310/BIOL7263_Genomics/array_project/reference_seq/GCA_002994745.2_RchiOBHm-V2_genomic.gff
