# Ensure the arguments file exists
if [ ! -f "arguments.txt" ]; then
  echo "Error: arguments.txt not found!"
  exit 1
fi

# Loop through each line in the arguments file
while read -r INPUT_SAM OUTPUT_BAM OUTPUT_SORT; do
  echo "Processing $INPUT_SAM..."

  # Step 1: Convert SAM to BAM
  samtools view -S -b /scratch/biol726310/BIOL7263_Genomics/array_project/alignment_data/$INPUT_SAM > /scratch/biol726310/BIOL7263_Genomics/array_project/alignment_data/$OUTPUT_BAM
  
  # Step 2: Sort the BAM file
  samtools sort /scratch/biol726310/BIOL7263_Genomics/array_project/alignment_data/$OUTPUT_BAM -o /scratch/biol726310/BIOL7263_Genomics/array_project/alignment_data/$OUTPUT_SORT
  
  # Step 3: Index the sorted BAM file for IGV
  samtools index /scratch/biol726310/BIOL7263_Genomics/array_project/alignment_data/$OUTPUT_SORT
  
  # Step 4: Generate alignment statistics with flagstat
  samtools flagstat /scratch/biol726310/BIOL7263_Genomics/array_project/alignment_data/$OUTPUT_SORT > /scratch/biol726310/BIOL7263_Genomics/array_project/alignment_data/${OUTPUT_SORT%.bam}_stats.txt

  # Step 5: Quality control with Qualimap
  qualimap bamqc -outdir /scratch/biol726310/BIOL7263_Genomics/array_project/alignment_data/bamqc_${OUTPUT_SORT%.bam}
  -bam /scratch/biol726310/BIOL7263_Genomics/array_project/alignment_data/$OUTPUT_SORT
  -gff /scratch/biol726310/BIOL7263_Genomics/array_project/reference_seq/GCA_002994745.2_RchiOBHm-V2_genomic.gff

  echo "Completed processing for $INPUT_SAM"

done < arguments.txt