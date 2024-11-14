## Using trimming software to initiate the trimming
trim_galore --paired --fastqc --gzip --cores 4 --length 100 /scratch/biol726310/BIOL7263_Genomics/array_project/raw_data/$1 /scratch/biol726310/BIOL7263_Genomics/array_project/raw_data/$2 --basename trimmed_reads -o /scratch/biol726310/BIOL7263_Genomics/array_project/trimmed_reads
