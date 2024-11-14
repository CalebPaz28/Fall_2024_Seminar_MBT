## Downloading the DNA genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/994/745/GCA_002994745.2_RchiOBHm-V2/GCA_002994745.2_RchiOBHm-V2_genomic.fna.gz -P /scratch/biol726310/BIOL7263_Genomics/array_project/reference_seq

## Downloading trascripts
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/994/745/GCA_002994745.2_RchiOBHm-V2/GCA_002994745.2_RchiOBHm-V2_genomic.gff.gz -P /scratch/biol726310/BIOL7263_Genomics/array_project/reference_seq

## Uncompressing the files
gunzip /scratch/biol726310/BIOL7263_Genomics/array_project/reference_seq/*.gz
