#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 32G
#SBATCH --output=nucleotide_db_%J_stdout.txt
#SBATCH --error=nucleotide_db_%J_stderr.txt
#SBATCH --job-name=nucleotide_db
# 

bash /scratch/biol726310/BIOL7263_Genomics/array_project/scripts/nucleotide_db_download.sh
