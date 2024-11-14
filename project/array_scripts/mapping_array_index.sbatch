#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 32G
#SBATCH --output=mapping_array_index_%J_stdout.txt
#SBATCH --error=mapping_array_index_%J_stderr.txt
#SBATCH --job-name=mapping_array_index
# 

bash /scratch/biol726310/BIOL7263_Genomics/array_project/scripts/mapping_array_index.sh
