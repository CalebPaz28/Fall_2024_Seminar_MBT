#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 32G
#SBATCH --output=mapping_array_%J_stdout.txt
#SBATCH --error=mapping_array_%J_stderr.txt
#SBATCH --job-name=mapping_array
#SBATCH --array=1-4
# 

bash /scratch/biol726310/BIOL7263_Genomics/array_project/scripts/mapping_array.sh $(sed -n "${SLURM_ARRAY_TASK_ID}p" /scratch/biol726310/BIOL7263_Genomics/array_project/scripts/mapping_array.args)

