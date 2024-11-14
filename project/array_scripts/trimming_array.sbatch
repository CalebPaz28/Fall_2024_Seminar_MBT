#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=trimming_array_%J_stdout.txt
#SBATCH --error=trimming_array_%J_stderr.txt
#SBATCH --job-name=trimming_array
#SBATCH --array=1-8
# 

bash /scratch/biol726310/BIOL7263_Genomics/array_project/scripts/trimming_array.sh $(sed -n "${SLURM_ARRAY_TASK_ID}p" /scratch/biol726310/BIOL7263_Genomics/array_project/scripts/trimming_array.args)

