#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=fastqc_array_%J_stdout.txt
#SBATCH --error=fastqc_array_%J_stderr.txt
#SBATCH --job-name=fastqc_array
#SBATCH --array=1-8
# 

bash /scratch/biol726310/BIOL7263_Genomics/array_project/scripts/fastqc_array.sh $(sed -n "${SLURM_ARRAY_TASK_ID}p" /scratch/biol726310/BIOL7263_Genomics/array_project/scripts/fastqc_array.args)
