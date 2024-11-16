#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem 64G
#SBATCH --output=assembly_array_%J_stdout.txt
#SBATCH --error=assembly_array_%J_stderr.txt
#SBATCH --job-name=assembly_array
#SBATCH --array=1-4
# 

bash /scratch/biol726310/BIOL7263_Genomics/array_project/scripts/assembly_array.sh $(sed -n "${SLURM_ARRAY_TASK_ID}p" /scratch/biol726310/BIOL7263_Genomics/array_project/scripts/assembly_array.args)
