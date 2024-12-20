#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem 32G
#SBATCH --output=protein_nt_blasting_%J_stdout.txt
#SBATCH --error=protein_nt_blasting_%J_stderr.txt
#SBATCH --job-name=protein_nt_blasting
#SBATCH --array=1-4
# 

bash /scratch/biol726310/BIOL7263_Genomics/array_project/scripts/protein_nt_blasting.sh $(sed -n "${SLURM_ARRAY_TASK_ID}p" /scratch/biol726310/BIOL7263_Genomics/array_project/scripts/protein_nt_blasting.args)