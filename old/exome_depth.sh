#!/usr/bin/env bash
# The name to show in queue lists for this job:
##SBATCH -J exomeDetph_R.sh

# Number of desired cpus:
#SBATCH --ntasks=1

# Amount of RAM needed for this job:
#SBATCH --mem=20gb

# The time the job will be running:
#SBATCH --time=16:00:00


# To use GPUs you have to request them:
##SBATCH --gres=gpu:1

# If you need nodes with special features uncomment the desired constraint line:
##SBATCH --constraint=slim
##SBATCH --constraint=bigmem
##SBATCH --constraint=cal

# Set output and error files
#SBATCH --error=job.exomeDetph_R.err
#SBATCH --output=job.exomeDetph_R.out



# programa a ejecutar, con sus argumentos:
module load R-Bioconductor/3.0.1


time Rscript exome_depth.R samples3.txt trusight_one_Covered2.bed
