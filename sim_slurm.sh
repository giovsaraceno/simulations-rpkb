#!/bin/bash
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc
#SBATCH --job-name=rpkb
#SBATCH --mem=256G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=gsaracen@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-50
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err

module load gcc
module load openmpi
module load r

Rscript parameters.R
split -l 33 parameters.txt split_files/parameters_
for file in split_files/*; do
  files_list+=($file)
done

ulimit -u 10000
srun ./sim_shell.sh ${files_list[$SLURM_ARRAY_TASK_ID - 1]}