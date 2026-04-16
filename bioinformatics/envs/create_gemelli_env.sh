#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1-00:00:00
#SBATCH --partition=day
#SBATCH --output=%x.%j.out

ENV_FILE="gemelli_env.yml"

### Create the mamba environment ###############################################

echo "Creating conda environment..."

if ! conda info --envs | grep -q 'gemelli_env'; then
  echo "Creating conda environment at: $(date)"
  mamba env create -f $ENV_FILE
else
  echo "Environment already exists."
fi

echo "Finished creating conda environment at: $(date)"

### List installed packages #####################################################

source activate gemelli_env

echo "Installed package versions:"
conda list

conda deactivate
