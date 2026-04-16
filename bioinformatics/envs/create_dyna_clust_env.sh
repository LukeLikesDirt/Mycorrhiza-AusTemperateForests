#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=1-00:00:00
#SBATCH --partition=day
#SBATCH --output=%x.%j.out

ENV_FILE="dyna_clust_env.yml"

### Create the mamba environment ###############################################

echo "Create the conda environment."

if ! conda info --envs | grep -q 'dyna_clust_env'; then
  echo "Creating conda environment at: $(date)"
  mamba env create -f $ENV_FILE
else
  echo "Environment already exists."
fi

### Install DADA2 ###############################################################

echo "Install DADA2."
Rscript -e "devtools::install_github('benjjneb/dada2', ref='v1.16')"

