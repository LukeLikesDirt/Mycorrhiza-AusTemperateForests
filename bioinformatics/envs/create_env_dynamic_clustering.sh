#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1-00:00:00
#SBATCH --partition=day
#SBATCH --output=%x.%j.out

# Set constants
ENV_FILE="./dynamic_clustering.yml"

# Create the conda environment #################################################

echo "Create the conda environment."

if ! conda info --envs | grep -q 'dynamic_clustering'; then
  echo "Creating conda environment at: $(date)"
  mamba env create -f $ENV_FILE
else
  echo "Environment already exists."
fi

# Download DADA2 ################################################################

echo "Install DADA2."
Rscript -e "devtools::install_github('benjjneb/dada2', ref='v1.16')"

# List installed packages #######################################################

# Activate the environment
source activate dynamic_clustering

# Print the versions of installed packages
echo "Installed package versions:"
conda list

# Deactivate the environment
conda deactivate