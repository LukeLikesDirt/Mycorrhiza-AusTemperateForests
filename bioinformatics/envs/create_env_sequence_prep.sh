#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1-00:00:00
#SBATCH --partition=day
#SBATCH --output=%x.%j.out

# Set constants
ENV_FILE="./sequence_prep.yml"

### Create the mamba environment ###############################################

echo "Creating the conda environment..."

if ! conda info --envs | grep -q 'sequence_prep'; then
  echo "Creating conda environment at: $(date)"
  mamba env create -f $ENV_FILE
else
  echo "Environment already exists!"
fi

echo "Finished creating conda environment at: $(date)"

# List installed packages #######################################################

# Activate the environment
source activate sequence_prep

# Print the versions of installed packages
echo "Installed package versions:"
conda list

# Deactivate the environment
conda deactivate