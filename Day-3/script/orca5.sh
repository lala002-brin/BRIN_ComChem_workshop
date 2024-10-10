#!/bin/bash
  
#SBATCH --job-name=OPTIMIZE
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --error=error.txt
#SBATCH --mem=16GB

## Required Modules
module load openmpi4/4.1.4
module load nuclear/orca/5.0.3

echo "start = date"

# Put your job here
/mgpfs/apps/nuclear/apps/orca/5.0.3/orca orca5.inp > hasil.out --oversubscribe



echo "end = date"
