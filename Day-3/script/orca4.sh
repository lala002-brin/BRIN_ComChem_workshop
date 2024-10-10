#!/bin/bash
  
#SBATCH --job-name=NAMA-PERHITUNGAN
#SBATCH --partition=workshop
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --error=error.txt
#SBATCH --mem=16GB

## Required Modules
module load openmpi4/4.1.4
module load nuclear/orca/4.3.2

echo "start = date"

# Put your job here
/mgpfs/apps/nuclear/apps/orca/4.3.2/orca orca4.inp > hasil.out --oversubscribe
