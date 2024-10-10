#!/bin/bash
  
#SBATCH --job-name=Quantum-Chemistry-Calculation
#SBATCH --partition=workshop
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --error=error.txt
#SBATCH --mem=16GB

## Required Modules
module load openmpi4/4.1.4
module load nuclear/orca/6.0.0


echo "start = date"

# Put your job here
/mgpfs/apps/nuclear/apps/orca/6.0.0/orca orca6.inp > hasil.out --oversubscribe



echo "end = date"
