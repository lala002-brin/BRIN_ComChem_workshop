#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=16GB
#SBATCH --partition=workshop

#SBATCH --output=md_analyze.out
#SBATCH --error=md_analyze.err


FILE_INPUT=md_analyze.mcr

/mgpfs/apps/sample/apps/ys/yasara -txt ${FILE_INPUT}
