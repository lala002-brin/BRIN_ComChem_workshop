#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=32GB
#SBATCH --partition=workshop

#SBATCH --output=md_runmembrane.out
#SBATCH --error=md_runmembrane.err


FILE_INPUT=md_runmembrane.mcr

yasara -txt ${FILE_INPUT}
