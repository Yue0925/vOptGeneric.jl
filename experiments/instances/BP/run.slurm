#!/bin/bash
#SBATCH --job-name=mipLib

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=mipLib_OUTPUT.txt
#SBATCH --error=mipLib_ERROR.txt


for file in ./mipLib/*.mps; do
    echo "$file"
    srun julia parserBP.jl "$file"
done