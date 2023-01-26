#!/bin/bash
#SBATCH --job-name=test_MDKP

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=test_MDKP_OUTPUT.txt
#SBATCH --error=test_MDKP_ERROR.txt



for file in ./ORLib/*.txt; do
    echo "$file"
    srun julia vOptMDKP.jl "$file"
done
