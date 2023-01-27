#!/bin/bash
#SBATCH --job-name=BOMDKP

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=test_BOMDKP_OUTPUT.txt
#SBATCH --error=test_BOMDKP_ERROR.txt



for file in ./ORLib/*.txt; do
    echo "$file"
    srun julia vOptMDKP.jl "$file"
done
