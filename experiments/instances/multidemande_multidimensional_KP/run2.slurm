#!/bin/bash
#SBATCH --job-name=BOMDMDKP

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=test_BOMDMDKP_OUTPUT.txt
#SBATCH --error=test_BOMDMDKP_ERROR.txt



for file in ./ORLib/*.txt; do
    echo "$file"
    srun julia vOptMDMDKP.jl "$file"
done
