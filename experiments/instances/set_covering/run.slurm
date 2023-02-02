#!/bin/bash
#SBATCH --job-name=BOSCP

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=BOSCP_OUTPUT.txt
#SBATCH --error=BOSCP_ERROR.txt



for file in ./ORLib/*.txt; do
    echo "$file"
    srun julia vOptSCP.jl "$file"
done
