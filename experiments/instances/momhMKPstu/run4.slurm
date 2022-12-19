#!/bin/bash
#SBATCH --job-name=testCPRootRelax

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=testCPRootRelax_OUTPUT.txt
#SBATCH --error=testCPRootRelax_ERROR.txt


for file in ./MOBKP/set3/*; do
    echo "$file"
    srun julia vOptMomkp4.jl "$file"
done

