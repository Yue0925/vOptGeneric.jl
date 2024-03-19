#!/bin/bash
#SBATCH --job-name=lambda_test

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=lambda_test.txt
#SBATCH --error=lambda_test.txt




for file in ./MOBKP/set3/*; do
    echo "$file"
    srun julia vOptMomkp5.jl "$file"
done
