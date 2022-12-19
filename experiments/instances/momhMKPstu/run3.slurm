#!/bin/bash
#SBATCH --job-name=test_cuts

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=test_cuts_OUTPUT.txt
#SBATCH --error=test_cuts_ERROR.txt


for file in ./MOBKP/set3/*; do
    echo "$file"
    srun julia vOptMomkp_complementaire.jl "$file"
done

