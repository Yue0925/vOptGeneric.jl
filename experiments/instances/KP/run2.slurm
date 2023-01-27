#!/bin/bash
#SBATCH --job-name=BOHKP

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=test_BOHKP_OUTPUT.txt
#SBATCH --error=test_BOHKP_ERROR.txt


for file in ./hardinstances_pisinger/*.csv; do
    echo "$file"
    srun julia vOptHKP.jl "$file"
done
