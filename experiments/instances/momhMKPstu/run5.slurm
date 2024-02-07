#!/bin/bash
#SBATCH --job-name=rootrelax_param2

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=rootrelax_param2.txt
#SBATCH --error=rootrelax_param2.txt


for file in ./MOBKP/set3/*; do
    echo "$file"
    srun julia vOptMomkp_complementaire.jl "$file"
done

# srun julia latexWriter.jl