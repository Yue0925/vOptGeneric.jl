#!/bin/bash
#SBATCH --job-name=rootrelax_param

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=rootrelax_param.txt
#SBATCH --error=rootrelax_param.txt


for file in ./MOBKP/set3/*; do
    echo "$file"
    srun ../../../../julia-1.7.3/bin/julia vOptMomkpMixed.jl "$file"
done

# srun julia latexWriter.jl