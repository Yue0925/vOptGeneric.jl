#!/bin/bash
#SBATCH --job-name=basic_param2

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=basic_param2_OUTPUT.txt
#SBATCH --error=basic_param2_ERROR.txt


# for file in ./MOBKP/set2/*; do
#     echo "$file"
#     srun julia vOptMomkp2.jl "$file"
# done

# srun julia latexWriter.jl


for file in ./MOBKP/set3/*; do
    echo "$file"
    srun julia vOptMomkp2.jl "$file"
done