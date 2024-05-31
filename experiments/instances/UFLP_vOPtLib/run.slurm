#!/bin/bash
#SBATCH --job-name=testUFLP

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=testUFLP_OUTPUT.txt
#SBATCH --error=testUFLP_ERROR.txt


for file in ./dataDidactic/*.txt; do
    echo "$file"
    srun julia uncapacitatedFacilityLocation.jl "$file"
done

for file in ./dataFernandez/*.txt; do
    echo "$file"
    srun julia uncapacitatedFacilityLocation.jl "$file"
done

# srun julia latexWriter.jl