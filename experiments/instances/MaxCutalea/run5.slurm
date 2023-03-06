#!/bin/bash
#SBATCH --job-name=MaxCut5
#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=MaxCut5_OUTPUT.txt
#SBATCH --error=MaxCut5_ERROR.tx
methodes=("epsilon")

for file in ./instances/*; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        srun julia vOptMaxCut.jl "$file" $mthd
    done
done
