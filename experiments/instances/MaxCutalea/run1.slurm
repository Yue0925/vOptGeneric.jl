#!/bin/bash
#SBATCH --job-name=MaxCut_basic_param
#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=MaxCut_basic_param.txt
#SBATCH --error=MaxCut_basic_param.txt

methodes=("epsilon") # "bb" "bc"

for file in ./instances/*; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        srun julia vOptMaxCut.jl "$file" $mthd
    done
done
