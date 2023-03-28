#!/bin/bash
#SBATCH --job-name=MaxCut1
#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --output=MaxCut1_OUTPUT.txt
#SBATCH --error=MaxCut1_ERROR.txt


methodes=("epsilon")

for file in ./instances/*; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        julia vOptMaxCut.jl "$file" $mthd
    done
done
