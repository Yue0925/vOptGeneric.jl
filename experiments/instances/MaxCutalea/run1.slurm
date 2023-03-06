#!/bin/bash
#SBATCH --job-name=MaxCut1
#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=MaxCut1_OUTPUT.txt
#SBATCH --error=MaxCut1_ERROR.txt

methodes=("bb" "bc")

for file in ./instances/*; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        srun julia vOptMaxCut.jl "$file" $mthd
    done
done
