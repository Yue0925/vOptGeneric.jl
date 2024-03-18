#!/bin/bash
#SBATCH --job-name=BOSCP_basic_param
#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=BOSCP_basic_param.txt
#SBATCH --error=BOSCP_basic_param.txt

methodes=("epsilon") # "bb" "bc"

for file in ./instances/*; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        srun julia vOptSCP.jl "$file" $mthd
    done
done
