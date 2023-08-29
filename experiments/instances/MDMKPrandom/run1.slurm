#!/bin/bash
#SBATCH --job-name=aleaMDMKP_basic_param
#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=aleaMDMKP_basic_param.txt
#SBATCH --error=aleaMDMKP_basic_param.txt

methodes=("epsilon" "bb" "bc")

for file in ./instances/PCC/*; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        srun julia vOptMDMKP.jl "$file" $mthd
    done
done
