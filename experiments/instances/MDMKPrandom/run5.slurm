#!/bin/bash
#SBATCH --job-name=aleaMDMKP_lambda_test2
#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=aleaMDMKP_lambda_test2.txt
#SBATCH --error=aleaMDMKP_lambda_test2.txt

methodes=("bc_rootRelax")

for file in ./instances/PCC/*; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        srun julia vOptMDMKP2.jl "$file" $mthd
    done
done