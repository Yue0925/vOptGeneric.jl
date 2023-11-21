#!/bin/bash
#SBATCH --job-name=BOSCP_lambda_test
#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=BOSCP_lambda_test.txt
#SBATCH --error=BOSCP_lambda_test.txt

methodes=("bc_rootRelax" "bc_rootRelaxEPB")

for file in ./instances/*; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        srun julia vOptSCP2.jl "$file" $mthd
    done
done
