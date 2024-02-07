#!/bin/bash
#SBATCH --job-name=aleaMDMKP_lambda_test
#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=aleaMDMKP_lambda_test.txt
#SBATCH --error=aleaMDMKP_lambda_test.txt

# methodes=("bc_rootRelax" "bc_rootRelaxEPB")

# for file in ./instances/PCC/*; do
#     for mthd in ${methodes[@]}; do
#         echo "$file ... " $mthd
#         srun julia vOptMDMKP2.jl "$file" $mthd
#     done
# done

methodes=("bc_rootRelaxCP" "bc_rootRelaxCPEPB")

for file in ./instances/PCC/*; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        srun julia vOptMDMKP.jl "$file" $mthd
    done
done
