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

# methodes=("bb" "bc")

# for file in ./instances/PCC/*; do
#     for mthd in ${methodes[@]}; do
#         echo "$file ... " $mthd
#         srun julia vOptMDMKP.jl "$file" $mthd
#     done
# done


methodes=("bb" "bc" "bb_EPB" "bc_EPB" "bc_rootRelax" "bc_rootRelaxEPB" "bc_rootRelaxCP" "bc_rootRelaxCPEPB")

for mthd in ${methodes[@]}; do
    echo " ... " $mthd
    srun julia vOptMDMKP.jl ./instances/PCC/MDMKP_n30_m2_q2_α0.25 $mthd
done

