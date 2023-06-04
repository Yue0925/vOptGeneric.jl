#!/bin/bash
#SBATCH --job-name=alea_MDMKP_basic_algo
#SBATCH --gres=gpu:1
#SBATCH --qos=qos_gpu-t4
#SBATCH --output=./alea_MDMKP_basic_algo.out
#SBATCH --error=./alea_MDMKP_basic_algo.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr

source ~/.bashrc

methodes=("epsilon" "bb" "bc" "bb_EPB" "bc_EPB")


for file in ./instances/PCC/*; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        julia vOptMDMKP.jl "$file" $mthd
    done
done

methodes=("bc_rootRelaxCPEPB" "bc_rootRelaxCP" "bc_rootRelaxEPB" "bc_rootRelax")

for file in ./instances/PCC/*; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        julia vOptMDMKP.jl "$file" $mthd
    done
done
