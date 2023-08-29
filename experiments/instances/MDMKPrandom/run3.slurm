#!/bin/bash
#SBATCH --job-name=aleaMDMKP_rootrelax_param
#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=aleaMDMKP_rootrelax_param.txt
#SBATCH --error=aleaMDMKP_rootrelax_param.txt

methodes=("bc_rootRelax" "bc_rootRelaxEPB" "bc_rootRelaxCP" "bc_rootRelaxCPEPB")

for file in ./instances/PCC/*; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        srun julia vOptMDMKP.jl "$file" $mthd
    done
done
