#!/bin/bash
#SBATCH --job-name=BOSCP_rootrelax_param
#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=BOSCP_rootrelax_param.txt
#SBATCH --error=BOSCP_rootrelax_param.txt

methodes=("bc_rootRelax" "bc_rootRelaxEPB")

for file in ./instances/*; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        srun julia vOptSCP.jl "$file" $mthd
    done
done
