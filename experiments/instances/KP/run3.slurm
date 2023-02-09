#!/bin/bash
#SBATCH --job-name=BOHKP_new_exp3

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=BOHKP_new_exp3_OUTPUT.txt
#SBATCH --error=BOHKP_new_exp3_ERROR.txt

methodes=("epsilon" "bc_rootRelaxEPB" "bc_rootRelax")


for file in ./hardinstances_pisinger/*.csv; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        srun julia vOptHKP.jl "$file" $mthd
    done
done
