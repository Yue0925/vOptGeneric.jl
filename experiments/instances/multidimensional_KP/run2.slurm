#!/bin/bash
#SBATCH --job-name=BOMDKP_new_exp2

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=BOMDKP_new_exp2_OUTPUT.txt
#SBATCH --error=BOMDKP_new_exp2_ERROR.txt

methodes=("bc_rootRelaxCPEPB" "bc_rootRelaxCP")

for file in ./ORLib/*.txt; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        srun julia vOptMDKP_alt.jl "$file" $mthd
    done
done