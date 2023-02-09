#!/bin/bash
#SBATCH --job-name=BOMDMDKP_new_exp

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=BOMDMDKP_new_exp_OUTPUT.txt
#SBATCH --error=BOMDMDKP_new_exp_ERROR.txt

methodes=("bb" "bc" "bb_EPB" "bc_EPB")

for file in ./ORLib/*.txt; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        srun julia vOptMDMDKP_alt.jl "$file" $mthd
    done
done
