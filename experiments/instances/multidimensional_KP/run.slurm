#!/bin/bash
#SBATCH --job-name=alt_BOMDKP

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=alt_BOMDKP_OUTPUT.txt
#SBATCH --error=alt_BOMDKP_ERROR.txt

methodes=("epsilon" "bb" "bc" "bc_rootRelax" "bb_EPB" "bc_EPB" "bc_rootRelaxEPB" "bc_rootRelaxCP" "bc_rootRelaxCPEPB")

for file in ./ORLib/*.txt; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        srun julia vOptMDKP_alt.jl "$file" $mthd
    done
done