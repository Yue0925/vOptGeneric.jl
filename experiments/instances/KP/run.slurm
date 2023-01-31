#!/bin/bash
#SBATCH --job-name=test_BOHKP

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=test_BOHKP_OUTPUT.txt
#SBATCH --error=test_BOHKP_ERROR.txt

methodes=("epsilon" "bb" "bc" "bc_rootRelax" "bb_EPB" "bc_EPB" "bc_rootRelaxEPB" "bc_rootRelaxCP" "bc_rootRelaxCPEPB")

for file in ./hardinstances_pisinger/*.csv; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        srun julia vOptHKP.jl "$file" $mthd
    done
done
