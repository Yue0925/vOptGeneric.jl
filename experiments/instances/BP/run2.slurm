#!/bin/bash
#SBATCH --job-name=BOBP

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=BOBP_OUTPUT.txt
#SBATCH --error=BOBP_ERROR.txt

methodes=("epsilon" "bb" "bc" "bc_rootRelax" "bb_EPB" "bc_EPB" "bc_rootRelaxEPB" "bc_rootRelaxCP" "bc_rootRelaxCPEPB")


for file in ./mipLib/*.mps; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        srun julia vOptBP.jl "$file" $mthd
    done
done
