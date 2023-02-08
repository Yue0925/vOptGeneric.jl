#!/bin/bash
#SBATCH --job-name=alt_BOMDKP4

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=alt_BOMDKP4_OUTPUT.txt
#SBATCH --error=alt_BOMDKP4_ERROR.txt

methodes=("bb" "bc" "bb_EPB" "bc_EPB")

for file in ./ORLib/*.txt; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        srun julia vOptMDKP_alt.jl "$file" $mthd
    done
done
