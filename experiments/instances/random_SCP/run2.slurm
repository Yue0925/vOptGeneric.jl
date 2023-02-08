#!/bin/bash
#SBATCH --job-name=alea_BOSCP2
#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=alea_BOSCP2_OUTPUT.txt
#SBATCH --error=alea_BOSCP2_ERROR.txt

methodes=("bb" "bc" "bb_EPB" "bc_EPB")

for file in ./instances/*; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        srun julia vOptSCP.jl "$file" $mthd
    done
done


for file in ./instances/*; do
    echo "$file ... "
    srun julia vOptSCP.jl "$file" epsilon
done
