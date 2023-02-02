#!/bin/bash
#SBATCH --job-name=BOBP2

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=BOBP2_OUTPUT.txt
#SBATCH --error=BOBP2_ERROR.txt


methodes=("bb" "bc" "bb_EPB" "bc_EPB")

for file in ./mipLib/*.mps; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        srun julia vOptBP.jl "$file" $mthd
    done
done



for file in ./mipLib/*.mps; do
    echo "$file ... " 
    srun julia vOptBP.jl "$file" epsilon
done
