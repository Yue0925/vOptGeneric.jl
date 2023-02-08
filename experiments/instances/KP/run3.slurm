#!/bin/bash
#SBATCH --job-name=BOHKP3

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=BOHKP3_OUTPUT.txt
#SBATCH --error=BOHKP3_ERROR.txt

# methodes=("bb" "bc" "bb_EPB" "bc_EPB")


# for file in ./hardinstances_pisinger/*.csv; do
#     for mthd in ${methodes[@]}; do
#         echo "$file ... " $mthd
#         srun julia vOptHKP.jl "$file" $mthd
#     done
# done


for file in ./hardinstances_pisinger/*.csv; do
    echo "$file ... "
    srun julia vOptHKP.jl "$file" epsilon
done
