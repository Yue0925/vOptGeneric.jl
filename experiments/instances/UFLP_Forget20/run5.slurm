#!/bin/bash
#SBATCH --job-name=UFLP5

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=UFLP5.txt
#SBATCH --error=UFLP5.txt

# srun julia latexWriter.jl


methodes=("epsilon")



for file in ./instance/*.raw; do
    for mthd in ${methodes[@]}; do
        echo "$file" " ... " $mthd
        srun julia vOptUFLP_Forget20.jl "$file" $mthd
    done
done
