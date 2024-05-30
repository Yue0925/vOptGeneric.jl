#!/bin/bash
#SBATCH --job-name=AP_test5

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=AP_test5.txt
#SBATCH --error=AP_test5.txt

# srun julia latexWriter.jl


methodes=("epsilon")


for file in ./AP/*.dat; do
    for mthd in ${methodes[@]}; do
        echo "$file" " ... " $mthd
        srun julia vOptBOAP.jl "$file" $mthd
    done
done
