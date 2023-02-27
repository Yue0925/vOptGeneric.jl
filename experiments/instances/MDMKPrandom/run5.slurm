#!/bin/bash
#SBATCH --job-name=alea_MDMKP5
#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=alea_MDMKP5_OUTPUT.txt
#SBATCH --error=alea_MDMKP5_ERROR.txt

methodes=("epsilon")

for file in ./instances/PCC/*; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        srun julia vOptMDMKP.jl "$file" $mthd
    done
done
