#!/bin/bash
#SBATCH --job-name=alea_MDMKP4
#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=alea_MDMKP4_OUTPUT.txt
#SBATCH --error=alea_MDMKP4_ERROR.txt

methodes=("bc_rootRelaxEPB" "bc_rootRelax")

for file in ./instances/PCC/*; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        srun julia vOptMDMKP.jl "$file" $mthd
    done
done
