#!/bin/bash
#SBATCH --job-name=UFLP4

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=UFLP4.txt
#SBATCH --error=UFLP4.txt

# srun julia latexWriter.jl


methodes=("bc_rootRelaxEPB")



for file in ./instance/*.raw; do
    for mthd in ${methodes[@]}; do
        echo "$file" " ... " $mthd
        srun ../../../../julia-1.7.3/bin/julia vOptUFLP_Forget20.jl "$file" $mthd
    done
done
