#!/bin/bash
#SBATCH --job-name=KPFgt3

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=KPFgt3.txt
#SBATCH --error=KPFgt3.txt



methodes=("bc_rootRelax")
#   "bc_rootRelaxCP" "bc_rootRelaxCPEPB" "bc_rootRelaxEPB" "epsilon"



for file in ./KP/*.raw; do
    for mthd in ${methodes[@]}; do
        echo "$file" " ... " $mthd
        srun ../../../../julia-1.7.3/bin/julia vOptKP_Forgetmixed.jl "$file" $mthd
    done
done
