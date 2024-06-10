#!/bin/bash
#SBATCH --job-name=KPFgt5

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=KPFgt5.txt
#SBATCH --error=KPFgt5.txt



methodes=("bc_rootRelax" "bc_rootRelaxEPB")
#   "bc_EPB"  "bc_rootRelaxCP" "bc_rootRelaxCPEPB" 



for file in ./KP/*.raw; do
    for mthd in ${methodes[@]}; do
        echo "$file" " ... " $mthd
        srun ../../../../julia-1.7.3/bin/julia vOptKP_Forget20_.jl "$file" $mthd
    done
done
