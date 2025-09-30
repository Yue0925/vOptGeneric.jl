#!/bin/bash
#SBATCH --job-name=KPFgt2

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=KPFgt2.txt
#SBATCH --error=KPFgt2.txt



methodes=("bc_EPB") # "bc" 
#  "bc_rootRelax" "bc_rootRelaxEPB" "bc_rootRelaxCP" "bc_rootRelaxCPEPB" 



for file in ./KP/*.raw; do
    for mthd in ${methodes[@]}; do
        echo "$file" " ... " $mthd
        srun julia vOptKP_Forget20.jl "$file" $mthd
    done
done
