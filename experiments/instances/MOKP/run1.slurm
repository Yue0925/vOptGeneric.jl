#!/bin/bash
#SBATCH --job-name=MOKP1

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=MOKP1Err.txt
#SBATCH --error=MOKP1Out.txt



methodes=("bb" "bb_EPB")
#   "bc_EPB" "bc_rootRelax" "bc_rootRelaxEPB" "bc_rootRelaxCP" "bc_rootRelaxCPEPB" 



for file in ./refKP/*.dat; do
    for mthd in ${methodes[@]}; do
        echo "$file" " ... " $mthd
        srun ../../../../julia-1.7.3/bin/julia vOptMOKP.jl "$file" $mthd
    done
done
