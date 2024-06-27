#!/bin/bash
#SBATCH --job-name=MOKP4

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=MOKP4Err.txt
#SBATCH --error=MOKP4Out.txt



methodes=("bc_rootRelaxCP" "bc_rootRelaxCPEPB")
# 



for file in ./refKP/*.dat; do
    for mthd in ${methodes[@]}; do
        echo "$file" " ... " $mthd
        srun ../../../../julia-1.7.3/bin/julia vOptMOKP.jl "$file" $mthd
    done
done
