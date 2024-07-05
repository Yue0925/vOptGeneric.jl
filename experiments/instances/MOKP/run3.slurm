#!/bin/bash
#SBATCH --job-name=MOKP3

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=MOKP3Err.txt
#SBATCH --error=MOKP3Out.txt



methodes=("bc_rootRelax")
# "bc_rootRelaxCP" "bc_rootRelaxCPEPB" "bc_rootRelaxEPB" "epsilon"



for file in ./refKP/*.dat; do
    for mthd in ${methodes[@]}; do
        echo "$file" " ... " $mthd
        srun ../../../../julia-1.7.3/bin/julia vOptMOKPmixed.jl "$file" $mthd
    done
done
