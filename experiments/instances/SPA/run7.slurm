#!/bin/bash
#SBATCH --job-name=SPA_lambda

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=SPA_lambdaOut.txt
#SBATCH --error=SPA_lambdaErr.txt





methodes=("bb" "bb_EPB" "bc_rootRelax" "bc_rootRelaxEPB")
files=("./BOSPA/biosppnw08.txt" "./BOSPA/biosppnw10.txt" "./BOSPA/biosppnw12.txt" "./BOSPA/biosppnw19.txt" "./BOSPA/biosppnw20.txt")


for file in ${files[@]}; do
    for mthd in ${methodes[@]}; do
        echo $file " ... " $mthd
        srun ../../../../julia-1.7.3/bin/julia vOptBOSPAheur.jl $file $mthd
    done
done
