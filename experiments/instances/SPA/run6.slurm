#!/bin/bash
#SBATCH --job-name=SPA_lambda

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=SPA_lambda.txt
#SBATCH --error=SPA_lambda.txt

# srun julia latexWriter.jl


methodes=("bc_rootRelax" "bc_rootRelaxEPB")


for file in ./BOSPA/*.txt; do
    for mthd in ${methodes[@]}; do
        echo "$file" " ... " $mthd
        srun ../../../../julia-1.7.3/bin/julia vOptBOSPA2.jl "$file" $mthd
    done
done
