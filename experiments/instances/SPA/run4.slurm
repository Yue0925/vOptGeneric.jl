#!/bin/bash
#SBATCH --job-name=SPA_test4

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=SPA_test4.txt
#SBATCH --error=SPA_test4.txt

# srun julia latexWriter.jl


methodes=("bc_rootRelaxCP" "bc_rootRelaxCPEPB")


for file in ./BOSPA/*.txt; do
    for mthd in ${methodes[@]}; do
        echo "$file" " ... " $mthd
        srun julia vOptBOSPA.jl "$file" $mthd
    done
done
