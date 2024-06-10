#!/bin/bash
#SBATCH --job-name=testSPA_heur

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=testSPA_heur_OUTPUT.txt
#SBATCH --error=testSPA_heur_ERROR.txt


# for file in ./BOSPA/*.txt; do
#     echo "$file"
#     srun julia vOptBOSPA.jl "$file"
# done

# srun julia latexWriter.jl




methodes=("bb" "bc_rootRelax" "bc_rootRelaxEPB")
# "bc" "bb_EPB" "bc_EPB" "bc_rootRelax" "bc_rootRelaxEPB" "bc_rootRelaxCP" "bc_rootRelaxCPEPB" 

for mthd in ${methodes[@]}; do
    echo " ... " $mthd
    srun ../../../../julia-1.7.3/bin/julia vOptBOSPA.jl ./BOSPA/biosppnw12.txt $mthd
done
