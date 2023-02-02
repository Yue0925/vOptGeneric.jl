#!/bin/bash
#SBATCH --job-name=alt_BOSCP
#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=alt_BOSCP_OUTPUT.txt
#SBATCH --error=alt_BOSCP_ERROR.txt

methodes=("bc_rootRelaxCPEPB" "bc_rootRelaxEPB" "bc_rootRelaxCP" "bc_rootRelax")

for file in ./ORLib/*.txt; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        srun julia vOptSCP_alt.jl "$file" $mthd
    done
done
