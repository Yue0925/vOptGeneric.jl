#!/bin/bash
#SBATCH --job-name=aleaBOSCP
#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=aleaBOSCP_OUTPUT.txt
#SBATCH --error=aleaBOSCP_ERROR.txt

methodes=("bc_rootRelaxCPEPB" "bc_rootRelaxEPB" "bc_rootRelaxCP" "bc_rootRelax")

for file in ./instances/*; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        srun julia vOptSCP.jl "$file" $mthd
    done
done
