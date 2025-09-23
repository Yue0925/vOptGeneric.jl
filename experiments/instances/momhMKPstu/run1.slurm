#!/bin/bash
#SBATCH --job-name=1run

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=1run_OUTPUT.txt
#SBATCH --error=1run_ERROR.txt


# for file in ./MOBKP/set3/*; do
#     echo "$file"
#     srun julia vOptMomkp.jl "$file"
# done

# srun julia latexWriter.jl
# srun julia vOptMomkp.jl ./MOBKP/set3/Wcollage-tube.DAT



for file in ./MOBKP/set3/*; do
    echo "$file"
    srun julia vOptMomkp.jl "$file"
done