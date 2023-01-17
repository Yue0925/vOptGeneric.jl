#!/bin/bash
#SBATCH --job-name=logN30

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=logN30.txt
#SBATCH --error=logN30.txt


srun julia vOptMomkp_complementaire.jl ./MOBKP/set3/Wcollage-tube.DAT 