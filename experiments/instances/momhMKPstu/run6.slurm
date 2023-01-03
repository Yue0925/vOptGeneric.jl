#!/bin/bash
#SBATCH --job-name=logN40

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=logN40.txt
#SBATCH --error=logN40.txt



srun julia vOptMomkp.jl ./MOBKP/set3/W7BI-tube1-asyn.DAT 
