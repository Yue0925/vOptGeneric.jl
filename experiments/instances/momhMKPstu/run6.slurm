#!/bin/bash
#SBATCH --job-name=inst1

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=inst1.txt
#SBATCH --error=inst1.txt





srun julia vOptMomkp_complementaire.jl ./MOBKP/set3/W7BI-tube1-asyn.DAT 

