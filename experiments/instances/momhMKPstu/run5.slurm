#!/bin/bash
#SBATCH --job-name=inst1_cplex

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=inst1_cplex.txt
#SBATCH --error=inst1_cplex.txt





# srun julia vOptMomkp3.jl ./MOBKP/set3/W7BI-tube1-asyn.DAT 

srun julia vOptMomkp3.jl ./MOBKP/set3/W7BI-tube1-1800.DAT 
