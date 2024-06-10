#!/bin/bash
#SBATCH --job-name=AP_lambda

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=AP_lambdaOut.txt
#SBATCH --error=AP_lambdaErr.txt





methodes=("bb" "bb_EPB" "bc_rootRelax" "bc_rootRelaxEPB")
files=("./AP/AP_p-3_n-5_ins-1.dat" "./AP/AP_p-3_n-5_ins-10.dat" "./AP/AP_p-3_n-10_ins-1.dat" "./AP/AP_p-3_n-10_ins-10.dat" "./AP/AP_p-3_n-10_ins-2.dat" "./AP/AP_p-3_n-15_ins-1.dat" "./AP/AP_p-3_n-15_ins-10.dat" "./AP/AP_p-3_n-15_ins-2.dat" "./AP/AP_p-3_n-20_ins-1.dat" "./AP/AP_p-3_n-20_ins-2.dat" "./AP/AP_p-3_n-20_ins-3.dat")


for file in ${files[@]}; do
    for mthd in ${methodes[@]}; do
        echo $file " ... " $mthd
        srun ../../../../julia-1.7.3/bin/julia vOptBOAPheur.jl $file $mthd
    done
done
