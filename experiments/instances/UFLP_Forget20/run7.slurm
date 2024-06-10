#!/bin/bash
#SBATCH --job-name=UFLP7

#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks=1 
#SBATCH --mail-type=end
#SBATCH --mail-user=yue.zhang@lipn.univ-paris13.fr
#SBATCH --partition=COMPUTE-SHORT
#SBATCH --output=UFLP7Out.txt
#SBATCH --error=UFLP7Err.txt





methodes=("bb" "bb_EPB" "bc_rootRelax" "bc_rootRelaxEPB")
files=("./instance/Forget20-UFLP_10_3_1-1000_1001-2000_2box_1_1.raw" "./instance/Forget20-UFLP_10_3_1-1000_1001-2000_2box_1_10.raw" "./instance/Forget20-UFLP_10_3_1-1000_1001-2000_2box_1_2.raw" "./instance/Forget20-UFLP_10_3_1-1000_1001-2000_2box_1_3.raw" "./instance/Forget20-UFLP_10_3_1-1000_1001-2000_2box_1_4.raw" "./instance/Forget20-UFLP_12_3_1-1000_1001-2000_2box_1_4.raw" "./instance/Forget20-UFLP_12_3_1-1000_1001-2000_2box_1_3.raw" "./instance/Forget20-UFLP_12_3_1-1000_1001-2000_2box_1_2.raw" "./instance/Forget20-UFLP_12_3_1-1000_1001-2000_2box_1_1.raw")


for file in ${files[@]}; do
    for mthd in ${methodes[@]}; do
        echo $file " ... " $mthd
        srun ../../../../julia-1.7.3/bin/julia vOptUFLP_Forget20heur.jl $file $mthd
    done
done
