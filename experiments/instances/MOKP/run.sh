#!/bin/bash


methodes=("bb" "bc" "bc_rootRelax")
#   "bc_EPB" "bc_rootRelax" "bc_rootRelaxEPB" "bc_rootRelaxCP" "bc_rootRelaxCPEPB" 



# for file in ./refKP/*.dat; do
    for mthd in ${methodes[@]}; do
        echo " ./refKP/KP_p-3_n-10_ins-1.dat ... " $mthd
        # srun ../../../../julia-1.7.3/bin/julia vOptKP_Forget20.jl "$file" $mthd
        julia vOptMOKP.jl "./refKP/KP_p-3_n-10_ins-1.dat" $mthd
    done
# done
