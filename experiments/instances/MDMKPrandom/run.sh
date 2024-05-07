#!/bin/bash


# for file in ./instances/PCC/*; do
#     echo "$file"
#     julia vOptMDMKP.jl "$file" epsilon
# done


methodes=("bb" "bc" "bb_EPB" "bc_EPB" "bc_rootRelax" "bc_rootRelaxEPB" "bc_rootRelaxCP" "bc_rootRelaxCPEPB")

for mthd in ${methodes[@]}; do
    echo " ... " $mthd
    julia vOptMDMKP.jl ./instances/PCC/MDMKP_n30_m2_q2_α0.25 $mthd
done
