#!/bin/bash


# methodes=("epsilon" "bb" "bc" "bb_EPB" "bc_EPB")

# for file in ./instances/*; do
#     for mthd in ${methodes[@]}; do
#         echo "$file ... " $mthd
#         julia vOptSCP.jl "$file" $mthd
#     done
# done

methodes=("bb" "bc" "bb_EPB" "bc_EPB" "bc_rootRelax" "bc_rootRelaxEPB" "bc_rootRelaxCP" "bc_rootRelaxCPEPB")

for mthd in ${methodes[@]}; do
    echo " ... " $mthd
    julia vOptSCP.jl ./instances/SCP_80_15_2 $mthd
done
