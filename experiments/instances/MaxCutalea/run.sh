#!/bin/bash


# for file in ./instances/*; do
#     echo "$file"
#     julia vOptMaxCut.jl "$file" epsilon
# done

methodes=("bc_rootRelaxCPEPB" "bc_rootRelaxCP" "bc_rootRelaxEPB" "bc_rootRelax" "bb_EPB" "bc_EPB" "bb" "bc" "epsilon")

# methodes=("bc_rootRelax")

for file in ./instances/*; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        julia vOptMaxCut.jl "$file" $mthd
    done
done
