#!/bin/bash


# for file in ./instances/*; do
#     echo "$file"
#     julia vOptMaxCut.jl "$file" epsilon
# done

methodes=("bc_rootRelax" "bc_rootRelaxCPEPB" "bc_rootRelaxCP" "bc_rootRelaxEPB")

# methodes=("bb" "bc_rootRelax" "bb_EPB" "bc_EPB" "bc")

for file in ./instances/*; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        julia vOptMaxCut.jl "$file" $mthd
    done
done
