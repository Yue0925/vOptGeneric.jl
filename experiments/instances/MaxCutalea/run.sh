#!/bin/bash


# for file in ./instances/*; do
#     echo "$file"
#     julia vOptMaxCut.jl "$file" epsilon
# done

methodes=("bb" "bc_rootRelax" "bc_rootRelaxCPEPB" "bc_rootRelaxCP" "bc_rootRelaxEPB" "bb_EPB" "bc_EPB" "bc")

# methodes=("bb" "bc_rootRelax")

for file in ./instances/*; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        julia vOptMaxCut.jl "$file" $mthd
    done
done
