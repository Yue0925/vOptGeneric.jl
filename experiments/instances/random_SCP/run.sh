#!/bin/bash


methodes=("epsilon" "bb" "bc" "bb_EPB" "bc_EPB")

for file in ./instances/*; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        julia vOptSCP.jl "$file" $mthd
    done
done