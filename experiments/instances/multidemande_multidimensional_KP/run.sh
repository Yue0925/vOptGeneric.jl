#!/bin/bash

# for file in ./ORLib/*.txt; do
#     echo "$file"
#     julia vOptMDMDKP.jl "$file"
# done



for file in ./ORLib/*.txt; do
    echo "$file ... "
    julia vOptMDMDKP_alt.jl "$file" epsilon
done
