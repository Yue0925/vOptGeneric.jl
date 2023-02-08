#!/bin/bash

# for file in ./ORLib/*.txt; do
#     echo "$file"
#     julia vOptMDKP.jl "$file"
# done




for file in ./ORLib/*.txt; do
    echo "$file ... "
    julia vOptMDKP_alt.jl "$file" epsilon
done
