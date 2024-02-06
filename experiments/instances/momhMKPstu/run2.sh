#!/bin/bash

# for file in ./MOBKP/set2/*; do
#     echo "$file"
#     julia vOptMomkp2.jl "$file"
# done


for file in ./MOBKP/set3/*; do
    echo "$file"
    julia vOptMomkp2.jl "$file"
done