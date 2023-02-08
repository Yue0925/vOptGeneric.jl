#!/bin/bash

# for file in ./mipLib/*.mps; do
#     echo "$file"
#     julia parserBP.jl "$file"
# done


for file in ./mipLib/*.mps; do
    echo "$file ... " 
    julia vOptBP.jl "$file" epsilon
done
