#!/bin/bash

# for file in ./mipLib/*.mps; do
#     echo "$file"
#     julia parserBP.jl "$file"
# done



methodes=("bb" "bc" "bb_EPB" "bc_EPB")

for file in ./mipLib/*.mps; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        srun julia vOptBP.jl "$file" $mthd
    done
done



for file in ./mipLib/*.mps; do
        echo "$file ... " 
        srun julia vOptBP.jl "$file" epsilon
done
