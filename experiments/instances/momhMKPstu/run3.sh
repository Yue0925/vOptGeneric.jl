#!/bin/bash



for file in ./MOBKP/set3/*; do
    echo "$file"
    julia vOptMomkp_complementaire.jl "$file"
done
