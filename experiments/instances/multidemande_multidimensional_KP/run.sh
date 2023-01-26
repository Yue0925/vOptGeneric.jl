#!/bin/bash

for file in ./ORLib/*.txt; do
    echo "$file"
    julia vOptMDMDKP.jl "$file"
done
