#!/bin/bash


for file in ./instances/PCC/*; do
    echo "$file"
    julia vOptMDMKP.jl "$file" epsilon
done
