#!/bin/bash


for file in ./ORLib/*.txt; do
    echo "$file"
    julia vOptSCP.jl "$file"
done
