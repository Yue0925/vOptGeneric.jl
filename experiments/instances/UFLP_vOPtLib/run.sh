#!/bin/bash


for file in ./dataDidactic/*.txt; do
    echo "$file"
    julia uncapacitatedFacilityLocation.jl "$file"
done

for file in ./dataFernandez/*.txt; do
    echo "$file"
    julia uncapacitatedFacilityLocation.jl "$file"
done

# julia latexWriter.jl