#!/bin/bash


for subfolder in ./MOBKP/*; do
    for file in $subfolder/*; do
        echo "$file"
        julia vOptMomkp.jl "$file"
    done
done


for file in ./MOMKP/*; do
    echo "$file"
    julia vOptMomkp.jl "$file"
done