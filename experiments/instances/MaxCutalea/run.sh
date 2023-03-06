#!/bin/bash


for file in ./instances/*; do
    echo "$file"
    julia vOptMaxCut.jl "$file" epsilon
done
