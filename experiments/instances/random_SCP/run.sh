#!/bin/bash


for file in ./instances/*; do
    echo "$file"
    julia vOptSCP.jl "$file"
done
