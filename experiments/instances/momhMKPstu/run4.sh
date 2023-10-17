#!/bin/bash

for file in ./MOBKP/set3/*; do
    echo "$file"
    julia vOptMomkp3.jl "$file"
done

