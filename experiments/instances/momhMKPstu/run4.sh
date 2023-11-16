#!/bin/bash

for file in ./MOBKP/set3/*; do
    echo "$file"
    julia vOptMomkp4.jl "$file"
done

