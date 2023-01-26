#!/bin/bash

for file in ./hardinstances_pisinger/*.csv; do
    echo "$file"
    julia vOptHKP.jl "$file"
done
