#!/bin/bash

for file in ./hardinstances_pisinger/*.csv; do
    echo "$file"
    julia vOptHKP2.jl "$file"
done
