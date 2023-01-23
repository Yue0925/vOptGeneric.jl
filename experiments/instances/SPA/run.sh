#!/bin/bash

for file in ./BOSPA/*.txt; do
    echo "$file"
    julia vOptBOSPA.jl "$file"
done

# julia latexWriter.jl