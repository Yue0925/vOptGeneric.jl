#!/bin/bash

for file in ./ALL_tsp/*.gz; do
    echo "$file"
    gzip -d $file
done