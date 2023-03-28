#!/bin/bash

density=(25 50 100)

for ((i=5;i<=20;i++)); do
    for d in ${density[@]}; do
        ./kqkp $i 100 $d 
    done
done