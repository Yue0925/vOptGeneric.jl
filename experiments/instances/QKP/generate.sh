#!/bin/bash

g++ -o kqkp kqkp_inst.cc

density=(25 50 75 100)

for ((i=5;i<=25;i++)); do
    for d in ${density[@]}; do
        ./kqkp $i 100 $d 333
    done
done