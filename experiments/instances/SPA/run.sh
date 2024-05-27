#!/bin/bash

# for file in ./BOSPA/*.txt; do
#     echo "$file"
#     julia vOptBOSPA.jl "$file"
# done

# julia latexWriter.jl



methodes=("bb" "bc_rootRelax" "bc_rootRelaxEPB")
# "bc" "bb_EPB" "bc_EPB" "bc_rootRelax" "bc_rootRelaxEPB" "bc_rootRelaxCP" "bc_rootRelaxCPEPB" 

for mthd in ${methodes[@]}; do
    echo " ... " $mthd
    julia vOptBOSPA.jl ./BOSPA/biosppnw12.txt $mthd
done
