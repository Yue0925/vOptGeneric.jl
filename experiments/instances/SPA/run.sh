#!/bin/bash

# for file in ./BOSPA/*.txt; do
#     echo "$file"
#     julia vOptBOSPA.jl "$file"
# done

# julia latexWriter.jl



# methodes=("bc_rootRelaxEPB")
# # "bc" "bb_EPB" "bc_EPB" "bc_rootRelax" "bc_rootRelaxEPB" "bc_rootRelaxCP" "bc_rootRelaxCPEPB" 

# for mthd in ${methodes[@]}; do
#     echo " ... " $mthd
#     julia vOptBOSPA.jl ./BOSPA/biosppnw10.txt $mthd
# done


# methodes=("bc_rootRelax" "bc_rootRelaxEPB")


# for file in ./BOSPA/*.txt; do
#     for mthd in ${methodes[@]}; do
#         echo "$file" " ... " $mthd
#         julia vOptBOSPA2.jl "$file" $mthd
#     done
# done


methodes=("bc_rootRelaxEPB") # "bb" "bc_rootRelax" 
# files=("./BOSPA/biosppnw24.txt")

# for file in ${files[@]}; do
    for mthd in ${methodes[@]}; do
        echo " ./BOSPA/biosppnw24.txt ... " $mthd
        julia vOptBOSPA.jl ./BOSPA/biosppnw24.txt $mthd
    done
# done

