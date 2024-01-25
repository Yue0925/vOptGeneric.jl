#!/bin/bash


# for file in ./instances/*; do
#     echo "$file"
#     julia vOptQKP.jl "$file" epsilon
# done

# methodes=("bc_rootRelaxCPEPB" "bc_rootRelaxCP" "bc_rootRelaxEPB" "bb_EPB" "bc_EPB")

# methodes=("bc_rootRelax" "bc_rootRelaxEPB" "bc_rootRelaxCP" "bc_rootRelaxCPEPB")
# # "bc_rootRelax" "bc"

# for file in ./instances/*; do
#     for mthd in ${methodes[@]}; do
#         echo "$file ... " $mthd
#         julia vOptQKP.jl "$file" $mthd
#     done
# done


# --------- porta 1 write .ieq file
# for file in ./instances/*; do
#     echo "$file ... "
#     julia plotObjectiveSpace.jl "$file"
# done

# -------- porta 2 generate feasible point  
# for file in ./porta/*.ieq; do
#     echo "$file ... "
#     vint "$file"
# done

# -------- porta 3 image feasible set Y
for file in ./instances/*; do
    echo "$file ... "
    julia plotObjectiveSpace.jl "$file"
done

