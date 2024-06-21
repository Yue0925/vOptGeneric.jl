methodes=("bc_rootRelax")
#  "bb_EPB" "bc_EPB" "bc_rootRelax" "bc_rootRelaxEPB" "bc_rootRelaxCP" "bc_rootRelaxCPEPB" 


for file in ./instance/*.raw; do
    for mthd in ${methodes[@]}; do
        echo "$file" " ... " $mthd
        julia vOptUFLP_Forget20mixed.jl "$file" $mthd
    done
done

