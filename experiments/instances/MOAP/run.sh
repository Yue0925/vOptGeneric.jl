
methodes=("bb" "bc")
#  "bb_EPB" "bc_EPB" "bc_rootRelax" "bc_rootRelaxEPB" "bc_rootRelaxCP" "bc_rootRelaxCPEPB" 


for file in ./AP/*.dat; do
    for mthd in ${methodes[@]}; do
        echo "$file" " ... " $mthd
        julia vOptBOAP.jl "$file" $mthd
    done
done
