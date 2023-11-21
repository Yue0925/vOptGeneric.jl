methodes=("bc_rootRelax" "bc_rootRelaxEPB")

for file in ./instances/*; do
    for mthd in ${methodes[@]}; do
        echo "$file ... " $mthd
        julia vOptSCP2.jl "$file" $mthd
    done
done
