


for file in ./KP/*.raw; do
    # for mthd in ${methodes[@]}; do
        echo "$file" " ... "
        # julia vOptUFLP_Forget20.jl "$file" $mthd
        julia parserKP_Forget20.jl "$file"
    # done
done

