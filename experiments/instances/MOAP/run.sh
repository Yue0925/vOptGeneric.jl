
# methodes=("bb" "bc")
#  "bb_EPB" "bc_EPB" "bc_rootRelax" "bc_rootRelaxEPB" "bc_rootRelaxCP" "bc_rootRelaxCPEPB" 


# for file in ./AP/*.dat; do
    # for mthd in ${methodes[@]}; do
    #     echo "$file" " ... " $mthd
    #     julia vOptBOAP.jl "$file" $mthd
    # done
# done



# methodes=("bc_rootRelax" "bc_rootRelaxEPB")


# for file in ./AP/*.dat; do
#     for mthd in ${methodes[@]}; do
#         echo "$file" " ... " $mthd
#         julia vOptBOAP2.jl "$file" $mthd
#     done
# done

methodes=("bb" "bb_EPB" "bc_rootRelax" "bc_rootRelaxEPB")
files=("./AP/AP_p-3_n-5_ins-1.dat" "./AP/AP_p-3_n-5_ins-10.dat" "./AP/AP_p-3_n-10_ins-1.dat" "./AP/AP_p-3_n-10_ins-10.dat" "./AP/AP_p-3_n-10_ins-2.dat" "./AP/AP_p-3_n-15_ins-1.dat" "./AP/AP_p-3_n-15_ins-10.dat" "./AP/AP_p-3_n-15_ins-2.dat" "./AP/AP_p-3_n-20_ins-1.dat" "./AP/AP_p-3_n-20_ins-2.dat" "./AP/AP_p-3_n-20_ins-3.dat")

for file in ${files[@]}; do
    for mthd in ${methodes[@]}; do
        echo $file " ... " $mthd
        julia vOptBOAPheur.jl $file $mthd
    done
done

