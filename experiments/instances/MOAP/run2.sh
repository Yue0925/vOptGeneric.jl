


methodes=("bb" "bc_rootRelax" "bc_rootRelaxEPB")
files=("./AP/AP_p-3_n-10_ins-10.dat" "./AP/AP_p-3_n-10_ins-2.dat" "./AP/AP_p-3_n-15_ins-1.dat" "./AP/AP_p-3_n-15_ins-10.dat" "./AP/AP_p-3_n-15_ins-2.dat" "./AP/AP_p-3_n-20_ins-1.dat" "./AP/AP_p-3_n-20_ins-2.dat" "./AP/AP_p-3_n-20_ins-3.dat")

# "./AP/AP_p-3_n-5_ins-1.dat" "./AP/AP_p-3_n-5_ins-10.dat" "./AP/AP_p-3_n-10_ins-1.dat"

for file in ${files[@]}; do
    for mthd in ${methodes[@]}; do
        echo $file " ... " $mthd
        julia vOptBOAPmixedCopy.jl $file $mthd
    done
done

