
function comparisons(instances::String)
    work_dir = "../../results/" * instances
    @assert isdir(work_dir) "This directory doesn't exist $work_dir !"

    fout = open(work_dir * "/comparisonTable.tex", "w")

    latex = raw"""\begin{table}[!ht]
    \centering
    \resizebox{\columnwidth}{!}{%
    \begin{tabular}{lcccccccccccc}
    \toprule
    \textbf{Instance} & \textbf{n} & \textbf{m} & \multicolumn{2}{c}{\textbf{$\mathbf{\epsilon}$-constraint}}  & \multicolumn{2}{c}{\textbf{B\&B}} & \multicolumn{2}{c}{\textbf{EPB B\&B}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} \\
    
    \cmidrule(r){4-5} \cmidrule(r){6-7} \cmidrule(r){8-9} \cmidrule(r){10-11} \cmidrule(r){12-13} 
    
    ~ & ~ & ~ & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} & \textbf{Time(s)} & \textbf{Node} & \textbf{Time(s)} & \textbf{Node} & \textbf{Time(s)} & \textbf{Node} & \textbf{Time(s)} & \textbf{Node} \\
    \midrule
    """
    println(fout, latex)
    methods = ["epsilon", "bb", "bb_EPB", "bc_rootRelax", "bc_rootRelaxEPB"] 


    n = 0
    count = 0
    avg_n = 0 ; avg_m = 0
    avgT = Dict(k => 0.0 for k in methods) ; avgY = Dict(k => 0.0 for k in methods) ; avgTO = Dict(k => 0 for k in methods)


    # ∀ file in dicho
    for file in readdir(work_dir * "/epsilon/")
        if split(file, ".")[end] == "png"
            continue
        end

        times = [] ; pts = []
        include(work_dir * "/epsilon/" * file)


        # new n folder 
        if n!= vars
            if n!= 0  
                avg_n = round(avg_n/count, digits = 2) ; avg_m = round(avg_m/count, digits = 2)
                for m in methods
                    avgT[m] = round(avgT[m]/count, digits = 2); avgY[m] = round(avgY[m]/count, digits = 2) 
                end

                print(fout, "\\hline avg & " * string(avg_n) * " & " * string(avg_m) )

                for m in methods
                    print(fout, " & " * string(avgT[m]) * "& " * string(avgY[m]))
                end

                println(fout, "\\\\ \\hline")
                println("n = $n , count = $count TO = $avgTO")

            end

            n = vars 
            count = 0
            avg_n = 0 ; avg_m = 0
            avgT = Dict(k => 0.0 for k in methods) ; avgY = Dict(k => 0.0 for k in methods) ; avgTO = Dict(k => 0 for k in methods) 

            # count += 1
            # avg_n += vars ; avg_m += constr
        end

        count += 1
        avg_n += vars ; avg_m += constr

        name_seg = split(file, "_")
        for s in name_seg[1:end-1]
            print(fout, s * "\\_")
        end
        print(fout, name_seg[end] * " & ")
        print(fout, string(vars) * " & " * string(constr) * " & ")
        

        # ∀ method 
        m = methods[1]
        if isfile(work_dir * "/" * m * "/" * file)
            include(work_dir * "/" * m * "/" * file)
            push!(times, total_times_used); push!(pts, size_Y_N)
            total_times_used > 3600.0 ? avgTO[m] += 1 : nothing

            avgT[m] += total_times_used ; avgY[m] += size_Y_N
        else
            push!(times, -1); push!(pts, -1)
        end
        for m in methods[2:end]
            if isfile(work_dir * "/" * m * "/" * file)
                include(work_dir * "/" * m * "/" * file)
                push!(times, total_times_used); push!(pts, total_nodes)
                total_times_used > 3600.0 ? avgTO[m] += 1 : nothing

                avgT[m] += total_times_used ; avgY[m] += total_nodes
            else
                push!(times, -1); push!(pts, -1)
            end
        end

        # ------------------
        for i=1:length(methods)-1
            if times[i] == -1
                print(fout, " - & ")
            else
                print(fout, string(times[i]) * " & ")
            end

            if pts[i] == -1
                print(fout, " - & ")
            else
                print(fout, string(pts[i]) * " & ")
            end

        end

        print(fout, string(times[end]) * " & ") 

        println(fout, string(pts[end]) * " \\\\")

    end

    avg_n = round(avg_n/count, digits = 2) ; avg_m = round(avg_m/count, digits = 2)
    for m in methods
        avgT[m] = round(avgT[m]/count, digits = 2); avgY[m] = round(avgY[m]/count, digits = 2) 
    end

    print(fout, "\\hline avg & " * string(avg_n) * " & " * string(avg_m) )

    for m in methods
        print(fout, " &  " * string(avgT[m]) * " & " * string(avgY[m]))
    end

    println(fout, " \\\\ \\hline")
    println("n = $n , count = $count TO = $avgTO")


    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout, latex)
    println(fout, "\\caption{Comparison of the different algorithms performances for instances $instances .}")
    println(fout, "%\\label{tab:table_compare_$instances }")
    println(fout, "\\end{table}")
    close(fout)

end


function Cut_Branch(instances::String)
    work_dir = "../../results/" * instances
    @assert isdir(work_dir * "/mixed2") "This directory doesn't exist $work_dir !"
 
    methods = ["bc_rootRelax"] 

    n = 0
    count = 0
    avg_n = 0 ; avg_m = 0
    avgT = Dict(k => 0.0 for k in methods) ; avgY = Dict(k => 0.0 for k in methods) ; avgTO = Dict(k => 0 for k in methods)


    # ∀ file in dicho
    for file in readdir(work_dir * "/epsilon/")
        if split(file, ".")[end] == "png"
            continue
        end

        times = [] ; pts = []
        include(work_dir * "/epsilon/" * file)

        # new n folder 
        if n!= vars
            if n!= 0  
                avg_n = round(avg_n/count, digits = 2) ; avg_m = round(avg_m/count, digits = 2)
                for m in methods
                    avgT[m] = round(avgT[m]/count, digits = 2); avgY[m] = round(avgY[m]/count, digits = 2) 
                end

                for m in methods
                    println("n = $n , count = $count TO = $avgTO time $(avgT[m]) , node $(avgY[m])")
                end

            end

            n = vars 
            count = 0
            avg_n = 0 ; avg_m = 0
            avgT = Dict(k => 0.0 for k in methods) ; avgY = Dict(k => 0.0 for k in methods) ; avgTO = Dict(k => 0 for k in methods) 
        end

        count += 1
        avg_n += vars ; avg_m += constr

        for m in methods
            if isfile(work_dir * "/mixed2"* "/" * m * "/" * file)
                include(work_dir * "/mixed2" * "/" * m * "/" * file)
                push!(times, total_times_used); push!(pts, total_nodes)
                total_times_used > 3600.0 ? avgTO[m] += 1 : nothing

                avgT[m] += total_times_used ; avgY[m] += total_nodes
            else
                push!(times, -1); push!(pts, -1)
            end
        end

    end

    avg_n = round(avg_n/count, digits = 2) ; avg_m = round(avg_m/count, digits = 2)
    for m in methods
        avgT[m] = round(avgT[m]/count, digits = 2); avgY[m] = round(avgY[m]/count, digits = 2) 
    end


    for m in methods
        println("n = $n , count = $count TO = $avgTO time $(avgT[m]) , node $(avgY[m])")
    end

end

# comparisons("MOAP/AP")
Cut_Branch("MOAP/AP")
# n = 100 , count = 10 TO = Dict("bc_rootRelax" => 0, "bb" => 0, "epsilon" => 0, "bc_rootRelaxEPB" => 0, "bb_EPB" => 0)
# n = 225 , count = 10 TO = Dict("bc_rootRelax" => 0, "bb" => 0, "epsilon" => 0, "bc_rootRelaxEPB" => 0, "bb_EPB" => 0)
# n = 400 , count = 10 TO = Dict("bc_rootRelax" => 7, "bb" => 0, "epsilon" => 0, "bc_rootRelaxEPB" => 0, "bb_EPB" => 6)
# n = 625 , count = 10 TO = Dict("bc_rootRelax" => 10, "bb" => 9, "epsilon" => 0, "bc_rootRelaxEPB" => 5, "bb_EPB" => 10)
# n = 900 , count = 10 TO = Dict("bc_rootRelax" => 10, "bb" => 10, "epsilon" => 0, "bc_rootRelaxEPB" => 10, "bb_EPB" => 10)
# n = 1225 , count = 10 TO = Dict("bc_rootRelax" => 10, "bb" => 10, "epsilon" => 0, "bc_rootRelaxEPB" => 10, "bb_EPB" => 10)
# n = 1600 , count = 10 TO = Dict("bc_rootRelax" => 10, "bb" => 10, "epsilon" => 0, "bc_rootRelaxEPB" => 10, "bb_EPB" => 10)
# n = 2025 , count = 10 TO = Dict("bc_rootRelax" => 10, "bb" => 10, "epsilon" => 0, "bc_rootRelaxEPB" => 10, "bb_EPB" => 10)
# n = 2500 , count = 10 TO = Dict("bc_rootRelax" => 2, "bb" => 10, "epsilon" => 0, "bc_rootRelaxEPB" => 1, "bb_EPB" => 10)
# n = 25 , count = 10 TO = Dict("bc_rootRelax" => 0, "bb" => 0, "epsilon" => 0, "bc_rootRelaxEPB" => 0, "bb_EPB" => 0)
