
# using PyPlot
# import PyPlot; const plt = PyPlot


function comparisons(instances::String)
    work_dir = "../../results/" * instances
    @assert isdir(work_dir) "This directory doesn't exist $work_dir !"

    fout = open(work_dir * "/comparisonTable.tex", "w")

    latex = raw"""\begin{sidewaystable}[!ht]
    \centering
    \resizebox{\columnwidth}{!}{%
    \begin{tabular}{lcccccccccccccccccccc}
    \toprule
    \textbf{Instance} & \textbf{n} & \textbf{m} & \multicolumn{2}{c}{\textbf{$\mathbf{\epsilon}$-constraint}}  & \multicolumn{2}{c}{\textbf{B\&B}} & \multicolumn{2}{c}{\textbf{B\&C(LP+CP)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&B}} & \multicolumn{2}{c}{\textbf{EPB B\&C(LP+CP)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex+CP)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex+CP)}} \\
    
    \cmidrule(r){4-5} \cmidrule(r){6-7} \cmidrule(r){8-9} \cmidrule(r){10-11} \cmidrule(r){12-13} \cmidrule(r){14-15} \cmidrule(r){16-17} \cmidrule(r){18-19} \cmidrule(r){20-21}
    ~ & ~ & ~ & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} \\
    \midrule
    """
    println(fout, latex)
    methods = ["epsilon", "bb", "bc", "bc_rootRelax", "bb_EPB", "bc_EPB", "bc_rootRelaxEPB", "bc_rootRelaxCP", "bc_rootRelaxCPEPB"]

    # ∀ file in dicho
    for file in readdir(work_dir * "/bc_rootRelaxCPEPB/")
        if split(file, ".")[end] == "png"
            continue
        end

        name_seg = split(file, "_")
        for s in name_seg[1:end-1]
            print(fout, s * "\\_")
        end
        print(fout, name_seg[end] * " & ")
        times = [] ; pts = []

        include(work_dir * "/bc_rootRelaxCPEPB/" * file)
        print(fout, string(vars) * " & " * string(constr) * " & ")

        # ∀ method 
        for m in methods
            if isfile(work_dir * "/" * m * "/" * file)
                include(work_dir * "/" * m * "/" * file)
                push!(times, total_times_used); push!(pts, size_Y_N)

            else
                push!(times, -1); push!(pts, -1)
            end
        end

        # ------------------
        for i=1:length(methods)-1
            if times[i] == -1
                print(fout, " - & ")                
            elseif times[i] == minimum(filter(x -> x > 0 ,times))
                times[i] >= 3600.0 ? print(fout, "TO & ") : print(fout, " \\textcolor{blue2}{" * string(times[i]) * "} & ")
            else 
                times[i] >= 3600.0 ? print(fout, "TO & ") : print(fout, string(times[i]) * " & ")
            end

            if pts[i] == -1
                print(fout, " - & ")
            else
                print(fout, string(pts[i]) * " & ")
            end

        end

        if times[end] == minimum(filter(x -> x > 0 ,times))
            times[end] >= 3600.0 ? print(fout, "TO & ") : print(fout, " \\textcolor{blue2}{" * string(times[end]) * " & ")
        else
            times[end] >= 3600.0 ? print(fout, "TO & ") : print(fout, string(times[end]) * " & ")
        end
        println(fout, string(pts[end]) * " \\\\")

    end

    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout, latex)
    println(fout, "\\caption{Comparison of the different algorithms performances for instances $instances .}")
    println(fout, "\\label{tab:table_compare_$instances }")
    println(fout, "\\end{sidewaystable}")
    close(fout)

end


comparisons("MDMDKP")