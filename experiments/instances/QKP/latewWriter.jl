
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
    for file in readdir(work_dir * "/bc_rootRelax/")
        if split(file, ".")[end] == "png" || split(file, ".")[end] == "tex"
            continue
        end

        name_seg = split(file, "_")
        for s in name_seg[1:end-1]
            print(fout, s * "\\_")
        end
        print(fout, name_seg[end] * " & ")
        times = [] ; pts = []

        # write dichotomy result 
        include(work_dir * "/bc_rootRelax/" * file)
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



function MOBB_perform(instances::String)
    methods = ["bb", "bc", "bc_rootRelax", "bb_EPB", "bc_EPB", "bc_rootRelaxEPB", "bc_rootRelaxCP", "bc_rootRelaxCPEPB"]

    for m in methods
        dir = "../../results/" * instances * "/" * m
        if !isdir(dir) 
            @warn "This directory doesn't exist $dir !"
            continue
        end

        fout = open(dir * "/" * m * ".tex", "w")

        latex = raw"""\begin{table}[!h]
        \centering
        \resizebox{\columnwidth}{!}{%
        \hspace*{-1cm}\begin{tabular}{lccccccccccc}
        \toprule
        \textbf{Instance} & \textbf{n} & \textbf{m} & \multicolumn{4}{c}{\textbf{Time(s)}} & \multicolumn{2}{c}{\textbf{Nodes}}  & \textbf{Tree(MB)} & \textbf{$|\mathcal{Y}_N|$} & \textbf{$|\mathcal{X}_E|$}
        \\
        \cmidrule(r){4-7} \cmidrule(r){8-9} 
        ~ & ~ & ~ & \textbf{total} &\textbf{relax} & \textbf{dominance} & \textbf{incumbent} & \textbf{total} & \textbf{pruned} & ~ & ~ & ~ \\
        \midrule
        """
        println(fout, latex)


        for file in readdir(dir * "/")
            if split(file, ".")[end] == "png" || split(file, ".")[end] == "tex"
                continue
            end

            name_seg = split(file, "_")
            for s in name_seg[1:end-1]
                print(fout, s * "\\_")
            end
            print(fout, name_seg[end] * " & ")

            include(dir * "/" * file)
            print(fout, string(vars) * " & " * string(constr) * " & ")

            print(fout, string(total_times_used) * " & " * string(relaxation_time) * " & " * string(test_dominance_time) * " & "*
                string( update_incumbent_time ) * " & " * string(total_nodes) * " & " * string(pruned_nodes)*
                " & " * string(tree_size) * " & " * string(size_Y_N) * " & " * string(size_X_E)
            )

            println(fout, "\\\\")

        end

        latex = raw"""\bottomrule
        \end{tabular}%
        }%
        \caption{The detailed experimental information about BB algorithm.}
        \label{tab:table_$m}
        \end{table}
        """
        println(fout, latex)
        close(fout)
    end
end



# MOBB_perform("QKP")
comparisons("QKP")