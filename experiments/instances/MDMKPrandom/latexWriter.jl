
# using PyPlot
# import PyPlot; const plt = PyPlot



"""
Comparison tree table time & nodes explored (all B&C methodes)
"""
function comparisons5(instances::String)
    work_dir = "../../results/" * instances
    @assert isdir(work_dir) "This directory doesn't exist $work_dir !"

    fout = open(work_dir * "/comparisonBCcplexTable.tex", "w")

    latex = raw"""\begin{table}[!ht]
    \centering
    \resizebox{\columnwidth}{!}{%
    \hspace*{-1cm}\begin{tabular}{lcccccccccccccc}
    \toprule
    \textbf{Instance} & \textbf{n} & \textbf{m} & \multicolumn{2}{c}{\textbf{B\&C}} & \multicolumn{2}{c}{\textbf{B\&C (cplex)}} & \multicolumn{2}{c}{\textbf{B\&C (cplex + CP)}}  & \multicolumn{2}{c}{\textbf{EPB B\&C}} & \multicolumn{2}{c}{\textbf{EPB B\&C (cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C (cplex + CP)}} \\
    
    \cmidrule(r){4-5} \cmidrule(r){6-7} \cmidrule(r){8-9} \cmidrule(r){10-11} \cmidrule(r){12-13} \cmidrule(r){14-15}  
    ~ & ~ & ~ & \textbf{Time(s)} &\textbf{Nodes} & \textbf{Time(s)} &\textbf{Nodes} & \textbf{Time(s)} &\textbf{Nodes} & \textbf{Time(s)} &\textbf{Nodes} & \textbf{Time(s)} &\textbf{Nodes} & \textbf{Time(s)} &\textbf{Nodes}  \\
    \midrule
    """
    println(fout, latex)
    methods = ["bc", "bc_rootRelax", "bc_rootRelaxCP", "bc_EPB", "bc_rootRelaxEPB", "bc_rootRelaxCPEPB"] 


    n = 0
    count = 0
    avg_n = 0 ; avg_m = 0
    avgT = Dict(k => 0.0 for k in methods) ; avgY = Dict(k => 0.0 for k in methods)


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
            end

            n = vars 
            count = 0
            avg_n = 0 ; avg_m = 0
            avgT = Dict(k => 0.0 for k in methods) ; avgY = Dict(k => 0.0 for k in methods)

            count += 1
            avg_n += vars ; avg_m += constr
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
        for m in methods
            if isfile(work_dir * "/" * m * "/" * file)
                include(work_dir * "/" * m * "/" * file)
                push!(times, total_times_used); push!(pts, total_nodes)

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

    avg_n = round(Int, avg_n/count) ; avg_m = round(Int, avg_m/count)
    for m in methods
        avgT[m] = round(avgT[m]/count, digits = 2); avgY[m] = round(avgY[m]/count, digits = 2) 
    end

    print(fout, "\\hline avg & " * string(avg_n) * " & " * string(avg_m) )

    for m in methods
        print(fout, " &  " * string(avgT[m]) * " & " * string(avgY[m]))
    end

    println(fout, " \\\\ \\hline")


    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout, latex)
    println(fout, "\\caption{Comparison of the different algorithms performances for instances $instances .}")
    println(fout, "%\\label{tab:table_compare_$instances }")
    println(fout, "\\end{sidewaystable}")
    close(fout)

end



# compare B&B/B&C vs EPB B&B B&C 
function comparisons4(instances::String)
    work_dir = "../../results/" * instances
    @assert isdir(work_dir) "This directory doesn't exist $work_dir !"

    fout = open(work_dir * "/comparisonBBBCTable.tex", "w")

    latex = raw"""\begin{table}[!ht]
    \centering
    % \resizebox{\columnwidth}{!}{%
    \hspace*{-1cm}\begin{tabular}{lcccccccccc}
    \toprule
    \textbf{Instance} & \textbf{n} & \textbf{m} & \multicolumn{2}{c}{\textbf{B\&B}} & \multicolumn{2}{c}{\textbf{B\&C}}  & \multicolumn{2}{c}{\textbf{EPB B\&B}} & \multicolumn{2}{c}{\textbf{EPB B\&C}} \\
    
    \cmidrule(r){4-5} \cmidrule(r){6-7} \cmidrule(r){8-9} \cmidrule(r){10-11} 
    ~ & ~ & ~ & \textbf{Time(s)} &\textbf{Nodes} & \textbf{Time(s)} &\textbf{Nodes} & \textbf{Time(s)} &\textbf{Nodes} & \textbf{Time(s)} &\textbf{Nodes}  \\
    \midrule
    """
    println(fout, latex)
    methods = ["bb", "bc", "bb_EPB", "bc_EPB"] 

    n = 0
    count = 0
    avg_n = 0 ; avg_m = 0
    avgT = Dict(k => 0.0 for k in methods) ; avgY = Dict(k => 0.0 for k in methods)

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
            end

            n = vars 
            count = 0
            avg_n = 0 ; avg_m = 0
            avgT = Dict(k => 0.0 for k in methods) ; avgY = Dict(k => 0.0 for k in methods)

            count += 1
            avg_n += vars ; avg_m += constr
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
        for m in methods
            if isfile(work_dir * "/" * m * "/" * file)
                include(work_dir * "/" * m * "/" * file)
                push!(times, total_times_used); push!(pts, total_nodes)

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

    avg_n = round(avg_n/count, digits = 2) ; avg_m = round( avg_m/count, digits = 2)
    for m in methods
        avgT[m] = round(avgT[m]/count, digits = 2); avgY[m] = round(avgY[m]/count, digits = 2) 
    end

    print(fout, "\\hline avg & " * string(avg_n) * " & " * string(avg_m) )

    for m in methods
        print(fout, " &  " * string(avgT[m]) * " & " * string(avgY[m]))
    end

    println(fout, " \\\\ \\hline")

    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout, latex)
    println(fout, "\\caption{Comparison of the different algorithms performances for instances $instances .}")
    println(fout, "%\\label{tab:table_compare_$instances }")
    println(fout, "\\end{sidewaystable}")
    close(fout)

end


# compare SP and MP cuts 
function comparisonsCP(instances::String)
    work_dir = "../../results/" * instances
    @assert isdir(work_dir) "This directory doesn't exist $work_dir !"

    fout = open(work_dir * "/comparisonBCTable.tex", "w")

    latex = raw"""\begin{table}[!ht]
    \centering
    % \resizebox{\columnwidth}{!}{%
    \hspace*{-1cm}\begin{tabular}{lcccccc}
    \toprule
    \textbf{Instance} & \textbf{n} & \textbf{m} & \multicolumn{2}{c}{\textbf{B\&C}}  & \multicolumn{2}{c}{\textbf{EPB B\&C}}
    \\
    \cmidrule(r){4-5} \cmidrule(r){6-7} 
    ~ & ~ & ~ & \textbf{SP cuts} &\textbf{MP cuts} & \textbf{SP cuts} &\textbf{MP cuts} \\
    \midrule
    """
    println(fout, latex)
    methods = ["bc", "bc_EPB"] ;  

    n = 0
    count = 0
    avg_n = 0 ; avg_m = 0
    avgSP = Dict(k => 0.0 for k in methods) ; avgMP = Dict(k => 0.0 for k in methods)


    # ∀ file in dicho
    for file in readdir(work_dir * "/epsilon/")
        if split(file, ".")[end] == "png"
            continue
        end

        sp = [] ; mp = []

        include(work_dir * "/epsilon/" * file)


        # new n folder 
        if n!= vars
            if n!= 0  
                avg_n = round( avg_n/count, digits = 2) ; avg_m = round(avg_m/count, digits = 2)
                for m in methods
                    avgSP[m] = round(avgSP[m]/count, digits = 2); avgMP[m] = round(avgMP[m]/count, digits = 2) 
                end

                print(fout, "\\cline{1-9} avg & " * string(avg_n) * " & " * string(avg_m) )

                for m in methods
                    print(fout, " & " * string(avgSP[m]) * "& " * string(avgMP[m]))
                end

                println(fout, "\\\\ \\hline")
            end

            n = vars 
            count = 0
            avg_n = 0 ; avg_m = 0
            avgSP = Dict(k => 0.0 for k in methods) ; avgMP = Dict(k => 0.0 for k in methods)

            count += 1
            avg_n += vars ; avg_m += constr
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
        for m in methods
            if isfile(work_dir * "/" * m * "/" * file)
                include(work_dir * "/" * m * "/" * file)
                push!(sp, sp_cuts); push!(mp, mp_cuts)

                avgSP[m] += sp_cuts ; avgMP[m] += mp_cuts
            else
                push!(sp, -1); push!(mp, -1)
            end
        end

        # ------------------
        for i=1:length(methods)-1
            if sp[i] == -1
                print(fout, " - & ")
            else
                print(fout, string(sp[i]) * " & ")
            end

            if mp[i] == -1
                print(fout, " - & ")
            else
                print(fout, string(mp[i]) * " & ")
            end

        end


        print(fout, string(sp[end]) * " & ") 
        
        println(fout, string(mp[end]) * " \\\\")

    end

    avg_n = round( avg_n/count, digits = 2) ; avg_m = round(avg_m/count, digits = 2)
    for m in methods
        avgSP[m] = round(avgSP[m]/count, digits = 2); avgMP[m] = round(avgMP[m]/count, digits = 2) 
    end

    print(fout, "\\cline{1-7} \\textbf{avg} & \\textbf{" * string(avg_n) * "} & \\textbf{" * string(avg_m) )

    for m in methods
        print(fout, "} & \\textbf{" * string(avgSP[m]) * "} & \\textbf{" * string(avgMP[m]))
    end

    println(fout, "} " * "\\\\ \\cline{1-7}")

    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout, latex)
    println(fout, "\\caption{Comparison of the B\\&B algorithms performances for instances $instances .}")
    println(fout, "%\\label{tab:table_compareBB_$instances }")
    println(fout, "\\end{table}")
    close(fout)

end



function comparisons_eps_BB_EPB(instances::String)
    work_dir = "../../results/" * instances
    @assert isdir(work_dir) "This directory doesn't exist $work_dir !"

    fout = open(work_dir * "/comparisonEpsBBTable.tex", "w")

    latex = raw"""\begin{table}[!ht]
    \centering
    \resizebox{\columnwidth}{!}{%
    \hspace*{-1cm}\begin{tabular}{lcccccccc}
    \toprule
    \textbf{Instance} & \textbf{n} & \textbf{m} & \textbf{$\epsilon$-constraint} & \multicolumn{2}{c}{\textbf{B\&B}}  & \multicolumn{2}{c}{\textbf{EPB B\&B}} & \textbf{$|\mathcal{Y}_N|$}
    \\
    \cmidrule(r){5-6} \cmidrule(r){7-8} 
    ~ & ~ & ~ & ~ & \textbf{Time(s)} &\textbf{Nodes} & \textbf{Time(s)} &\textbf{Nodes} & ~ \\
    \midrule
    """
    println(fout, latex)
    methods = ["epsilon", "bb", "bb_EPB"] 

    n = 0
    count = 0
    avg_n = 0 ; avg_m = 0
    avgT = Dict(k => 0.0 for k in methods) ; avgY = Dict(k => 0.0 for k in methods)
    countY_N = 0


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
                avg_n = round( avg_n/count, digits = 2) ; avg_m = round(avg_m/count, digits = 2)
                for m in methods
                    avgT[m] = round(avgT[m]/count, digits = 2); avgY[m] = round(avgY[m]/count, digits = 2) 
                end

                print(fout, "\\cline{1-9} avg & " * string(avg_n) * " & " * string(avg_m) )

                print(fout, " & " * string(avgT[methods[1]]) )
                for m in methods[2:end]
                    print(fout, " & " * string(avgT[m]) * "& " * string(avgY[m]))
                end

                println(fout, " & " * string(round(countY_N/count, digits = 2))* "\\\\ \\cline{1-9}")
            end

            n = vars 
            count = 0
            avg_n = 0 ; avg_m = 0
            avgT = Dict(k => 0.0 for k in methods) ; avgY = Dict(k => 0.0 for k in methods)
            countY_N = 0

            count += 1
            avg_n += vars ; avg_m += constr
            countY_N += size_Y_N
        end
        count += 1
        avg_n += vars ; avg_m += constr
        countY_N += size_Y_N

        name_seg = split(file, "_")
        for s in name_seg[1:end-1]
            print(fout, s * "\\_")
        end
        print(fout, name_seg[end] * " & ")
        print(fout, string(vars) * " & " * string(constr) * " & ")


        # ∀ method 
        for m in methods
            if isfile(work_dir * "/" * m * "/" * file)
                include(work_dir * "/" * m * "/" * file)
                push!(times, total_times_used); avgT[m] += total_times_used
                
                if m == "epsilon"
                    push!(pts, 0) ; avgY[m] += 0
                else
                    push!(pts, total_nodes) ; avgY[m] += total_nodes
                end
                
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

            if i==1 continue end 
            if pts[i] == -1
                print(fout, " - & ")
            else
                print(fout, string(pts[i]) * " & ")
            end

        end

        if times[end] == minimum(filter(x -> x > 0 ,times))
            times[end] >= 3600.0 ? print(fout, "TO & ") : print(fout, " \\textcolor{blue2}{" * string(times[end]) * "} & ")
        else
            times[end] >= 3600.0 ? print(fout, "TO & ") : print(fout, string(times[end]) * " & ")
        end
        
        println(fout, string(pts[end]) * " & " * string(size_Y_N) * " \\\\")

    end

    avg_n = round( avg_n/count, digits = 2) ; avg_m = round(avg_m/count, digits = 2)
    for m in methods
        avgT[m] = round(avgT[m]/count, digits = 2); avgY[m] = round(avgY[m]/count, digits = 2) 
    end

    print(fout, "\\cline{1-9} avg & " * string(avg_n) * " & " * string(avg_m) )

    print(fout, " & " * string(avgT[methods[1]]) )
    for m in methods[2:end]
        print(fout, " & " * string(avgT[m]) * "& " * string(avgY[m]))
    end

    println(fout, " & " * string(round(countY_N/count, digits = 2))* "\\\\ \\cline{1-9}")


    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout, latex)
    println(fout, "\\caption{Comparison of the B\\&B algorithms performances for instances $instances .}")
    println(fout, "\\label{tab:table_EPSILONvsBBvsEPBBB_$instances }")
    println(fout, "\\end{table}")
    close(fout)
end


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
        if split(file, ".")[end] == "png"
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
            times[end] >= 3600.0 ? print(fout, "TO & ") : print(fout, " \\textcolor{blue2}{" * string(times[end]) * "} & ")
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




# comparisons("MDMKPrandom")
# comparisons_eps_BB_EPB("MDMKPrandom")
# comparisonsCP("MDMKPrandom")

# comparisons4("MDMKPrandom")
comparisons5("MDMKPrandom")
