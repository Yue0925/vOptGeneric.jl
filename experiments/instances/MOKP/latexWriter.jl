
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

    n = 0
    avg_n = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    avg_count = Dict(n=> 0 for n in avg_n)
    avg_m = Dict(n => 0.0 for n in avg_n)
    avgT = Dict(n => Dict(k => 0.0 for k in methods) for n in avg_n ); 
    avgY = Dict(n=> Dict(k => 0.0 for k in methods) for n in avg_n); 
    avgTO = Dict(n=> Dict(k => 0 for k in methods) for n in avg_n)


    # ∀ file in epsilon
    for file in readdir(work_dir * "/epsilon/")
        if split(file, ".")[end] == "png"
            continue
        end

        times = [] ; pts = []
        include(work_dir * "/epsilon/" * file)

        for seg in avg_n
            if vars == seg
                avg_count[seg] += 1 ; avg_m[seg] += constr
                break
            end
        end

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

            for seg in avg_n
                if vars == seg
                    total_times_used > 3600.0 ? avgTO[seg][m] += 1 : nothing
                    avgT[seg][m] += total_times_used ; avgY[seg][m] += size_Y_N
                    break
                end
            end

        else
            push!(times, -1); push!(pts, -1)
        end

        for m in methods[2:end]
            if isfile(work_dir * "/" * m * "/" * file)
                include(work_dir * "/" * m * "/" * file)
                push!(times, total_times_used); push!(pts, total_nodes)

                for seg in avg_n
                    if vars == seg
                        total_times_used > 3600.0 ? avgTO[seg][m] += 1 : nothing
                        avgT[seg][m] += total_times_used ; avgY[seg][m] += total_nodes
                        break
                    end
                end
                
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

    for seg in avg_n
        avg_m[seg] = round(avg_m[seg]/avg_count[seg], digits = 2)
        for m in methods
            avgT[seg][m] = round(avgT[seg][m]/avg_count[seg], digits = 2); avgY[seg][m] = round(avgY[seg][m]/avg_count[seg], digits = 2) 
        end
    end

    for seg in avg_n
        print(fout, "\\hline avg & " * string(seg) * " & " * string(avg_m[seg]) )

        for m in methods
            print(fout, " &  " * string(avgT[seg][m]) * " & " * string(avgY[seg][m]))
            println("n = $seg , count = $(avg_count[seg]), ($m) TO = $(avgTO[seg][m])")

        end
    
        println(fout, " \\\\ \\hline")
    end

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
    avg_n = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    avg_count = Dict(n=> 0 for n in avg_n)
    avg_m = Dict(n => 0.0 for n in avg_n)

    avgSP = Dict(n => Dict(k => 0.0 for k in methods) for n in avg_n ); 
    avgMP = Dict(n => Dict(k => 0.0 for k in methods) for n in avg_n )


    # ∀ file in dicho
    for file in readdir(work_dir * "/epsilon/")
        if split(file, ".")[end] == "png"
            continue
        end

        sp = [] ; mp = []

        include(work_dir * "/epsilon/" * file)
        for seg in avg_n
            if vars == seg
                avg_count[seg] += 1 ; avg_m[seg] += constr
                break
            end
        end

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

                for seg in avg_n
                    if vars == seg
                        avgSP[seg][m] += sp_cuts ; avgMP[seg][m] += mp_cuts
                        break
                    end
                end
    
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

    for seg in avg_n
        avg_m[seg] = round(avg_m[seg]/avg_count[seg], digits = 2)
        for m in methods
            avgSP[seg][m] = round(avgSP[seg][m]/avg_count[seg], digits = 2); avgMP[seg][m] = round(avgMP[seg][m]/avg_count[seg], digits = 2) 
        end
    end

    for seg in avg_n
        print(fout, "\\hline avg & " * string(seg) * " & " * string(avg_m[seg]) )

        for m in methods
            print(fout, " &  " * string(avgSP[seg][m]) * " & " * string(avgMP[seg][m]))
        end
    
        println(fout, " \\\\ \\hline")
    end


    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout, latex)
    println(fout, "\\caption{Comparison of the B\\&B algorithms performances for instances $instances .}")
    println(fout, "%\\label{tab:table_compareBB_$instances }")
    println(fout, "\\end{table}")
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
    avg_n = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    avg_count = Dict(n=> 0 for n in avg_n)
    avg_m = Dict(n => 0.0 for n in avg_n)
    avgT = Dict(n => Dict(k => 0.0 for k in methods) for n in avg_n ); 
    avgY = Dict(n=> Dict(k => 0.0 for k in methods) for n in avg_n); 
    avgTO = Dict(n=> Dict(k => 0 for k in methods) for n in avg_n)

    # ∀ file in dicho
    for file in readdir(work_dir * "/epsilon/")
        if split(file, ".")[end] == "png"
            continue
        end

        times = [] ; pts = []
        include(work_dir * "/epsilon/" * file)

        for seg in avg_n
            if vars == seg
                avg_count[seg] += 1 ; avg_m[seg] += constr
                break
            end
        end
 
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

                for seg in avg_n
                    if vars == seg
                        total_times_used > 3600.0 ? avgTO[seg][m] += 1 : nothing
                        avgT[seg][m] += total_times_used ; avgY[seg][m] += total_nodes
                        break
                    end
                end

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
    
    for seg in avg_n
        avg_m[seg] = round(avg_m[seg]/avg_count[seg], digits = 2)
        for m in methods
            avgT[seg][m] = round(avgT[seg][m]/avg_count[seg], digits = 2); avgY[seg][m] = round(avgY[seg][m]/avg_count[seg], digits = 2) 
        end
    end

    for seg in avg_n
        print(fout, "\\hline avg & " * string(seg) * " & " * string(avg_m[seg]) )

        for m in methods
            print(fout, " &  " * string(avgT[seg][m]) * " & " * string(avgY[seg][m]))
            println("n = $seg , count = $(avg_count[seg]), ($m) TO = $(avgTO[seg][m])")

        end
    
        println(fout, " \\\\ \\hline")
    end


    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout, latex)
    println(fout, "\\caption{Comparison of the different algorithms performances for instances $instances .}")
    println(fout, "%\\label{tab:table_compare_$instances }")
    println(fout, "\\end{sidewaystable}")
    close(fout)

end


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
    avg_n = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    avg_count = Dict(n=> 0 for n in avg_n)
    avg_m = Dict(n => 0.0 for n in avg_n)
    avgT = Dict(n => Dict(k => 0.0 for k in methods) for n in avg_n ); 
    avgY = Dict(n=> Dict(k => 0.0 for k in methods) for n in avg_n); 
    avgTO = Dict(n=> Dict(k => 0 for k in methods) for n in avg_n)

    # ∀ file in dicho
    for file in readdir(work_dir * "/epsilon/")
        if split(file, ".")[end] == "png"
            continue
        end

        times = [] ; pts = []
        include(work_dir * "/epsilon/" * file)

        for seg in avg_n
            if vars == seg
                avg_count[seg] += 1 ; avg_m[seg] += constr
                break
            end
        end

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

                for seg in avg_n
                    if vars == seg
                        total_times_used > 3600.0 ? avgTO[seg][m] += 1 : nothing
                        avgT[seg][m] += total_times_used ; avgY[seg][m] += total_nodes
                        break
                    end
                end

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


    for seg in avg_n
        avg_m[seg] = round(avg_m[seg]/avg_count[seg], digits = 2)
        for m in methods
            avgT[seg][m] = round(avgT[seg][m]/avg_count[seg], digits = 2); avgY[seg][m] = round(avgY[seg][m]/avg_count[seg], digits = 2) 
        end
    end

    for seg in avg_n
        print(fout, "\\hline avg & " * string(seg) * " & " * string(avg_m[seg]) )

        for m in methods
            print(fout, " &  " * string(avgT[seg][m]) * " & " * string(avgY[seg][m]))
            println("n = $seg , count = $(avg_count[seg]), ($m) TO = $(avgTO[seg][m])")

        end
    
        println(fout, " \\\\ \\hline")
    end

    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout, latex)
    println(fout, "\\caption{Comparison of the different algorithms performances for instances $instances .}")
    println(fout, "%\\label{tab:table_compare_$instances }")
    println(fout, "\\end{sidewaystable}")
    close(fout)

end


function Cut_Branch(instances::String)
    work_dir = "../../results/" * instances
    @assert isdir(work_dir * "/mixed2") "This directory doesn't exist $work_dir !"

    methods = ["bc_rootRelax"] 

    n = 0
    avg_n = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    avg_count = Dict(n=> 0 for n in avg_n)
    avg_m = Dict(n => 0.0 for n in avg_n)
    avgT = Dict(n => Dict(k => 0.0 for k in methods) for n in avg_n ); 
    avgY = Dict(n=> Dict(k => 0.0 for k in methods) for n in avg_n); 
    avgTO = Dict(n=> Dict(k => 0 for k in methods) for n in avg_n)


    # ∀ file in epsilon
    for file in readdir(work_dir * "/epsilon/")
        if split(file, ".")[end] == "png"
            continue
        end

        times = [] ; pts = []
        include(work_dir * "/epsilon/" * file)

        for seg in avg_n
            if vars == seg
                avg_count[seg] += 1 ; avg_m[seg] += constr
                break
            end
        end      

        for m in methods
            if isfile(work_dir *"/mixed2" * "/" * m * "/" * file)
                include(work_dir *"/mixed2" * "/" * m * "/" * file)
                push!(times, total_times_used); push!(pts, total_nodes)

                for seg in avg_n
                    if vars == seg
                        total_times_used > 3600.0 ? avgTO[seg][m] += 1 : nothing
                        avgT[seg][m] += total_times_used ; avgY[seg][m] += total_nodes
                        break
                    end
                end
                
            else
                push!(times, -1); push!(pts, -1)
            end
        end

    end

    for seg in avg_n
        avg_m[seg] = round(avg_m[seg]/avg_count[seg], digits = 2)
        for m in methods
            avgT[seg][m] = round(avgT[seg][m]/avg_count[seg], digits = 2); avgY[seg][m] = round(avgY[seg][m]/avg_count[seg], digits = 2) 
        end
    end

    for seg in avg_n

        for m in methods
            println("n = $seg , count = $(avg_count[seg]), ($m) TO = $(avgTO[seg][m]) time $(avgT[seg][m]) , node $(avgY[seg][m])")

        end
    
    end

end





# comparisons("MOKP")
# comparisonsCP("MOKP")
# comparisons4("MOKP")
# comparisons5("MOKP")
Cut_Branch("MOKP")

# n = 10 , count = 30, (epsilon) TO = 0
# n = 10 , count = 30, (bb) TO = 0
# n = 10 , count = 30, (bc) TO = 0
# n = 10 , count = 30, (bc_rootRelax) TO = 0
# n = 10 , count = 30, (bb_EPB) TO = 0
# n = 10 , count = 30, (bc_EPB) TO = 0
# n = 10 , count = 30, (bc_rootRelaxEPB) TO = 0
# n = 10 , count = 30, (bc_rootRelaxCP) TO = 0
# n = 10 , count = 30, (bc_rootRelaxCPEPB) TO = 0
# n = 20 , count = 30, (epsilon) TO = 0
# n = 20 , count = 30, (bb) TO = 0
# n = 20 , count = 30, (bc) TO = 0
# n = 20 , count = 30, (bc_rootRelax) TO = 0
# n = 20 , count = 30, (bb_EPB) TO = 0
# n = 20 , count = 30, (bc_EPB) TO = 0
# n = 20 , count = 30, (bc_rootRelaxEPB) TO = 0
# n = 20 , count = 30, (bc_rootRelaxCP) TO = 0
# n = 20 , count = 30, (bc_rootRelaxCPEPB) TO = 0
# n = 30 , count = 20, (epsilon) TO = 0
# n = 30 , count = 20, (bb) TO = 0
# n = 30 , count = 20, (bc) TO = 0
# n = 30 , count = 20, (bc_rootRelax) TO = 0
# n = 30 , count = 20, (bb_EPB) TO = 0
# n = 30 , count = 20, (bc_EPB) TO = 0
# n = 30 , count = 20, (bc_rootRelaxEPB) TO = 0
# n = 30 , count = 20, (bc_rootRelaxCP) TO = 0
# n = 30 , count = 20, (bc_rootRelaxCPEPB) TO = 0
# n = 40 , count = 20, (epsilon) TO = 0
# n = 40 , count = 20, (bb) TO = 0
# n = 40 , count = 20, (bc) TO = 0
# n = 40 , count = 20, (bc_rootRelax) TO = 0
# n = 40 , count = 20, (bb_EPB) TO = 0
# n = 40 , count = 20, (bc_EPB) TO = 0
# n = 40 , count = 20, (bc_rootRelaxEPB) TO = 0
# n = 40 , count = 20, (bc_rootRelaxCP) TO = 0
# n = 40 , count = 20, (bc_rootRelaxCPEPB) TO = 0
# n = 50 , count = 10, (epsilon) TO = 0
# n = 50 , count = 10, (bb) TO = 0
# n = 50 , count = 10, (bc) TO = 5
# n = 50 , count = 10, (bc_rootRelax) TO = 0
# n = 50 , count = 10, (bb_EPB) TO = 5
# n = 50 , count = 10, (bc_EPB) TO = 5
# n = 50 , count = 10, (bc_rootRelaxEPB) TO = 0
# n = 50 , count = 10, (bc_rootRelaxCP) TO = 0
# n = 50 , count = 10, (bc_rootRelaxCPEPB) TO = 0
# n = 60 , count = 10, (epsilon) TO = 0
# n = 60 , count = 10, (bb) TO = 9
# n = 60 , count = 10, (bc) TO = 10
# n = 60 , count = 10, (bc_rootRelax) TO = 0
# n = 60 , count = 10, (bb_EPB) TO = 10
# n = 60 , count = 10, (bc_EPB) TO = 10
# n = 60 , count = 10, (bc_rootRelaxEPB) TO = 0
# n = 60 , count = 10, (bc_rootRelaxCP) TO = 0
# n = 60 , count = 10, (bc_rootRelaxCPEPB) TO = 0
# n = 70 , count = 10, (epsilon) TO = 0
# n = 70 , count = 10, (bb) TO = 10
# n = 70 , count = 10, (bc) TO = 10
# n = 70 , count = 10, (bc_rootRelax) TO = 0
# n = 70 , count = 10, (bb_EPB) TO = 10
# n = 70 , count = 10, (bc_EPB) TO = 10
# n = 70 , count = 10, (bc_rootRelaxEPB) TO = 0
# n = 70 , count = 10, (bc_rootRelaxCP) TO = 0
# n = 70 , count = 10, (bc_rootRelaxCPEPB) TO = 0
# n = 80 , count = 10, (epsilon) TO = 0
# n = 80 , count = 10, (bb) TO = 10
# n = 80 , count = 10, (bc) TO = 10
# n = 80 , count = 10, (bc_rootRelax) TO = 3
# n = 80 , count = 10, (bb_EPB) TO = 10
# n = 80 , count = 10, (bc_EPB) TO = 10
# n = 80 , count = 10, (bc_rootRelaxEPB) TO = 0
# n = 80 , count = 10, (bc_rootRelaxCP) TO = 3
# n = 80 , count = 10, (bc_rootRelaxCPEPB) TO = 0
# n = 90 , count = 10, (epsilon) TO = 0
# n = 90 , count = 10, (bb) TO = 10
# n = 90 , count = 10, (bc) TO = 10
# n = 90 , count = 10, (bc_rootRelax) TO = 0
# n = 90 , count = 10, (bb_EPB) TO = 10
# n = 90 , count = 10, (bc_EPB) TO = 10
# n = 90 , count = 10, (bc_rootRelaxEPB) TO = 0
# n = 90 , count = 10, (bc_rootRelaxCP) TO = 0
# n = 90 , count = 10, (bc_rootRelaxCPEPB) TO = 0
# n = 100 , count = 10, (epsilon) TO = 0
# n = 100 , count = 10, (bb) TO = 10
# n = 100 , count = 10, (bc) TO = 10
# n = 100 , count = 10, (bc_rootRelax) TO = 7
# n = 100 , count = 10, (bb_EPB) TO = 10
# n = 100 , count = 10, (bc_EPB) TO = 10
# n = 100 , count = 10, (bc_rootRelaxEPB) TO = 0
# n = 100 , count = 10, (bc_rootRelaxCP) TO = 7
# n = 100 , count = 10, (bc_rootRelaxCPEPB) TO = 0
