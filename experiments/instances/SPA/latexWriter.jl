
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
    avg_n = [100, 500, 1000, 1500, 2000, 2500, 3000]
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
            if vars <= seg
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
                if vars <= seg
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
                    if vars <= seg
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
            println("n = $seg , count = $(avg_count[seg]) $m TO = $(avgTO[seg][m])")

        end
    
        println(fout, " \\\\ \\hline")
    end





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
    @assert isdir(work_dir *"/mixed2") "This directory doesn't exist $work_dir !"

    methods = ["bc_rootRelax"] 


    n = 0
    avg_n = [100, 500, 1000, 1500, 2000, 2500, 3000]
    avg_N = Dict(n => 0.0 for n in avg_n)
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
            if vars <= seg
                avg_count[seg] += 1 ; avg_m[seg] += constr ; avg_N[seg] += vars
                break
            end
        end

        # ∀ method 
        for m in methods
            if isfile(work_dir *"/mixed2"* "/" * m * "/" * file)
                include(work_dir *"/mixed2"* "/" * m * "/" * file)
                push!(times, total_times_used); push!(pts, total_nodes)

                for seg in avg_n
                    if vars <= seg
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
        avg_N[seg] = round(avg_N[seg]/avg_count[seg], digits = 2)
        
        for m in methods
            avgT[seg][m] = round(avgT[seg][m]/avg_count[seg], digits = 2); avgY[seg][m] = round(avgY[seg][m]/avg_count[seg], digits = 2) 
        end
    end

    for seg in avg_n

        for m in methods
            println("n = $seg $(avg_N[seg]), count = $(avg_count[seg]) $m TO = $(avgTO[seg][m]) , time $(avgT[seg][m]) , node $(avgY[seg][m])")

        end
    
    end


end


function final_table(instances::String)
    work_dir = "../../results/" * instances
    @assert isdir(work_dir) "This directory doesn't exist $work_dir !"

    fout = open(work_dir * "/FinalTable.tex", "w")

    latex = raw"""\begin{center}
\begin{scriptsize}
\begin{longtable}{lrrrrrrrrr}
\caption{A performance comparison between the $\epsilon$-constraint and BOBLB\&B\&C algorithms.} \\

\toprule
~& ~& ~ & \textbf{\tiny{$\epsilon$-constraint}}  & \multicolumn{2}{c}{\textbf{B\&B}} & \multicolumn{2}{c}{\textbf{EPB B\&C (ISC) ($|\Lambda|$)}} &  \multicolumn{2}{c}{\textbf{Cut\&Branch}} 
\\
\cmidrule(r){4-4} \cmidrule(r){5-6} \cmidrule(r){7-8} \cmidrule(r){9-10} 
\textbf{Instance} & \textbf{n} & \textbf{m}  &\textbf{Time(s)}  & \textbf{Time(s)} &\textbf{\#Nodes} & \textbf{Time(s)} &\textbf{\#Nodes}& \textbf{Time(s)} &\textbf{\#Nodes} \\
\midrule

\endfirsthead
    """
    println(fout, latex)

    # ∀ file in epsilon
    for file in readdir(work_dir * "/epsilon/")
        if split(file, ".")[end] == "png"
            continue
        end

        include(work_dir * "/epsilon/" * file)

        name_seg = split(file, "_")
        for s in name_seg[1:end-1]
            print(fout, s * "\\_")
        end
        print(fout, name_seg[end] * " & ")
        print(fout, string(vars) * " & " * string(constr) * " & ")
        


        # epsilon method 
        if isfile(work_dir * "/epsilon/" * file)
            include(work_dir * "/epsilon/" * file)
            if total_times_used >= 3600.0
                print(fout, " TL & ")
            else
                print(fout, total_times_used, " & ")
            end

        else
            print(fout, " - & ")
        end

        # bb method 
        if isfile(work_dir * "/bb/" * file)
            include(work_dir * "/bb/" * file)
            if total_times_used >= 3600.0
                print(fout, " TL & ")
            else
                print(fout, total_times_used, " & ")
            end
            print(fout, total_nodes, " & ")
            
        else
            print(fout, " - & - & ")
        end
    

        # BC λ
        if isfile(work_dir * "/lambda_limit/2/bc_rootRelaxEPB/"  * file)
            include(work_dir * "/lambda_limit/2/bc_rootRelaxEPB/" * file)

            if total_times_used >= 3600.0
                print(fout, " TL & ")
            else
                print(fout, total_times_used, " & ")
            end
            print(fout, total_nodes, " & ")

        else
            print(fout, " - & - & ")
        end
    
        # cut&branch 
        if isfile(work_dir * "/mixed2/bc_rootRelax/"  * file)
            include(work_dir * "/mixed2/bc_rootRelax/" * file)
            if total_times_used >= 3600.0
                print(fout, " TL & ")
            else
                print(fout, total_times_used, " & ")
            end
            print(fout, total_nodes)

        else
            print(fout, " - & -")
        end

        println(fout, "\\\\")



    end

    latex = raw"""\bottomrule
\end{longtable}
\end{scriptsize}
\end{center}
"""
    println(fout, latex)

    
    close(fout)

end


# comparisons("SPA/BOSPA")
# Cut_Branch("SPA/BOSPA")


final_table("SPA/BOSPA")