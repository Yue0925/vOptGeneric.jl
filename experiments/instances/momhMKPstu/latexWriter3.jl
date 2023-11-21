
function comparisonsLambdaLimits(instances::String)
    work_dir = "../../results/" * instances * "/lambda_limit"
    @assert isdir(work_dir) "This directory doesn't exist $work_dir !"

    fout = open(work_dir * "/comparisonLambdaLimitTable.tex", "w")

    latex = raw"""\begin{sidewaystable}[!ht]
    \centering
    \resizebox{\columnwidth}{!}{%
    \begin{tabular}{lccccccccccccccccccc}
    \toprule
    % line 1 
    \textbf{Instance} & \textbf{n} & \textbf{m} & \textbf{$\mathcal{Y}_N$}
    & \multicolumn{4}{c}{\textbf{$|\lambda| = Inf$}} & \multicolumn{4}{c}{\textbf{$|\lambda| = 2$}} & \multicolumn{4}{c}{\textbf{$|\lambda| = 4$}} & \multicolumn{4}{c}{\textbf{$|\lambda| = 6$}} \\
    
    %line 2 
    \cmidrule(r){5-8} \cmidrule(r){9-12} \cmidrule(r){13-16} \cmidrule(r){17-20}
    ~ & ~ & ~ & ~ & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} \\
    
    %line 3
    \cmidrule(r){5-6} \cmidrule(r){7-8} \cmidrule(r){9-10} \cmidrule(r){11-12} \cmidrule(r){13-14} \cmidrule(r){15-16} \cmidrule(r){17-18} \cmidrule(r){19-20}
    ~ & ~ & ~ & ~ & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} \\
    \midrule
    """
    println(fout, latex)

    methods = ["bc_rootRelax", "bc_rootRelaxEPB"] ; record_n = []

    λ_limits = [] ; 
    for folder_λ in readdir(work_dir)
        if split(folder_λ, ".")[end] == "tex" continue end
        push!(λ_limits, folder_λ)
    end
    keys_combo = [] ; 
    for λ in λ_limits
        for m in methods
            push!(keys_combo, λ * "_" * m)
        end
    end


    # ∀ folder_n  
    for folder_n in readdir(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]) ) 
        avgT = Dict(k => 0.0 for k in keys_combo) ; avgY = Dict(k => 0.0 for k in keys_combo)

        count = 0
        avg_n = 0 ; avg_m = 0 ; avg_Y_N = 0

        # ∀ file each line 
        for file in readdir(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]) * "/" * string(folder_n))
            if split(file, ".")[end] == "png"
                continue
            end
            times = [] ; pts = []

            print(fout, file * " & ")

            # write dichotomy result 
            include(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]) * "/" * string(folder_n) * "/" * file)
            print(fout, string(vars) * " & " * string(constr) * " & " * string(size_Y_N) * " & ")

            count += 1
            avg_n += vars ; avg_m += constr ; avg_Y_N += size_Y_N

            for folder_λ in λ_limits
                for m in methods
                    if isfile(work_dir * "/" * string(folder_λ) * "/" * m * "/" * string(folder_n) * "/" * file)
                        include(work_dir * "/" * string(folder_λ) * "/" * m * "/" * string(folder_n) * "/" * file)
                        push!(times, total_times_used); push!(pts, total_nodes)
        
                        avgT[folder_λ * "_" * m] += total_times_used ; avgY[folder_λ * "_" * m] += total_nodes
                    else
                        push!(times, -1); push!(pts, -1)
                    end
                end
            end


            # each line 
            for i=1:length(keys_combo)-1
                if times[i] == -1
                    print(fout, " - & ")
                elseif times[i] >= 3600.0
                    print(fout, " TO & ")
                elseif times[i] == minimum(filter(x -> x > 0 ,times))
                    print(fout, " \\textcolor{blue2}{" * string(times[i]) * "} & ")
                else
                    print(fout, string(times[i]) * " & ")
                end
    
                if pts[i] == -1
                    print(fout, " - & ")
                else
                    print(fout, string(pts[i]) * " & ")
                end
    
            end

            if times[end] == minimum(filter(x -> x > 0 ,times))
                print(fout, " \\textcolor{blue2}{" * string(times[end]) * "} & ")
            elseif times[end] >= 3600.0
                print(fout, " TO & ")
            else
                print(fout, string(times[end]) * " & ") 
            end
            
            println(fout, string(pts[end]) * " \\\\")
    
        end

            
        avg_n = round(Int, avg_n/count) ; avg_m = round(Int, avg_m/count) ; avg_Y_N = round(Int, avg_Y_N/count) 
        for m in keys_combo
            avgT[m] = round(avgT[m]/count, digits = 2); avgY[m] = round(avgY[m]/count, digits = 2) 
        end

        print(fout, "\\cline{1-20} \\textbf{avg} & \\textbf{" * string(avg_n) * "} & \\textbf{" * string(avg_m) * "} & \\textbf{" * string(avg_Y_N) )

        for m in keys_combo
            print(fout, "} & \\textbf{" * string(avgT[m]) * "} & \\textbf{" * string(avgY[m]))
        end

        println(fout, "} \\\\ \\cline{1-20}")


    end

    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout, latex)
    println(fout, "\\caption{cplex cutting LBS non-exhaustive dichotomic concave-convex like algo on instances $instances (\$\\lambda\$ fixed) .}")
    println(fout, "\\label{tab:table_lambda_limits_$instances }")
    println(fout, "\\end{sidewaystable}")
    close(fout)

end




function comparisonsLBSFactor2(instances::String)
    work_dir = "../../results/" * instances * "/rootLBS_factor2"
    @assert isdir(work_dir) "This directory doesn't exist $work_dir !"

    fout = open(work_dir * "/comparisonLBSFactor2.tex", "w")

    latex = raw"""\begin{sidewaystable}[!ht]
    \centering
    \resizebox{\columnwidth}{!}{%
    \begin{tabular}{lccccccccccccccc}
    \toprule
    % line 1 
    \textbf{Instance} & \textbf{n} & \textbf{m} & \textbf{$L_{root}$}
    & \multicolumn{4}{c}{\textbf{$|\lambda| = L/4$}} & \multicolumn{4}{c}{\textbf{$|\lambda| = L/8$}} & \multicolumn{4}{c}{\textbf{$|\lambda| = L/12$}} \\
    
    %line 2 
    \cmidrule(r){5-8} \cmidrule(r){9-12} \cmidrule(r){13-16} 
    ~ & ~ & ~ & ~ & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} \\
    
    %line 3
    \cmidrule(r){5-6} \cmidrule(r){7-8} \cmidrule(r){9-10} \cmidrule(r){11-12} \cmidrule(r){13-14} \cmidrule(r){15-16}
    ~ & ~ & ~ & ~ & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes}\\
    \midrule
    """
    println(fout, latex)

    methods = ["bc_rootRelax", "bc_rootRelaxEPB"] ; record_n = []

    λ_limits = [] ; 
    for folder_λ in readdir(work_dir)
        if split(folder_λ, ".")[end] == "tex" continue end
        push!(λ_limits, folder_λ)
    end
    keys_combo = [] ; 
    for λ in λ_limits
        for m in methods
            push!(keys_combo, λ * "_" * m)
        end
    end


    # ∀ folder_n  
    for folder_n in readdir(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]) ) 
        avgT = Dict(k => 0.0 for k in keys_combo) ; avgY = Dict(k => 0.0 for k in keys_combo)

        count = 0
        avg_n = 0 ; avg_m = 0 ; avg_Y_N = 0

        # ∀ file each line 
        for file in readdir(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]) * "/" * string(folder_n))
            if split(file, ".")[end] == "png"
                continue
            end
            times = [] ; pts = []

            print(fout, file * " & ")

            # write dichotomy result 
            include(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]) * "/" * string(folder_n) * "/" * file)
            print(fout, string(vars) * " & " * string(constr) * " & " * string(rootLBS) * " & ")

            count += 1
            avg_n += vars ; avg_m += constr ; avg_Y_N += rootLBS

            for folder_λ in λ_limits
                for m in methods
                    if isfile(work_dir * "/" * string(folder_λ) * "/" * m * "/" * string(folder_n) * "/" * file)
                        include(work_dir * "/" * string(folder_λ) * "/" * m * "/" * string(folder_n) * "/" * file)
                        push!(times, total_times_used); push!(pts, total_nodes)
        
                        avgT[folder_λ * "_" * m] += total_times_used ; avgY[folder_λ * "_" * m] += total_nodes
                    else
                        push!(times, -1); push!(pts, -1)
                    end
                end
            end


            # each line 
            for i=1:length(keys_combo)-1
                if times[i] == -1
                    print(fout, " - & ")
                elseif times[i] >= 3600.0
                    print(fout, " TO & ")
                elseif times[i] == minimum(filter(x -> x > 0 ,times))
                    print(fout, " \\textcolor{blue2}{" * string(times[i]) * "} & ")
                else
                    print(fout, string(times[i]) * " & ")
                end
    
                if pts[i] == -1
                    print(fout, " - & ")
                else
                    print(fout, string(pts[i]) * " & ")
                end
    
            end

            if times[end] == minimum(filter(x -> x > 0 ,times))
                print(fout, " \\textcolor{blue2}{" * string(times[end]) * "} & ")
            elseif times[end] >= 3600.0
                print(fout, " TO & ")
            else
                print(fout, string(times[end]) * " & ") 
            end
            
            println(fout, string(pts[end]) * " \\\\")
    
        end

            
        avg_n = round(Int, avg_n/count) ; avg_m = round(Int, avg_m/count) ; avg_Y_N = round(Int, avg_Y_N/count) 
        for m in keys_combo
            avgT[m] = round(avgT[m]/count, digits = 2); avgY[m] = round(avgY[m]/count, digits = 2) 
        end

        print(fout, "\\cline{1-16} \\textbf{avg} & \\textbf{" * string(avg_n) * "} & \\textbf{" * string(avg_m) * "} & \\textbf{" * string(avg_Y_N) )

        for m in keys_combo
            print(fout, "} & \\textbf{" * string(avgT[m]) * "} & \\textbf{" * string(avgY[m]))
        end

        println(fout, "} \\\\ \\cline{1-16}")


    end

    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout, latex)
    println(fout, "\\caption{cplex cutting LBS non-exhaustive dichotomic concave-convex like algo on instances $instances (\$\\lambda = L_{root}/K\$ fixed, with lower bound 2, except EPB)  .}")
    println(fout, "\\label{tab:table_LBS_factor2_$instances }")
    println(fout, "\\end{sidewaystable}")
    close(fout)

end




function comparisonsLBSFactor3(instances::String)
    work_dir = "../../results/" * instances * "/rootLBS_factor3"
    @assert isdir(work_dir) "This directory doesn't exist $work_dir !"

    fout = open(work_dir * "/comparisonLBSFactor2.tex", "w")

    latex = raw"""\begin{sidewaystable}[!ht]
    \centering
    \resizebox{\columnwidth}{!}{%
    \begin{tabular}{lccccccccccccccc}
    \toprule
    % line 1 
    \textbf{Instance} & \textbf{n} & \textbf{m} & \textbf{$L_{root}$}
    & \multicolumn{4}{c}{\textbf{$|\lambda| = L/4$}} & \multicolumn{4}{c}{\textbf{$|\lambda| = L/8$}} & \multicolumn{4}{c}{\textbf{$|\lambda| = L/12$}} \\
    
    %line 2 
    \cmidrule(r){5-8} \cmidrule(r){9-12} \cmidrule(r){13-16} 
    ~ & ~ & ~ & ~ & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} \\
    
    %line 3
    \cmidrule(r){5-6} \cmidrule(r){7-8} \cmidrule(r){9-10} \cmidrule(r){11-12} \cmidrule(r){13-14} \cmidrule(r){15-16}
    ~ & ~ & ~ & ~ & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes}\\
    \midrule
    """
    println(fout, latex)

    methods = ["bc_rootRelax", "bc_rootRelaxEPB"] ; record_n = []

    λ_limits = [] ; 
    for folder_λ in readdir(work_dir)
        if split(folder_λ, ".")[end] == "tex" continue end
        push!(λ_limits, folder_λ)
    end
    keys_combo = [] ; 
    for λ in λ_limits
        for m in methods
            push!(keys_combo, λ * "_" * m)
        end
    end


    # ∀ folder_n  
    for folder_n in readdir(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]) ) 
        avgT = Dict(k => 0.0 for k in keys_combo) ; avgY = Dict(k => 0.0 for k in keys_combo)

        count = 0
        avg_n = 0 ; avg_m = 0 ; avg_Y_N = 0

        # ∀ file each line 
        for file in readdir(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]) * "/" * string(folder_n))
            if split(file, ".")[end] == "png"
                continue
            end
            times = [] ; pts = []

            print(fout, file * " & ")

            # write dichotomy result 
            include(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]) * "/" * string(folder_n) * "/" * file)
            print(fout, string(vars) * " & " * string(constr) * " & " * string(rootLBS) * " & ")

            count += 1
            avg_n += vars ; avg_m += constr ; avg_Y_N += rootLBS

            for folder_λ in λ_limits
                for m in methods
                    if isfile(work_dir * "/" * string(folder_λ) * "/" * m * "/" * string(folder_n) * "/" * file)
                        include(work_dir * "/" * string(folder_λ) * "/" * m * "/" * string(folder_n) * "/" * file)
                        push!(times, total_times_used); push!(pts, total_nodes)
        
                        avgT[folder_λ * "_" * m] += total_times_used ; avgY[folder_λ * "_" * m] += total_nodes
                    else
                        push!(times, -1); push!(pts, -1)
                    end
                end
            end


            # each line 
            for i=1:length(keys_combo)-1
                if times[i] == -1
                    print(fout, " - & ")
                elseif times[i] >= 3600.0
                    print(fout, " TO & ")
                elseif times[i] == minimum(filter(x -> x > 0 ,times))
                    print(fout, " \\textcolor{blue2}{" * string(times[i]) * "} & ")
                else
                    print(fout, string(times[i]) * " & ")
                end
    
                if pts[i] == -1
                    print(fout, " - & ")
                else
                    print(fout, string(pts[i]) * " & ")
                end
    
            end

            if times[end] == minimum(filter(x -> x > 0 ,times))
                print(fout, " \\textcolor{blue2}{" * string(times[end]) * "} & ")
            elseif times[end] >= 3600.0
                print(fout, " TO & ")
            else
                print(fout, string(times[end]) * " & ") 
            end
            
            println(fout, string(pts[end]) * " \\\\")
    
        end

            
        avg_n = round(Int, avg_n/count) ; avg_m = round(Int, avg_m/count) ; avg_Y_N = round(Int, avg_Y_N/count) 
        for m in keys_combo
            avgT[m] = round(avgT[m]/count, digits = 2); avgY[m] = round(avgY[m]/count, digits = 2) 
        end

        print(fout, "\\cline{1-16} \\textbf{avg} & \\textbf{" * string(avg_n) * "} & \\textbf{" * string(avg_m) * "} & \\textbf{" * string(avg_Y_N) )

        for m in keys_combo
            print(fout, "} & \\textbf{" * string(avgT[m]) * "} & \\textbf{" * string(avgY[m]))
        end

        println(fout, "} \\\\ \\cline{1-16}")


    end

    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout, latex)
    println(fout, "\\caption{cplex cutting LBS non-exhaustive dichotomic concave-convex like algo on instances $instances (\$\\lambda = L_{root}/K\$ fixed, with lower bound 2, except EPB)  .}")
    println(fout, "\\label{tab:table_LBS_factor3_$instances }")
    println(fout, "\\end{sidewaystable}")
    close(fout)

end


function comparisonsLambdaChordal(instances::String)
    work_dir = "../../results/" * instances * "/lambda_chordal"
    @assert isdir(work_dir) "This directory doesn't exist $work_dir !"

    fout = open(work_dir * "/comparisonLambdaChordalTable.tex", "w")

    latex = raw"""\begin{sidewaystable}[!ht]
    \centering
    \resizebox{\columnwidth}{!}{%
    \begin{tabular}{lccccccccccccccc}
    \toprule
    % line 1 
    \textbf{Instance} & \textbf{n} & \textbf{m} & \textbf{$\mathcal{Y}_N$}
    & \multicolumn{4}{c}{\textbf{$|\lambda| = 2$}} & \multicolumn{4}{c}{\textbf{$|\lambda| = 4$}} & \multicolumn{4}{c}{\textbf{$|\lambda| = 6$}} \\
    
    %line 2 
    \cmidrule(r){5-8} \cmidrule(r){9-12} \cmidrule(r){13-16} 
    ~ & ~ & ~ & ~ & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} \\
    
    %line 3
    \cmidrule(r){5-6} \cmidrule(r){7-8} \cmidrule(r){9-10} \cmidrule(r){11-12} \cmidrule(r){13-14} \cmidrule(r){15-16}
    ~ & ~ & ~ & ~ & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} \\
    \midrule
    """
    println(fout, latex)

    methods = ["bc_rootRelax", "bc_rootRelaxEPB"] ; record_n = []

    λ_limits = [] ; 
    for folder_λ in readdir(work_dir)
        if split(folder_λ, ".")[end] == "tex" continue end
        push!(λ_limits, folder_λ)
    end
    keys_combo = [] ; 
    for λ in λ_limits
        for m in methods
            push!(keys_combo, λ * "_" * m)
        end
    end


    # ∀ folder_n  
    for folder_n in readdir(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]) ) 
        avgT = Dict(k => 0.0 for k in keys_combo) ; avgY = Dict(k => 0.0 for k in keys_combo)

        count = 0
        avg_n = 0 ; avg_m = 0 ; avg_Y_N = 0

        # ∀ file each line 
        for file in readdir(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]) * "/" * string(folder_n))
            if split(file, ".")[end] == "png"
                continue
            end
            times = [] ; pts = []

            print(fout, file * " & ")

            # write dichotomy result 
            include(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]) * "/" * string(folder_n) * "/" * file)
            print(fout, string(vars) * " & " * string(constr) * " & " * string(size_Y_N) * " & ")

            count += 1
            avg_n += vars ; avg_m += constr ; avg_Y_N += size_Y_N

            for folder_λ in λ_limits
                for m in methods
                    if isfile(work_dir * "/" * string(folder_λ) * "/" * m * "/" * string(folder_n) * "/" * file)
                        include(work_dir * "/" * string(folder_λ) * "/" * m * "/" * string(folder_n) * "/" * file)
                        push!(times, total_times_used); push!(pts, total_nodes)
        
                        avgT[folder_λ * "_" * m] += total_times_used ; avgY[folder_λ * "_" * m] += total_nodes
                    else
                        push!(times, -1); push!(pts, -1)
                    end
                end
            end


            # each line 
            for i=1:length(keys_combo)-1
                if times[i] == -1
                    print(fout, " - & ")
                elseif times[i] >= 3600.0
                    print(fout, " TO & ")
                elseif times[i] == minimum(filter(x -> x > 0 ,times))
                    print(fout, " \\textcolor{blue2}{" * string(times[i]) * "} & ")
                else
                    print(fout, string(times[i]) * " & ")
                end
    
                if pts[i] == -1
                    print(fout, " - & ")
                else
                    print(fout, string(pts[i]) * " & ")
                end
    
            end

            if times[end] == minimum(filter(x -> x > 0 ,times))
                print(fout, " \\textcolor{blue2}{" * string(times[end]) * "} & ")
            elseif times[end] >= 3600.0
                print(fout, " TO & ")
            else
                print(fout, string(times[end]) * " & ") 
            end
            
            println(fout, string(pts[end]) * " \\\\")
    
        end

            
        avg_n = round(Int, avg_n/count) ; avg_m = round(Int, avg_m/count) ; avg_Y_N = round(Int, avg_Y_N/count) 
        for m in keys_combo
            avgT[m] = round(avgT[m]/count, digits = 2); avgY[m] = round(avgY[m]/count, digits = 2) 
        end

        print(fout, "\\cline{1-16} \\textbf{avg} & \\textbf{" * string(avg_n) * "} & \\textbf{" * string(avg_m) * "} & \\textbf{" * string(avg_Y_N) )

        for m in keys_combo
            print(fout, "} & \\textbf{" * string(avgT[m]) * "} & \\textbf{" * string(avgY[m]))
        end

        println(fout, "} \\\\ \\cline{1-16}")


    end

    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout, latex)
    println(fout, "\\caption{cplex cutting LBS non-exhaustive chordal like search algo on instances $instances (\$\\lambda\$ fixed) .}")
    println(fout, "\\label{tab:table_lambda_chordal_$instances }")
    println(fout, "\\end{sidewaystable}")
    close(fout)

end


function comparisonsLambdaDynamic(instances::String)
    work_dir = "../../results/" * instances * "/lambda_dynamic"
    @assert isdir(work_dir) "This directory doesn't exist $work_dir !"

    fout = open(work_dir * "/comparisonLambdaDynamicTable.tex", "w")

    latex = raw"""\begin{sidewaystable}[!ht]
    \centering
    \resizebox{\columnwidth}{!}{%
    \begin{tabular}{lccccccccccccccc}
    \toprule
    % line 1 
    \textbf{Instance} & \textbf{n} & \textbf{m} & \textbf{$\mathcal{Y}_N$}
    & \multicolumn{4}{c}{\textbf{$|\lambda| = 2$}} & \multicolumn{4}{c}{\textbf{$|\lambda| = 4$}} & \multicolumn{4}{c}{\textbf{$|\lambda| = 6$}} \\
    
    %line 2 
    \cmidrule(r){5-8} \cmidrule(r){9-12} \cmidrule(r){13-16} 
    ~ & ~ & ~ & ~ & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} \\
    
    %line 3
    \cmidrule(r){5-6} \cmidrule(r){7-8} \cmidrule(r){9-10} \cmidrule(r){11-12} \cmidrule(r){13-14} \cmidrule(r){15-16}
    ~ & ~ & ~ & ~ & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} \\
    \midrule
    """
    println(fout, latex)

    methods = ["bc_rootRelax", "bc_rootRelaxEPB"] ; record_n = []

    λ_limits = [] ; 
    for folder_λ in readdir(work_dir)
        if split(folder_λ, ".")[end] == "tex" continue end
        push!(λ_limits, folder_λ)
    end
    keys_combo = [] ; 
    for λ in λ_limits
        for m in methods
            push!(keys_combo, λ * "_" * m)
        end
    end


    # ∀ folder_n  
    for folder_n in readdir(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]) ) 
        avgT = Dict(k => 0.0 for k in keys_combo) ; avgY = Dict(k => 0.0 for k in keys_combo)

        count = 0
        avg_n = 0 ; avg_m = 0 ; avg_Y_N = 0

        # ∀ file each line 
        for file in readdir(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]) * "/" * string(folder_n))
            if split(file, ".")[end] == "png"
                continue
            end
            times = [] ; pts = []

            print(fout, file * " & ")

            # write dichotomy result 
            include(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]) * "/" * string(folder_n) * "/" * file)
            print(fout, string(vars) * " & " * string(constr) * " & " * string(size_Y_N) * " & ")

            count += 1
            avg_n += vars ; avg_m += constr ; avg_Y_N += size_Y_N

            for folder_λ in λ_limits
                for m in methods
                    if isfile(work_dir * "/" * string(folder_λ) * "/" * m * "/" * string(folder_n) * "/" * file)
                        include(work_dir * "/" * string(folder_λ) * "/" * m * "/" * string(folder_n) * "/" * file)
                        push!(times, total_times_used); push!(pts, total_nodes)
        
                        avgT[folder_λ * "_" * m] += total_times_used ; avgY[folder_λ * "_" * m] += total_nodes
                    else
                        push!(times, -1); push!(pts, -1)
                    end
                end
            end


            # each line 
            for i=1:length(keys_combo)-1
                if times[i] == -1
                    print(fout, " - & ")
                elseif times[i] >= 3600.0
                    print(fout, " TO & ")
                elseif times[i] == minimum(filter(x -> x > 0 ,times))
                    print(fout, " \\textcolor{blue2}{" * string(times[i]) * "} & ")
                else
                    print(fout, string(times[i]) * " & ")
                end
    
                if pts[i] == -1
                    print(fout, " - & ")
                else
                    print(fout, string(pts[i]) * " & ")
                end
    
            end

            if times[end] == minimum(filter(x -> x > 0 ,times))
                print(fout, " \\textcolor{blue2}{" * string(times[end]) * "} & ")
            elseif times[end] >= 3600.0
                print(fout, " TO & ")
            else
                print(fout, string(times[end]) * " & ") 
            end
            
            println(fout, string(pts[end]) * " \\\\")
    
        end

            
        avg_n = round(Int, avg_n/count) ; avg_m = round(Int, avg_m/count) ; avg_Y_N = round(Int, avg_Y_N/count) 
        for m in keys_combo
            avgT[m] = round(avgT[m]/count, digits = 2); avgY[m] = round(avgY[m]/count, digits = 2) 
        end

        print(fout, "\\cline{1-16} \\textbf{avg} & \\textbf{" * string(avg_n) * "} & \\textbf{" * string(avg_m) * "} & \\textbf{" * string(avg_Y_N) )

        for m in keys_combo
            print(fout, "} & \\textbf{" * string(avgT[m]) * "} & \\textbf{" * string(avgY[m]))
        end

        println(fout, "} \\\\ \\cline{1-16}")


    end

    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout, latex)
    println(fout, "\\caption{cplex cutting LBS non-exhaustive dynamic search algo on instances $instances (\$\\lambda\$ fixed) .}")
    println(fout, "\\label{tab:table_lambda_dynamic_$instances }")
    println(fout, "\\end{sidewaystable}")
    close(fout)

end


function comparisonPredLBSFactor(instances::String)
    work_dir = "../../results/" * instances * "/predLBS_factor"
    @assert isdir(work_dir) "This directory doesn't exist $work_dir !"

    fout = open(work_dir * "/comparisonPredLBSFactor.tex", "w")

    latex = raw"""\begin{sidewaystable}[!ht]
    \centering
    \resizebox{\columnwidth}{!}{%
    \begin{tabular}{lccccccccccccccc}
    \toprule
    % line 1 
    \textbf{Instance} & \textbf{n} & \textbf{m} & \textbf{$L_{pred}$}
    & \multicolumn{4}{c}{\textbf{$|\lambda| = L/4$}} & \multicolumn{4}{c}{\textbf{$|\lambda| = L/8$}} & \multicolumn{4}{c}{\textbf{$|\lambda| = L/12$}} \\
    
    %line 2 
    \cmidrule(r){5-8} \cmidrule(r){9-12} \cmidrule(r){13-16} 
    ~ & ~ & ~ & ~ & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} \\
    
    %line 3
    \cmidrule(r){5-6} \cmidrule(r){7-8} \cmidrule(r){9-10} \cmidrule(r){11-12} \cmidrule(r){13-14} \cmidrule(r){15-16}
    ~ & ~ & ~ & ~ & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes}\\
    \midrule
    """
    println(fout, latex)

    methods = ["bc_rootRelax", "bc_rootRelaxEPB"] ; record_n = []

    λ_limits = [] ; 
    for folder_λ in readdir(work_dir)
        if split(folder_λ, ".")[end] == "tex" continue end
        push!(λ_limits, folder_λ)
    end
    keys_combo = [] ; 
    for λ in λ_limits
        for m in methods
            push!(keys_combo, λ * "_" * m)
        end
    end


    # ∀ folder_n  
    for folder_n in readdir(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]) ) 
        avgT = Dict(k => 0.0 for k in keys_combo) ; avgY = Dict(k => 0.0 for k in keys_combo)

        count = 0
        avg_n = 0 ; avg_m = 0 ; avg_Y_N = 0

        # ∀ file each line 
        for file in readdir(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]) * "/" * string(folder_n))
            if split(file, ".")[end] == "png"
                continue
            end
            times = [] ; pts = []

            print(fout, file * " & ")

            # write dichotomy result 
            include(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]) * "/" * string(folder_n) * "/" * file)
            print(fout, string(vars) * " & " * string(constr) * " & "  * " ~  & ")

            count += 1
            avg_n += vars ; avg_m += constr ; #avg_Y_N += rootLBS

            for folder_λ in λ_limits
                for m in methods
                    if isfile(work_dir * "/" * string(folder_λ) * "/" * m * "/" * string(folder_n) * "/" * file)
                        include(work_dir * "/" * string(folder_λ) * "/" * m * "/" * string(folder_n) * "/" * file)
                        push!(times, total_times_used); push!(pts, total_nodes)
        
                        avgT[folder_λ * "_" * m] += total_times_used ; avgY[folder_λ * "_" * m] += total_nodes
                    else
                        push!(times, -1); push!(pts, -1)
                    end
                end
            end


            # each line 
            for i=1:length(keys_combo)-1
                if times[i] == -1
                    print(fout, " - & ")
                elseif times[i] >= 3600.0
                    print(fout, " TO & ")
                elseif times[i] == minimum(filter(x -> x > 0 ,times))
                    print(fout, " \\textcolor{blue2}{" * string(times[i]) * "} & ")
                else
                    print(fout, string(times[i]) * " & ")
                end
    
                if pts[i] == -1
                    print(fout, " - & ")
                else
                    print(fout, string(pts[i]) * " & ")
                end
    
            end

            if times[end] == minimum(filter(x -> x > 0 ,times))
                print(fout, " \\textcolor{blue2}{" * string(times[end]) * "} & ")
            elseif times[end] >= 3600.0
                print(fout, " TO & ")
            else
                print(fout, string(times[end]) * " & ") 
            end
            
            println(fout, string(pts[end]) * " \\\\")
    
        end

            
        avg_n = round(Int, avg_n/count) ; avg_m = round(Int, avg_m/count) ; #avg_Y_N = round(Int, avg_Y_N/count) 
        for m in keys_combo
            avgT[m] = round(avgT[m]/count, digits = 2); avgY[m] = round(avgY[m]/count, digits = 2) 
        end

        print(fout, "\\cline{1-16} \\textbf{avg} & \\textbf{" * string(avg_n) * "} & \\textbf{" * string(avg_m) * "} & " )

        for m in keys_combo
            print(fout, "~  & \\textbf{" * string(avgT[m]) * "} & \\textbf{" * string(avgY[m]))
        end

        println(fout, "} \\\\ \\cline{1-16}")


    end

    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout, latex)
    println(fout, "\\caption{cplex cutting LBS non-exhaustive dichotomic concave-convex like algo on instances $instances (\$\\lambda = L_{pred}/K\$ fixed, with lower bound 3, except EPB)  .}")
    println(fout, "\\label{tab:table_predLBS_factor3_$instances }")
    println(fout, "\\end{sidewaystable}")
    close(fout)

end



function comparisonsEPBLambdaLimits(instances::String)
    work_dir = "../../results/" * instances * "/lambda_limit"
    @assert isdir(work_dir) "This directory doesn't exist $work_dir !"

    fout = open(work_dir * "/comparisonEPBLambdaLimitTable.tex", "w")

    latex = raw"""\begin{sidewaystable}[!ht]
    \centering
    \resizebox{\columnwidth}{!}{%
    \begin{tabular}{lccccccccccccccccccccc}
    \toprule
    % line 1 
    \textbf{Instance} & \textbf{n} & \textbf{m} & \textbf{$\mathcal{Y}_N$}
    & \multicolumn{3}{c}{\textbf{$|\lambda| = Inf$}} & \multicolumn{3}{c}{\textbf{$|\lambda| = 2$}} & \multicolumn{3}{c}{\textbf{$|\lambda| = 3$}} & \multicolumn{3}{c}{\textbf{$|\lambda| = 4$}} & \multicolumn{3}{c}{\textbf{$|\lambda| = 6$}} & \multicolumn{3}{c}{\textbf{$|\lambda| = 8$}} \\
    
    %line 2
    \cmidrule(r){5-7} \cmidrule(r){8-10} \cmidrule(r){11-13} \cmidrule(r){14-16} \cmidrule(r){17-19} \cmidrule(r){20-22} 
    ~ & ~ & ~ & ~ & \textbf{Time(s)} & \textbf{\#VB} & \textbf{\#EPB} & \textbf{Time(s)} & \textbf{\#VB} & \textbf{\#EPB} & \textbf{Time(s)} & \textbf{\#VB} & \textbf{\#EPB} & \textbf{Time(s)} & \textbf{\#VB} & \textbf{\#EPB} & \textbf{Time(s)} & \textbf{\#VB} & \textbf{\#EPB} & \textbf{Time(s)} & \textbf{\#VB} & \textbf{\#EPB} \\
    \midrule
    """
    println(fout, latex)

    methods = ["bc_rootRelaxEPB"]

    λ_limits = [] ; 
    for folder_λ in readdir(work_dir)
        if split(folder_λ, ".")[end] == "tex" continue end
        push!(λ_limits, folder_λ)
    end
    keys_combo = [] ; 
    for λ in λ_limits
        for m in methods
            push!(keys_combo, λ * "_" * m)
        end
    end


    # ∀ folder_n  
    for folder_n in readdir(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]) ) 
        avgT = Dict(k => 0.0 for k in keys_combo) ; avgVB = Dict(k => 0.0 for k in keys_combo) ; avgEPB = Dict(k => 0.0 for k in keys_combo)

        count = 0
        avg_n = 0 ; avg_m = 0 ; avg_Y_N = 0

        # ∀ file each line 
        for file in readdir(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]) * "/" * string(folder_n))
            if split(file, ".")[end] == "png"
                continue
            end
            times = [] ; pts = [] ; pts2 = []

            print(fout, file * " & ")

            # write dichotomy result 
            include(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]) * "/" * string(folder_n) * "/" * file)
            print(fout, string(vars) * " & " * string(constr) * " & " * string(size_Y_N) * " & ")

            count += 1
            avg_n += vars ; avg_m += constr ; avg_Y_N += size_Y_N

            for folder_λ in λ_limits
                for m in methods
                    if isfile(work_dir * "/" * string(folder_λ) * "/" * m * "/" * string(folder_n) * "/" * file)
                        include(work_dir * "/" * string(folder_λ) * "/" * m * "/" * string(folder_n) * "/" * file)
                        push!(times, total_times_used); push!(pts, nb_nodes_VB) ; push!(pts2, nb_nodes_EPB)
        
                        avgT[folder_λ * "_" * m] += total_times_used ; avgVB[folder_λ * "_" * m] += nb_nodes_VB ; avgEPB[folder_λ * "_" * m] += nb_nodes_EPB
                    else
                        push!(times, -1); push!(pts, -1) ; push!(pts2, -1)
                    end
                end
            end


            # each line 
            for i=1:length(keys_combo)-1
                if times[i] == -1
                    print(fout, " - & ")
                elseif times[i] >= 3600.0
                    print(fout, " TO & ")
                elseif times[i] == minimum(filter(x -> x > 0 ,times))
                    print(fout, " \\textcolor{blue2}{" * string(times[i]) * "} & ")
                else
                    print(fout, string(times[i]) * " & ")
                end
    
                if pts[i] == -1
                    print(fout, " - & ")
                else
                    print(fout, string(pts[i]) * " & ")
                end

                if pts2[i] == -1
                    print(fout, " - & ")
                else
                    print(fout, string(pts2[i]) * " & ")
                end
    
            end

            if times[end] == minimum(filter(x -> x > 0 ,times))
                print(fout, " \\textcolor{blue2}{" * string(times[end]) * "} & ")
            elseif times[end] >= 3600.0
                print(fout, " TO & ")
            else
                print(fout, string(times[end]) * " & ") 
            end
            
            println(fout, string(pts[end]) * " & " * string(pts2[end]) * " \\\\")
    
        end

            
        avg_n = round(Int, avg_n/count) ; avg_m = round(Int, avg_m/count) ; avg_Y_N = round(Int, avg_Y_N/count) 
        for m in keys_combo
            avgT[m] = round(avgT[m]/count, digits = 2); avgVB[m] = round(avgVB[m]/count, digits = 2) ; avgEPB[m] = round(avgEPB[m]/count, digits = 2) 
        end

        print(fout, "\\cline{1-22} \\textbf{avg} & \\textbf{" * string(avg_n) * "} & \\textbf{" * string(avg_m) * "} & \\textbf{" * string(avg_Y_N) )

        for m in keys_combo
            print(fout, "} & \\textbf{" * string(avgT[m]) * "} & \\textbf{" * string(avgVB[m]) * "} & \\textbf{" * string(avgEPB[m]) )
        end

        println(fout, "} \\\\ \\cline{1-22}")


    end

    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout, latex)
    println(fout, "\\caption{EPB B\\&C(cplex) LBS non-exhaustive dichotomic concave-convex like algo on instances $instances (\$\\lambda\$ fixed) .}")
    println(fout, "\\label{tab:table_lambda_EPB_$instances }")
    println(fout, "\\end{sidewaystable}")
    close(fout)

end






# comparisonsLambdaLimits("momhMKPstu/MOBKP/set3")
# comparisonsLBSFactor2("momhMKPstu/MOBKP/set3")
# comparisonsLBSFactor3("momhMKPstu/MOBKP/set3")

# comparisonsLambdaChordal("momhMKPstu/MOBKP/set3")
# comparisonsLambdaDynamic("momhMKPstu/MOBKP/set3")


comparisonsEPBLambdaLimits("momhMKPstu/MOBKP/set3")


# comparisonPredLBSFactor("momhMKPstu/MOBKP/set3")