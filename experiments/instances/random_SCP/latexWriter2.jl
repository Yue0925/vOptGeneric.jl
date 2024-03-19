

using PyPlot
import PyPlot; const plt = PyPlot

function comparisonsEPBLambdaLimits(instances::String)
    work_dir = "../../results/" * instances * "/lambda_limit"
    @assert isdir(work_dir) "This directory doesn't exist $work_dir !"

    fout = open(work_dir * "/comparisonEPBLambdaLimitTable.tex", "w")

    latex = raw"""\begin{sidewaystable}[!ht]
    \centering
    \resizebox{\columnwidth}{!}{%
    \begin{tabular}{lccccccccccccccc}
    \toprule
    % line 1 
    \textbf{Instance} & \textbf{n} & \textbf{m} & \textbf{$\mathcal{Y}_N$}
    & \multicolumn{3}{c}{\textbf{$|\lambda| = Inf$}} & \multicolumn{3}{c}{\textbf{$|\lambda| = 2$}} & \multicolumn{3}{c}{\textbf{$|\lambda| = 3$}} & \multicolumn{3}{c}{\textbf{$|\lambda| = 6$}} \\
    
    %line 2
    \cmidrule(r){5-7} \cmidrule(r){8-10} \cmidrule(r){11-13} \cmidrule(r){14-16} 
    ~ & ~ & ~ & ~ & \textbf{Time(s)} & \textbf{\#VB} & \textbf{\#EPB} & \textbf{Time(s)} & \textbf{\#VB} & \textbf{\#EPB} & \textbf{Time(s)} & \textbf{\#VB} & \textbf{\#EPB} & \textbf{Time(s)} & \textbf{\#VB} & \textbf{\#EPB} \\
    \midrule
    """
    println(fout, latex)


    fout2 = open(work_dir * "/comparisonEPBLambdaLimitTable2.tex", "w")

    latex = raw"""\begin{sidewaystable}[!ht]
    \centering
    \resizebox{\columnwidth}{!}{%
    \begin{tabular}{lccccccccccccccc}
    \toprule
    % line 1 
    \textbf{Instance} & \textbf{n} & \textbf{m} & \textbf{$\mathcal{Y}_N$}
    & \multicolumn{3}{c}{\textbf{$|\lambda| = 8$}} & \multicolumn{3}{c}{\textbf{$|\lambda| = 10$}} & \multicolumn{3}{c}{\textbf{$|\lambda| = 13$}} & \multicolumn{3}{c}{\textbf{$|\lambda| = 17$}} \\
    
    %line 2
    \cmidrule(r){5-7} \cmidrule(r){8-10} \cmidrule(r){11-13} \cmidrule(r){14-16} 
    ~ & ~ & ~ & ~ & \textbf{Time(s)} & \textbf{\#VB} & \textbf{\#EPB} & \textbf{Time(s)} & \textbf{\#VB} & \textbf{\#EPB} & \textbf{Time(s)} & \textbf{\#VB} & \textbf{\#EPB} & \textbf{Time(s)} & \textbf{\#VB} & \textbf{\#EPB} \\
    \midrule
    """
    println(fout2, latex)


    methods = ["bc_rootRelaxEPB"]

    λ_limits = [] ; 
    for folder_λ in readdir(work_dir)
        if split(folder_λ, ".")[end] == "tex" continue end
        push!(λ_limits, parse(Int64 , folder_λ ) )
    end
    sort!(λ_limits)
    half = round(Int64, length(λ_limits)/2 )

    keys_combo = [] ; 
    for λ in λ_limits
        for m in methods
            push!(keys_combo, string(λ) * "_" * m)
        end
    end

    # ∀ file each line 
    for file in readdir(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]))
        if split(file, ".")[end] == "png"
            continue
        end

        name_seg = split(file, "_")
        for s in name_seg[1:end-1]
            print(fout, s * "\\_") ; print(fout2, s * "\\_")
        end

        print(fout, name_seg[end] * " & ") ; print(fout2, name_seg[end] * " & ")
        times = [] ; pts = [] ; pts2 = []

        # write dichotomy result 
        include(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]) * "/" * file)
        print(fout, string(vars) * " & " * string(constr) * " & " * string(size_Y_N) * " & ")
        print(fout2, string(vars) * " & " * string(constr) * " & " * string(size_Y_N) * " & ")

        for folder_λ in λ_limits
            for m in methods
                if isfile(work_dir * "/" * string(folder_λ) * "/" * m * "/" * file)
                    include(work_dir * "/" * string(folder_λ) * "/" * m * "/" * file)
                    push!(times, total_times_used); push!(pts, nb_nodes_VB) ; push!(pts2, nb_nodes_EPB)
    
                else
                    push!(times, -1); push!(pts, -1) ; push!(pts2, -1)
                end
            end
        end


        # each line 
        for i=1:(round(Int64, length(keys_combo)/2) -1 )
            if times[i] == -1
                print(fout, " - & ")
            elseif times[i] >= 3600.0
                print(fout, " TO & ")
            elseif times[i] == minimum(filter(x -> x >= 0.0 ,times))
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

        if times[round(Int64, length(keys_combo)/2)] == minimum(filter(x -> x >= 0.0 ,times))
            print(fout, " \\textcolor{blue2}{" * string(times[round(Int64, length(keys_combo)/2)]) * "} & ")
        elseif times[round(Int64, length(keys_combo)/2)] >= 3600.0
            print(fout, " TO & ")
        else
            print(fout, string(times[round(Int64, length(keys_combo)/2)]) * " & ") 
        end
        
        println(fout, string(pts[round(Int64, length(keys_combo)/2)]) * " & " * string(pts2[round(Int64, length(keys_combo)/2)]) * " \\\\")



        # each line 
        for i=(round(Int64, length(keys_combo)/2) +1 ):length(keys_combo)-1
            if times[i] == -1
                print(fout2, " - & ")
            elseif times[i] >= 3600.0
                print(fout2, " TO & ")
            elseif times[i] == minimum(filter(x -> x >= 0.0 ,times))
                print(fout2, " \\textcolor{blue2}{" * string(times[i]) * "} & ")
            else
                print(fout2, string(times[i]) * " & ")
            end

            if pts[i] == -1
                print(fout2, " - & ")
            else
                print(fout2, string(pts[i]) * " & ")
            end

            if pts2[i] == -1
                print(fout2, " - & ")
            else
                print(fout2, string(pts2[i]) * " & ")
            end

        end

        if times[end] == minimum(filter(x -> x >= 0.0 ,times))
            print(fout2, " \\textcolor{blue2}{" * string(times[end]) * "} & ")
        elseif times[end] >= 3600.0
            print(fout2, " TO & ")
        else
            print(fout2, string(times[end]) * " & ") 
        end
        
        println(fout2, string(pts[end]) * " & " * string(pts2[end]) * " \\\\")

    end


    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout, latex)
    println(fout, "\\caption{\\textbf{EPB B\\&C(cplex) }LBS non-exhaustive dichotomic concave-convex like algo on instances $instances .}")
    println(fout, "\\label{tab:table_lambda_EPB_$instances }")
    println(fout, "\\end{sidewaystable}")
    close(fout)


    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout2, latex)
    println(fout2, "\\caption{\\textbf{EPB B\\&C(cplex) }LBS non-exhaustive dichotomic concave-convex like algo on instances $instances .}")
    println(fout2, "\\label{tab:table2_lambda_EPB_$instances }")
    println(fout2, "\\end{sidewaystable}")
    close(fout2)
end


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
    & \multicolumn{4}{c}{\textbf{$|\lambda| = Inf$}} & \multicolumn{4}{c}{\textbf{$|\lambda| = 2$}} & \multicolumn{4}{c}{\textbf{$|\lambda| = 3$}} & \multicolumn{4}{c}{\textbf{$|\lambda| = 6$}} \\
    
    %line 2 
    \cmidrule(r){5-8} \cmidrule(r){9-12} \cmidrule(r){13-16} \cmidrule(r){17-20}
    ~ & ~ & ~ & ~ & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} \\
    
    %line 3
    \cmidrule(r){5-6} \cmidrule(r){7-8} \cmidrule(r){9-10} \cmidrule(r){11-12} \cmidrule(r){13-14} \cmidrule(r){15-16} \cmidrule(r){17-18} \cmidrule(r){19-20}
    ~ & ~ & ~ & ~ & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} \\
    \midrule
    """
    println(fout, latex)

    fout2 = open(work_dir * "/comparisonLambdaLimitTable2.tex", "w")

    latex = raw"""\begin{sidewaystable}[!ht]
    \centering
    \resizebox{\columnwidth}{!}{%
    \begin{tabular}{lccccccccccccccccccc}
    \toprule
    % line 1 
    \textbf{Instance} & \textbf{n} & \textbf{m} & \textbf{$\mathcal{Y}_N$}
    & \multicolumn{4}{c}{\textbf{$|\lambda| = 8$}} & \multicolumn{4}{c}{\textbf{$|\lambda| = 10$}} & \multicolumn{4}{c}{\textbf{$|\lambda| = 13$}} & \multicolumn{4}{c}{\textbf{$|\lambda| = 17$}} \\
    
    %line 2 
    \cmidrule(r){5-8} \cmidrule(r){9-12} \cmidrule(r){13-16} \cmidrule(r){17-20}
    ~ & ~ & ~ & ~ & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} \\
    
    %line 3
    \cmidrule(r){5-6} \cmidrule(r){7-8} \cmidrule(r){9-10} \cmidrule(r){11-12} \cmidrule(r){13-14} \cmidrule(r){15-16} \cmidrule(r){17-18} \cmidrule(r){19-20}
    ~ & ~ & ~ & ~ & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} \\
    \midrule
    """
    println(fout2, latex)


    methods = ["bc_rootRelaxEPB"] # "bc_rootRelax", 
    # method => n => (λ -> time)                     method => n => (λ -> nodes)
    avg_time = Dict{String , Dict{Int64, Dict{Int64, Float64}}}()
    avg_node = Dict{String , Dict{Int64, Dict{Int64, Float64}}}()
    count_per_n = Dict{Int64, Int64}()

    for m in methods
        avg_time[m] = Dict{Int64, Dict{Int64, Float64}}()
        avg_node[m] = Dict{Int64, Dict{Int64, Float64}}()
    end

    λ_limits = [] ;
    for folder_λ in readdir(work_dir)
        if split(folder_λ, ".")[end] == "tex" continue end
        push!(λ_limits, parse(Int64 , folder_λ ) )
    end
    sort!(λ_limits)
    half = round(Int64, length(λ_limits)/2 )

    keys_combo = [] ; 
    for λ in λ_limits
        for m in methods
            push!(keys_combo, string(λ) * "_" * m)
        end
    end

    # ∀ file each line 
    for file in readdir(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]) )
        if split(file, ".")[end] == "png"
            continue
        end

        name_seg = split(file, "_")
        for s in name_seg[1:end-1]
            print(fout, s * "\\_") ; print(fout2, s * "\\_")
        end

        print(fout, name_seg[end] * " & ") ; print(fout2, name_seg[end] * " & ")
        times = [] ; pts = []

        # write dichotomy result 
        include(work_dir * "/" * string(λ_limits[1]) * "/" * string(methods[1]) * "/" * file)
        print(fout, string(vars) * " & " * string(constr) * " & " * string(size_Y_N) * " & ")
        print(fout2, string(vars) * " & " * string(constr) * " & " * string(size_Y_N) * " & ")

        if !haskey(count_per_n, vars)
            count_per_n[vars] = 0
        end

        for m in methods
            if !haskey(avg_time[m], vars)
                avg_time[m][vars] = Dict{Int64, Float64}()
            end
            if !haskey(avg_node[m], vars)
                avg_node[m][vars] = Dict{Int64, Float64}()
            end
        end

        count_per_n[vars] += 1

        for folder_λ in λ_limits
            for m in methods
                if !haskey(avg_time[m][vars], folder_λ)
                    avg_time[m][vars][folder_λ] = 0.0
                end
                if !haskey(avg_node[m][vars], folder_λ)
                    avg_node[m][vars][folder_λ] = 0.0
                end
            end

            for m in methods

                if isfile(work_dir * "/" * string(folder_λ) * "/" * m * "/" * file)
                    include(work_dir * "/" * string(folder_λ) * "/" * m * "/" * file)
                    push!(times, total_times_used); push!(pts, total_nodes)
    
                    avg_time[m][vars][folder_λ] += total_times_used
                    avg_node[m][vars][folder_λ] += total_nodes

                else
                    push!(times, -1); push!(pts, -1)
                end
            end
        end


        # each line 
        for i=1:(round(Int64, length(keys_combo)/2) -1 )
            if times[i] == -1
                print(fout, " - & ")
            elseif times[i] >= 3600.0
                print(fout, " TO & ")
            elseif times[i] == minimum(filter(x -> x >= 0.0 ,times))
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

        if times[round(Int64, length(keys_combo)/2)] == minimum(filter(x -> x >= 0.0 ,times))
            print(fout, " \\textcolor{blue2}{" * string(times[round(Int64, length(keys_combo)/2)]) * "} & ")
        elseif times[round(Int64, length(keys_combo)/2)] >= 3600.0
            print(fout, " TO & ")
        else
            print(fout, string(times[round(Int64, length(keys_combo)/2)]) * " & ") 
        end
        
        println(fout, string(pts[round(Int64, length(keys_combo)/2)]) * " \\\\")



        # each line 
        for i=(round(Int64, length(keys_combo)/2) +1 ):length(keys_combo)-1
            if times[i] == -1
                print(fout2, " - & ")
            elseif times[i] >= 3600.0
                print(fout2, " TO & ")
            elseif times[i] == minimum(filter(x -> x >= 0.0 ,times))
                print(fout2, " \\textcolor{blue2}{" * string(times[i]) * "} & ")
            else
                print(fout2, string(times[i]) * " & ")
            end

            if pts[i] == -1
                print(fout2, " - & ")
            else
                print(fout2, string(pts[i]) * " & ")
            end

        end

        if times[end] == minimum(filter(x -> x >= 0.0 ,times))
            print(fout2, " \\textcolor{blue2}{" * string(times[end]) * "} & ")
        elseif times[end] >= 3600.0
            print(fout2, " TO & ")
        else
            print(fout2, string(times[end]) * " & ") 
        end
        
        println(fout2, string(pts[end]) * " \\\\")

    end


    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout, latex)
    println(fout, "\\caption{cplex cutting LBS non-exhaustive dichotomic concave-convex like algo on instances $instances .}")
    println(fout, "\\label{tab:table_lambda_limits_$instances }")
    println(fout, "\\end{sidewaystable}")
    close(fout)

    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout2, latex)
    println(fout2, "\\caption{cplex cutting LBS non-exhaustive dichotomic concave-convex like algo on instances $instances .}")
    println(fout2, "\\label{tab:table2_lambda_limits_$instances }")
    println(fout2, "\\end{sidewaystable}")
    close(fout2)

    for m in methods
        for n in keys(count_per_n)
            for λ in λ_limits
                avg_node[m][n][λ] = round(avg_node[m][n][λ]/count_per_n[n], digits = 1)
                avg_time[m][n][λ] = round(avg_time[m][n][λ]/count_per_n[n], digits = 1)
            end
            # println("$m  $n  nodes " , [p[2] for p in sort(collect(avg_node[m][n]), by = x->x[1])] )
            # println("$m  $n  times " , [p[2] for p in sort(collect(avg_time[m][n]), by = x->x[1])] )

        end
    end

    println(" --------------------   ")
    println("λ_limits : ", λ_limits)
    println("avg_node : ", avg_node)
    println("avg_time : ", avg_time)
    println(" --------------------   ")



    # plot for each method, for each n 
    for m in methods
        for n in keys(count_per_n)
    
            fig, ax = plt.subplots()
            ax.plot(λ_limits,  [p[2] for p in sort(collect(avg_time[m][n]), by = x->x[1])] ,
                    color="red", marker="o", linewidth=2.0, linestyle="--"
            )
            ax.set_xlabel("|λ|", fontsize=14)
            ax.set_ylabel("Average time(s)", color="red", fontsize=14)
            ax.set_title(instances * " avg n = $n ($m)", fontsize=14)
        
            ax2=ax.twinx()
            ax2.plot(λ_limits, [p[2] for p in sort(collect(avg_node[m][n]), by = x->x[1])], 
                    color="blue", marker="o", linewidth=2.0, linestyle="--"
            )
            ax2.set_ylabel("Average explored nodes", color="blue", fontsize=14)
            # plt.show()
        
            savefig(work_dir * "/$(m)_$(n)_times_nodes.png")
            plt.close()

        end
    end



end



comparisonsLambdaLimits("SCPrandom")
# comparisonsEPBLambdaLimits("SCPrandom")

