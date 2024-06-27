

using PyPlot
import PyPlot; const plt = PyPlot

function comparisonsLambdaLimits(instances::String)
    work_dir = "../../results/" * instances * "/lambda_limit"
    @assert isdir(work_dir) "This directory doesn't exist $work_dir !"

    fout = open(work_dir * "/comparisonLambdaLimitTable.tex", "w")

    latex = raw"""\begin{sidewaystable}[!ht]
    \centering
    \resizebox{\columnwidth}{!}{%
    \begin{tabular}{lccccccccccccccc}
    \toprule
    % line 1 
    \textbf{Instance} & \textbf{n} & \textbf{m} & \textbf{$\mathcal{Y}_N$}
    & \multicolumn{4}{c}{\textbf{$|\lambda| = Inf$}} & \multicolumn{4}{c}{\textbf{$|\lambda| = 2$}} & \multicolumn{4}{c}{\textbf{$|\lambda| = 3$}} \\
    
    %line 2 
    \cmidrule(r){5-8} \cmidrule(r){9-12} \cmidrule(r){13-16} 
    ~ & ~ & ~ & ~ & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} \\
    
    %line 3
    \cmidrule(r){5-6} \cmidrule(r){7-8} \cmidrule(r){9-10} \cmidrule(r){11-12} \cmidrule(r){13-14} \cmidrule(r){15-16} 
    ~ & ~ & ~ & ~ & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} \\
    \midrule
    """
    println(fout, latex)

    fout2 = open(work_dir * "/comparisonLambdaLimitTable2.tex", "w")

    latex = raw"""\begin{sidewaystable}[!ht]
    \centering
    \resizebox{\columnwidth}{!}{%
    \begin{tabular}{lccccccccccc}
    \toprule
    % line 1 
    \textbf{Instance} & \textbf{n} & \textbf{m} & \textbf{$\mathcal{Y}_N$}
    & \multicolumn{4}{c}{\textbf{$|\lambda| = 6$}} & \multicolumn{4}{c}{\textbf{$|\lambda| = 8$}} \\
    
    %line 2 
    \cmidrule(r){5-8} \cmidrule(r){9-12} 
    ~ & ~ & ~ & ~ & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} \\
    
    %line 3
    \cmidrule(r){5-6} \cmidrule(r){7-8} \cmidrule(r){9-10} \cmidrule(r){11-12}
    ~ & ~ & ~ & ~ & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes}  & \textbf{Time(s)} & \textbf{Nodes} \\
    \midrule
    """
    println(fout2, latex)


    methods = ["bc_rootRelax", "bc_rootRelaxEPB"] #    

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

    # println("λ_limits : ", λ_limits) ; println("keys_combo : ", keys_combo)


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
        for i=1:5
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

        if times[6] == minimum(filter(x -> x >= 0.0 ,times))
            print(fout, " \\textcolor{blue2}{" * string(times[6]) * "} & ")
        elseif times[6] >= 3600.0
            print(fout, " TO & ")
        else
            print(fout, string(times[6]) * " & ") 
        end
        
        println(fout, string(pts[6]) * " \\\\")



        # each line 
        for i=7:9
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
    println(fout, "\\caption{cplex cutting LBS non-exhaustive dichotomic concave-convex like algo on instances $instances (\$\\lambda\$ fixed except EPBranched nodes) .}")
    println(fout, "% \\label{tab:table_lambda_limits_$instances }")
    println(fout, "\\end{sidewaystable}")
    close(fout)

    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout2, latex)
    println(fout2, "\\caption{cplex cutting LBS non-exhaustive dichotomic concave-convex like algo on instances $instances (\$\\lambda\$ fixed except EPBranched nodes) .}")
    println(fout2, "% \\label{tab:table2_lambda_limits_$instances }")
    println(fout2, "\\end{sidewaystable}")
    close(fout2)



    for m in methods
        for n in keys(count_per_n)
            for λ in λ_limits
                avg_node[m][n][λ] = round(avg_node[m][n][λ]/count_per_n[n], digits = 2)
                avg_time[m][n][λ] = round(avg_time[m][n][λ]/count_per_n[n], digits = 2)
            end
            println("$m  $n  nodes " ,  avg_node[m][n])  
            println("$m  $n  times " ,  avg_time[m][n]) 

        end
    end

    println("   --------------------   ")
    println("λ_limits : ", λ_limits)
    println("avg_node : ", avg_node)
    println("avg_time : ", avg_time)
    println("   --------------------   ")



    # plot for each method, for each n 
    for m in methods
        for n in keys(count_per_n)
    
            fig, ax = plt.subplots()
            ax.plot(λ_limits,  [p[2] for p in sort(collect(avg_time[m][n]), by = x->x[1])] ,
                    color="red", marker="o", linewidth=2.0, linestyle="--"
            )
            ax.set_xlabel("|λ|", fontsize=14)
            ax.set_ylabel("Average time(s)", color="red", fontsize=14)
            ax.set_title("KP_Forget avg n = $n", fontsize=14)
        
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


comparisonsLambdaLimits("KP_Forget")



# bc_rootRelax  15  nodes Dict(0 => 350.12, 6 => 350.12, 2 => 350.12, 8 => 350.12, 3 => 350.12)
# bc_rootRelax  15  times Dict(0 => 4.61, 6 => 1.61, 2 => 1.62, 8 => 1.6, 3 => 1.62)
# bc_rootRelax  20  nodes Dict(0 => 754.08, 6 => 754.08, 2 => 754.08, 8 => 754.09, 3 => 754.08)
# bc_rootRelax  20  times Dict(0 => 8.19, 6 => 5.17, 2 => 5.19, 8 => 5.17, 3 => 5.18)
# bc_rootRelax  25  nodes Dict(0 => 1692.87, 6 => 1692.87, 2 => 1692.87, 8 => 1692.87, 3 => 1692.87)
# bc_rootRelax  25  times Dict(0 => 18.6, 6 => 15.48, 2 => 15.53, 8 => 15.47, 3 => 15.44)
# bc_rootRelax  10  nodes Dict(0 => 72.63, 6 => 72.63, 2 => 72.63, 8 => 72.63, 3 => 72.63)
# bc_rootRelax  10  times Dict(0 => 2.94, 6 => 0.22, 2 => 0.22, 8 => 0.22, 3 => 0.22)
# bc_rootRelax  35  nodes Dict(0 => 6410.53, 6 => 6410.53, 2 => 6410.53, 8 => 6410.53, 3 => 6410.53)
# bc_rootRelax  35  times Dict(0 => 98.05, 6 => 94.59, 2 => 95.06, 8 => 94.17, 3 => 94.9)
# bc_rootRelax  30  nodes Dict(0 => 3088.27, 6 => 3088.27, 2 => 3088.27, 8 => 3088.27, 3 => 3088.27)
# bc_rootRelax  30  times Dict(0 => 38.85, 6 => 36.1, 2 => 35.79, 8 => 36.24, 3 => 35.79)
# bc_rootRelax  40  nodes Dict(0 => 8781.2, 6 => 8781.2, 2 => 8781.2, 8 => 8781.2, 3 => 8781.2)
# bc_rootRelax  40  times Dict(0 => 162.25, 6 => 158.68, 2 => 159.94, 8 => 157.66, 3 => 159.11)
# bc_rootRelaxEPB  15  nodes Dict(0 => 134.11, 6 => 134.11, 2 => 134.11, 8 => 134.11, 3 => 134.11)
# bc_rootRelaxEPB  15  times Dict(0 => 4.03, 6 => 0.62, 2 => 0.61, 8 => 0.61, 3 => 0.62)
# bc_rootRelaxEPB  20  nodes Dict(0 => 226.37, 6 => 226.37, 2 => 226.37, 8 => 226.37, 3 => 226.37)
# bc_rootRelaxEPB  20  times Dict(0 => 4.74, 6 => 1.31, 2 => 1.31, 8 => 1.3, 3 => 1.32)
# bc_rootRelaxEPB  25  nodes Dict(0 => 671.5, 6 => 671.5, 2 => 671.5, 8 => 671.5, 3 => 671.5)
# bc_rootRelaxEPB  25  times Dict(0 => 7.69, 6 => 4.12, 2 => 4.16, 8 => 4.11, 3 => 4.14)
# bc_rootRelaxEPB  10  nodes Dict(0 => 42.28, 6 => 42.28, 2 => 42.28, 8 => 42.28, 3 => 42.28)
# bc_rootRelaxEPB  10  times Dict(0 => 3.25, 6 => 0.14, 2 => 0.14, 8 => 0.14, 3 => 0.14)
# bc_rootRelaxEPB  35  nodes Dict(0 => 2805.8, 6 => 2805.8, 2 => 2805.83, 8 => 2805.83, 3 => 2805.83)
# bc_rootRelaxEPB  35  times Dict(0 => 24.29, 6 => 20.62, 2 => 20.7, 8 => 20.63, 3 => 20.66)
# bc_rootRelaxEPB  30  nodes Dict(0 => 1395.53, 6 => 1395.53, 2 => 1395.53, 8 => 1395.53, 3 => 1395.53)
# bc_rootRelaxEPB  30  times Dict(0 => 12.75, 6 => 9.16, 2 => 9.21, 8 => 9.16, 3 => 9.18)
# bc_rootRelaxEPB  40  nodes Dict(0 => 4659.57, 6 => 4659.53, 2 => 4659.53, 8 => 4659.57, 3 => 4659.57)
# bc_rootRelaxEPB  40  times Dict(0 => 41.71, 6 => 37.99, 2 => 37.9, 8 => 37.98, 3 => 37.84)
#    --------------------   
# λ_limits : Any[0, 2, 3, 6, 8]
# avg_node : Dict("bc_rootRelax" => Dict(15 => Dict(0 => 350.12, 6 => 350.12, 2 => 350.12, 8 => 350.12, 3 => 350.12), 20 => Dict(0 => 754.08, 6 => 754.08, 2 => 754.08, 8 => 754.09, 3 => 754.08), 25 => Dict(0 => 1692.87, 6 => 1692.87, 2 => 1692.87, 8 => 1692.87, 3 => 1692.87), 10 => Dict(0 => 72.63, 6 => 72.63, 2 => 72.63, 8 => 72.63, 3 => 72.63), 35 => Dict(0 => 6410.53, 6 => 6410.53, 2 => 6410.53, 8 => 6410.53, 3 => 6410.53), 30 => Dict(0 => 3088.27, 6 => 3088.27, 2 => 3088.27, 8 => 3088.27, 3 => 3088.27), 40 => Dict(0 => 8781.2, 6 => 8781.2, 2 => 8781.2, 8 => 8781.2, 3 => 8781.2)), "bc_rootRelaxEPB" => Dict(15 => Dict(0 => 134.11, 6 => 134.11, 2 => 134.11, 8 => 134.11, 3 => 134.11), 20 => Dict(0 => 226.37, 6 => 226.37, 2 => 226.37, 8 => 226.37, 3 => 226.37), 25 => Dict(0 => 671.5, 6 => 671.5, 2 => 671.5, 8 => 671.5, 3 => 671.5), 10 => Dict(0 => 42.28, 6 => 42.28, 2 => 42.28, 8 => 42.28, 3 => 42.28), 35 => Dict(0 => 2805.8, 6 => 2805.8, 2 => 2805.83, 8 => 2805.83, 3 => 2805.83), 30 => Dict(0 => 1395.53, 6 => 1395.53, 2 => 1395.53, 8 => 1395.53, 3 => 1395.53), 40 => Dict(0 => 4659.57, 6 => 4659.53, 2 => 4659.53, 8 => 4659.57, 3 => 4659.57)))
# avg_time : Dict("bc_rootRelax" => Dict(15 => Dict(0 => 4.61, 6 => 1.61, 2 => 1.62, 8 => 1.6, 3 => 1.62), 20 => Dict(0 => 8.19, 6 => 5.17, 2 => 5.19, 8 => 5.17, 3 => 5.18), 25 => Dict(0 => 18.6, 6 => 15.48, 2 => 15.53, 8 => 15.47, 3 => 15.44), 10 => Dict(0 => 2.94, 6 => 0.22, 2 => 0.22, 8 => 0.22, 3 => 0.22), 35 => Dict(0 => 98.05, 6 => 94.59, 2 => 95.06, 8 => 94.17, 3 => 94.9), 30 => Dict(0 => 38.85, 6 => 36.1, 2 => 35.79, 8 => 36.24, 3 => 35.79), 40 => Dict(0 => 162.25, 6 => 158.68, 2 => 159.94, 8 => 157.66, 3 => 159.11)), "bc_rootRelaxEPB" => Dict(15 => Dict(0 => 4.03, 6 => 0.62, 2 => 0.61, 8 => 0.61, 3 => 0.62), 20 => Dict(0 => 4.74, 6 => 1.31, 2 => 1.31, 8 => 1.3, 3 => 1.32), 25 => Dict(0 => 7.69, 6 => 4.12, 2 => 4.16, 8 => 4.11, 3 => 4.14), 10 => Dict(0 => 3.25, 6 => 0.14, 2 => 0.14, 8 => 0.14, 3 => 0.14), 35 => Dict(0 => 24.29, 6 => 20.62, 2 => 20.7, 8 => 20.63, 3 => 20.66), 30 => Dict(0 => 12.75, 6 => 9.16, 2 => 9.21, 8 => 9.16, 3 => 9.18), 40 => Dict(0 => 41.71, 6 => 37.99, 2 => 37.9, 8 => 37.98, 3 => 37.84)))
