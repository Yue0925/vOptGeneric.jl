using PyPlot
import PyPlot; const plt = PyPlot

function comparisonThreeMethods(instances::String)
    work_dir = "../../results/" * instances
    @assert isdir(work_dir) "This directory doesn't exist $work_dir !"

    fout = open(work_dir * "/comparisonTable.tex", "w")

    latex = raw"""\begin{table}[h]
    \centering
    \resizebox{\columnwidth}{!}{%
    \hspace*{-1cm}\begin{tabular}{lcccccccccccc}
    \toprule
    \textbf{Instance} & \textbf{n} & \textbf{m} & \multicolumn{2}{c}{\textbf{Dichotomy}} & \multicolumn{2}{c}{\textbf{$\mathbf{\epsilon}$-constraint}}  & \multicolumn{3}{c}{\textbf{Branch-and-bound}} & \multicolumn{3}{c}{\textbf{Branch-and-cut}}
    \\
    \cmidrule(r){4-5} \cmidrule(r){6-7} \cmidrule(r){8-10} \cmidrule(r){11-13}
    ~ & ~ & ~ & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} & \textbf{$|\mathcal{X}_E|$} & \textbf{Time(s)} & \textbf{$|\mathcal{Y}_N|$} & \textbf{$|\mathcal{X}_E|$} \\
    \midrule
    """
    println(fout, latex)
    methods = ["epsilon", "B&B", "B&C"] ; record_n = []
    record_times = Dict(k => [] for k in methods) ; record_nodes = Dict(k => [] for k in methods[2:end])

    for folder_n in readdir(work_dir * "/dicho") # ∀ filder_n
        count = 0
        avg_n = 0
        avg_m = 0
        avg_dicho_T = 0.0
        avg_dicho_Y = 0.0
        avg_ϵ_T = 0.0
        avg_ϵ_Y = 0.0
        avg_bb_T = 0.0
        avg_bb_Y = 0.0
        avg_bb_X = 0.0
        avg_bb_node = 0.0
        avg_bc_T = 0.0
        avg_bc_Y = 0.0
        avg_bc_X = 0.0
        avg_bc_node = 0.0

        for file in readdir(work_dir * "/dicho/" * string(folder_n) * "/") # ∀ file in dicho
            if split(file, ".")[end] == "png"
                continue
            end

            print(fout, file * " & ")
            times = []
            pts = [] ; X = []

            # write dichotomy result 
            include(work_dir * "/dicho/" * string(folder_n) * "/" * file)
            print(fout, string(vars) * " & " * string(constr) * " & ")
            # print(fout, string(total_times_used)* " & " * string(size_Y_N) * " & ")
            push!(times, total_times_used); push!(pts, size_Y_N)

            count += 1
            avg_n += vars
            avg_m += constr
            avg_dicho_T += total_times_used
            avg_dicho_Y += size_Y_N

            # write ϵ-constraint result (ϵ = 0.5 by default)
            if isfile(work_dir * "/epsilon/" * string(folder_n) * "/" * file)
                include(work_dir * "/epsilon/" * string(folder_n) * "/" * file)
                # print(fout, string(total_times_used)* " & " * string(size_Y_N) * " & ")
                push!(times, total_times_used); push!(pts, size_Y_N)

                avg_ϵ_T += total_times_used
                avg_ϵ_Y += size_Y_N
            else
                # print(fout, "- & - & ")
                push!(times, -1); push!(pts, -1)
            end

            # write B&B result 
            if isfile(work_dir * "/bb/" * string(folder_n) * "/" * file)
                include(work_dir * "/bb/" * string(folder_n) * "/" * file)
                # print(fout, string(total_times_used)* " & " * string(size_Y_N) * " & ")
                push!(times, total_times_used); push!(pts, size_Y_N) 
                push!(X, size_X_E)

                avg_bb_T += total_times_used
                avg_bb_Y += size_Y_N
                avg_bb_X += size_X_E
                avg_bb_node += total_nodes
            else
                # print(fout, "- & - & ")
                push!(times, -1); push!(pts, -1)
                push!(X, -1)
            end

            # write B&C result 
            if isfile(work_dir * "/bc/" * string(folder_n) * "/" * file)
                include(work_dir * "/bc/" * string(folder_n) * "/" * file)
                # print(fout, string(total_times_used)* " & " * string(size_Y_N) * " & ")
                push!(times, total_times_used); push!(pts, size_Y_N) 
                push!(X, size_X_E)

                avg_bc_T += total_times_used
                avg_bc_Y += size_Y_N
                avg_bc_X += size_X_E
                avg_bc_node += total_nodes
            else
                # print(fout, "- & - & ")
                push!(times, -1); push!(pts, -1)
                push!(X, -1)
            end

            for i=1:4
                if times[i] == -1
                    print(fout, " - & ")
                elseif times[i] == minimum(times)
                    print(fout, " \\textcolor{blue2}{" * string(times[i]) * "} & ")
                else
                    print(fout, string(times[i]) * " & ")
                end
    
                if pts[i] == -1
                    print(fout, " - & ")
                elseif pts[i] == maximum(pts)
                    print(fout, " \\textbf{" * string(pts[i]) * "} & ")
                else
                    print(fout, string(pts[i]) * " & ")
                end

                if i==3
                    X[1]==-1 ? print(fout, " - & ") : print(fout, string(X[1]) * " & ")
                elseif i==4
                    X[2]==-1 ? println(fout, " - \\\\") : println(fout, string(X[2]) * " \\\\")
                end
            end
    
        end

        avg_n = round(Int, avg_n/count)
        avg_m = round(Int, avg_m/count)
        avg_dicho_T = round(avg_dicho_T/count, digits = 2)
        avg_dicho_Y = round(avg_dicho_Y/count, digits = 2)
        avg_ϵ_T = round(avg_ϵ_T/count, digits = 2)
        avg_ϵ_Y = round(avg_ϵ_Y/count, digits = 2)
        avg_bb_T = round(avg_bb_T/count, digits = 2)
        avg_bb_Y = round(avg_bb_Y/count, digits = 2)
        avg_bb_X = round(avg_bb_X/count, digits = 2)
        avg_bb_node = round(avg_bb_node/count, digits = 2)
        avg_bc_T = round(avg_bc_T/count, digits = 2)
        avg_bc_Y = round(avg_bc_Y/count, digits = 2)
        avg_bc_X = round(avg_bc_X/count, digits = 2)
        avg_bc_node = round(avg_bc_node/count, digits = 2)

        append!(record_n, avg_n)
        append!(record_times["epsilon"], avg_ϵ_T) ; append!(record_times["B&B"], avg_bb_T) ; append!(record_times["B&C"], avg_bc_T) 
        append!(record_nodes["B&B"], avg_bb_node) ; append!(record_nodes["B&C"], avg_bc_node)


        println(fout, "\\cline{1-13} \\textbf{avg} & \\textbf{" * string(avg_n) * "} & \\textbf{" * string(avg_m) * "} & \\textbf{" *
            string(avg_dicho_T) * "} & \\textbf{" * string(avg_dicho_Y) * "} & \\textbf{" * string(avg_ϵ_T) * "} & \\textbf{" *
            string(avg_ϵ_Y) * "} & \\textbf{" * string(avg_bb_T) * "} & \\textbf{" * string(avg_bb_Y) * "} & \\textbf{" * string(avg_bb_X) *
            "} & \\textbf{" * string(avg_bc_T) * "} & \\textbf{" * string(avg_bc_Y) * "} & \\textbf{" * string(avg_bc_X) * "} \\\\ \\cline{1-13}")

    end

    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout, latex)
    println(fout, "\\caption{Comparison of the different algorithms performances for instances $instances .}")
    println(fout, "\\label{tab:table_compare_$instances }")
    println(fout, "\\end{table}")
    close(fout)

    # labels = [10, 20, 30, 40] ; loc = [0, 1, 2, 3] ; width = 0.3 # the width of the bars

    # plt.bar(loc .- width, record_times["epsilon"], width, label="epsilon")
    # plt.bar(loc, record_times["B&B"], width, label="B&B")
    # plt.bar(loc .+ width, record_times["B&C"], width, label="B&C")
    # plt.xticks(loc, labels)
    # plt.xlabel("Number of variables")
    # plt.ylabel("Computation time(s)", fontsize=14)
    # plt.legend(methods)


    # pos = loc .+ (width/2)
    # for i =1:4
    #     plt.text(x = pos[i]-0.6 , y = record_times["epsilon"][i]+0.5, s = record_times["epsilon"][i], size = 7)
    # end

    # pos = loc .+ (width) .+ (width/2)
    # for i =1:4
    #     plt.text(x = pos[i]-0.6 , y = record_times["B&B"][i]+0.5, s = record_times["B&B"][i], size = 7)
    # end

    # pos = loc .+ (2 * width) .+ (width/2)
    # for i =1:4
    #     plt.text(x = pos[i]-0.6 , y = record_times["B&C"][i]+0.5, s = record_times["B&C"][i], size = 7)
    # end
    # title("Influence of instance size on different algorithms' performance", fontsize=14)
    # savefig(work_dir * "/comparisonTable.png")
    # plt.close()

    # # ---------------
    # width = 0.4
    # plt.bar(loc .- width/2, record_nodes["B&B"], width, label="B&B", color="darkorange")
    # plt.bar(loc .+ width/2, record_nodes["B&C"], width, label="B&C", color="forestgreen")
    # plt.xticks(loc, labels)
    # plt.xlabel("Number of variables")
    # plt.ylabel("Number of explored nodes", fontsize=14)
    # plt.legend(methods[2:end])

    # pos = loc .+ (width)/2# .+ (width/2)
    # for i =1:4
    #     plt.text(x = pos[i]-0.6 , y = record_nodes["B&B"][i]+0.5, s = record_nodes["B&B"][i], size = 7)
    # end

    # pos = loc .+ (3 *width)/2# .+ (width/2)
    # for i =1:4
    #     plt.text(x = pos[i]-0.6 , y = record_nodes["B&C"][i]+0.5, s = record_nodes["B&C"][i], size = 7)
    # end
    # title("Influence of instance size on different algorithms' tree size", fontsize=14)
    # savefig(work_dir * "/comparisonNodes.png")
    # plt.close()
end

# todo : pic ... BB and EPB  
function detailedMOBB_perform(instances::String)
    dir = "../../results/" * instances
    @assert isdir(dir) "This directory doesn't exist $dir !"

    labelT = ["BOLP relaxation", "Dominance test", "UBS update"] ; record_n = []
    labelNode = ["total", "pruned"]
    record_times = Dict(k => [] for k in labelT) ; record_nodes = Dict(k => [] for k in labelNode)

    for folder_n in readdir(dir * "/bb_EPB/")
        count = 0
        avg_n = 0
        avg_m = 0
        avg_totalT = 0.0
        avg_relaxT = 0.0
        avg_dominanceT = 0.0
        avg_incumbentT = 0.0
        avg_totalN= 0
        avg_prunedN = 0
        avg_treeSize = 0.0
        avg_YN = 0
        avg_XE = 0

        for file in readdir(dir * "/bb_EPB/" * string(folder_n) * "/")
            if split(file, ".")[end] == "png"
                continue
            end
    
            times = []
            pts = []
            include(dir * "/bb_EPB/" * string(folder_n) * "/" * file)

            count += 1
            avg_n += vars
            avg_m += constr
            avg_totalT += total_times_used
            avg_relaxT += relaxation_time
            avg_dominanceT += test_dominance_time
            avg_incumbentT += update_incumbent_time
            avg_totalN += total_nodes
            avg_prunedN += pruned_nodes
            avg_treeSize += tree_size
            avg_YN += size_Y_N
            avg_XE += size_X_E
        end

        avg_n = round(Int, avg_n/count)
        avg_m = round(Int, avg_m/count)
        avg_totalT = round(avg_totalT/count, digits = 2)
        avg_relaxT = round(avg_relaxT/count, digits = 2)
        avg_dominanceT = round(avg_dominanceT/count, digits = 2)
        avg_incumbentT = round(avg_incumbentT/count, digits = 2)
        avg_totalN = round(avg_totalN/count, digits = 2)
        avg_prunedN = round(avg_prunedN/count, digits = 2)
        avg_treeSize = round(avg_treeSize/count, digits = 2)
        avg_YN = round(avg_YN/count, digits = 2)
        avg_XE = round(avg_XE/count, digits = 2)

        append!(record_n, avg_n)
        append!(record_times["BOLP relaxation"], avg_relaxT) ; append!(record_times["Dominance test"], avg_dominanceT) ; append!(record_times["UBS update"], avg_incumbentT)
        append!(record_nodes["total"], avg_totalN) ; append!(record_nodes["pruned"], avg_prunedN)
    end

    labels = [10, 20, 30, 40] ; loc = [0, 1, 2, 3] ; width = 0.3 # the width of the bars

    plt.bar(loc .- width, record_times["BOLP relaxation"], width, label="BOLP relaxation")
    plt.bar(loc, record_times["Dominance test"], width, label="Dominance test")
    plt.bar(loc .+ width, record_times["UBS update"], width, label="UBS update")
    plt.xticks(loc, labels)
    plt.xlabel("Number of variables", fontsize=12)
    plt.ylabel("Computation time(s)", fontsize=12)
    plt.legend(labelT)

    pos = loc .+ (width/2)
    for i =1:4
        plt.text(x = pos[i]-0.6 , y = record_times["BOLP relaxation"][i]+0.5, s = record_times["BOLP relaxation"][i], size = 7)
    end

    pos = loc .+ (width) .+ (width/2)
    for i =1:4
        plt.text(x = pos[i]-0.6 , y = record_times["Dominance test"][i]+0.5, s = record_times["Dominance test"][i], size = 7)
    end

    pos = loc .+ (2 * width) .+ (width/2)
    for i =1:4
        plt.text(x = pos[i]-0.6 , y = record_times["UBS update"][i]+0.5, s = record_times["UBS update"][i], size = 7)
    end
    title("The influence of instance size on EPB BO01B&B", fontsize=10)
    savefig(dir * "/EPBBBperformTimes.png")
    plt.close()


    width = 0.4
    plt.bar(loc .- width/2, record_nodes["total"], width, label="total")
    plt.bar(loc .+ width/2, record_nodes["pruned"], width, label="pruned")
    plt.xticks(loc, labels)
    plt.xlabel("Number of variables", fontsize=12)
    plt.ylabel("Number of nodes", fontsize=12)
    # plt.yscale("log")
    plt.legend(labelNode)

    pos = loc .+ (width)/2# .+ (width/2)
    for i =1:4
        plt.text(x = pos[i]-0.6 , y = record_nodes["total"][i]+0.5, s = record_nodes["total"][i], size = 7)
    end

    pos = loc .+ (3 *width)/2# .+ (width/2)
    for i =1:4
        plt.text(x = pos[i]-0.6 , y = record_nodes["pruned"][i]+0.5, s = record_nodes["pruned"][i], size = 7)
    end
    title("The influence of instance size on EPB B&B tree", fontsize=10)
    savefig(dir * "/EPBBBperformNodes.png")
    plt.close()
end


# todo : bc, bc_EPB 
function MOBC_perform(instances::String)
    bc = "/bc_EPB"
    dir = "../../results/" * instances * bc
    @assert isdir(dir) "This directory doesn't exist $dir !"

    fout = open(dir * "/MOBC_bc_EPB.tex", "w")

    latex = raw"""\begin{sidewaystable}[h]
    \centering
    \resizebox{\columnwidth}{!}{%
    \begin{tabular}{lccccccccccccccccc}
    \toprule
    \textbf{Instance} & \textbf{n} & \textbf{m} & \multicolumn{2}{c}{\textbf{Nodes}} & \multicolumn{2}{c}{\textbf{CP iterations}} & \multicolumn{3}{c}{\textbf{Cuts applied}} & \textbf{Cuts} & \multicolumn{5}{c}{\textbf{CP Time(s)}} & \textbf{B\&C Time(s)} & \textbf{$|\mathcal{Y}_N|$}
    \\
    \cmidrule(r){4-5} \cmidrule(r){6-7} \cmidrule(r){8-10} \cmidrule(r){12-16}
    ~ & ~ & ~ & \textbf{total} & \textbf{pruned} & \textbf{total} & \textbf{average} & \textbf{total} & \textbf{sp} & \textbf{mp} & ~ & \textbf{total} & \textbf{dichotomy} & \textbf{pool oper} & \textbf{separators} & \textbf{cuts oper} & ~ & ~ \\
    \midrule
    """
    println(fout, latex)

    labelT = ["BOLP", "cutpool", "separator"] ;# record_n = []
    labelNode = ["total", "pruned"]
    record_times = Dict(k => [] for k in labelT) ; record_nodes = Dict(k => [] for k in labelNode)


    for folder_n in readdir(dir)
        if !isdir(dir * "/" * string(folder_n) ) continue end 
        count = 0
        avg_n = 0 ; avg_m = 0
        avg_totalNodes= 0 ; avg_prunedNodes = 0
        avg_totalIte = 0 ; avg_pernodesIte = 0
        avg_totalCuts = 0 ; avg_spCuts = 0 ; avg_mpCuts = 0 ; avg_Cuts = 0
        avg_totalT = 0.0 ; avg_dichoT = 0.0 ; avg_poolT = 0.0 ; avg_sepaT = 0.0 ; avg_cutsT = 0.0
        avg_BCtime = 0.0 ; avg_YN = 0

        for file in readdir(dir * "/" * string(folder_n) * "/")
            if split(file, ".")[end] == "png"
                continue
            end
    
            print(fout, file * " & ")
    
            include(dir * "/" * string(folder_n) * "/" * file)
            print(fout, string(vars) * " & " * string(constr) * " & ")
            print(fout, string(total_nodes) * " & " * string(pruned_nodes) * " & " * string(ite_total) * " & "*
                string( round(ite_total/total_nodes, digits = 2) ) * " & " * string(cuts_applied) * " & " * string(sp_cuts)*
                " & " * string(mp_cuts) * " & " * string(cuts_total) * " & " * string(times_total_for_cuts) * " & " *
                string(times_calling_dicho) * " & " * string(times_oper_cutPool) * " & " *string(times_calling_separators) * " & "*
                string(times_add_retrieve_cuts) * " & " *string(total_times_used) * " & " *string(size_Y_N)
            )

            println(fout, "\\\\")

            count += 1
            avg_n += vars ; avg_m += constr
            avg_totalNodes += total_nodes ; avg_prunedNodes += pruned_nodes
            avg_totalIte += ite_total ; avg_pernodesIte += round(ite_total/total_nodes, digits = 2)
            avg_totalCuts += cuts_applied ; avg_spCuts += sp_cuts ; avg_mpCuts += mp_cuts ; avg_Cuts += cuts_total
            avg_totalT += times_total_for_cuts ; avg_dichoT += times_calling_dicho ; avg_poolT += times_oper_cutPool ; avg_sepaT += times_calling_separators ; avg_cutsT += times_add_retrieve_cuts
            avg_BCtime += total_times_used ; avg_YN += size_Y_N
        end

        avg_n = round(Int, avg_n/count) ; avg_m = round(Int, avg_m/count)
        avg_totalNodes = round(avg_totalNodes/count, digits = 2)
        avg_prunedNodes = round(avg_prunedNodes/count, digits = 2)
        avg_totalIte = round(avg_totalIte/count, digits = 2)
        avg_pernodesIte = round(avg_pernodesIte/count, digits = 2)
        avg_totalCuts = round(avg_totalCuts/count, digits = 2)
        avg_spCuts = round(avg_spCuts/count, digits = 2)
        avg_mpCuts = round(avg_mpCuts/count, digits = 2)
        avg_Cuts = round(avg_Cuts/count, digits = 2)
        avg_dichoT = round(avg_dichoT/count, digits = 2)
        avg_totalT = round(avg_totalT/count, digits = 2)
        avg_poolT = round(avg_poolT/count, digits = 2)
        avg_sepaT = round(avg_sepaT/count, digits = 2)
        avg_cutsT = round(avg_cutsT/count, digits = 2)
        avg_BCtime = round(avg_BCtime/count, digits = 2)
        avg_YN = round(avg_YN/count, digits = 2)

        append!(record_times["BOLP"], avg_dichoT)
        append!(record_times["cutpool"], avg_poolT)
        append!(record_times["separator"], avg_sepaT)
        append!(record_nodes["total"], avg_totalNodes)
        append!(record_nodes["pruned"], avg_prunedNodes)


        println(fout, "\\cline{1-18} \\textbf{avg} & \\textbf{" * string(avg_n) * "} & \\textbf{" * string(avg_m) * "} & \\textbf{" *
            string(avg_totalNodes) * "} & \\textbf{" * string(avg_prunedNodes) * "} & \\textbf{" * string(avg_totalIte) * "} & \\textbf{" *
            string(avg_pernodesIte) * "} & \\textbf{" * string(avg_totalCuts) * "} & \\textbf{" * string(avg_spCuts) * "} & \\textbf{" *
            string(avg_mpCuts) * "} & \\textbf{" * string(avg_Cuts) * "} & \\textbf{" * string(avg_dichoT) *"} & \\textbf{" *
            string(avg_totalT) * "} & \\textbf{" * string(avg_poolT) * "} & \\textbf{" * string(avg_sepaT) * "} & \\textbf{" *
            string(avg_cutsT) * "} & \\textbf{" * string(avg_BCtime) * "} & \\textbf{" * string(avg_YN) * "} " * "\\\\ \\cline{1-18}")
    end


    latex = raw"""\bottomrule
    \end{tabular}%
    }%
    \caption{.}
    \label{tab:table_bc_EPB}
    \end{sidewaystable}
    """
    println(fout, latex)
    close(fout)

    # ----------------------------------------------------------------
    # ------------------ pics of times -------------------------------
    # ----------------------------------------------------------------
    # labels = [10, 20, 30, 40] ; loc = [0, 1, 2, 3] ; width = 0.3 # the width of the bars

    # plt.bar(loc .- width, record_times["BOLP"], width, label="BOLP")
    # plt.bar(loc, record_times["cutpool"], width, label="cutpool")
    # plt.bar(loc .+ width, record_times["separator"], width, label="separator")
    # plt.xticks(loc, labels)
    # plt.xlabel("Number of variables")
    # plt.ylabel("Computation time(s)", fontsize=14)
    # plt.legend(labelT)


    # pos = loc .+ (width/2)
    # for i =1:4
    #     plt.text(x = pos[i]-0.6 , y = record_times["BOLP"][i]+0.5, s = record_times["BOLP"][i], size = 7)
    # end

    # pos = loc .+ (width) .+ (width/2)
    # for i =1:4
    #     plt.text(x = pos[i]-0.6 , y = record_times["cutpool"][i]+0.5, s = record_times["cutpool"][i], size = 7)
    # end

    # pos = loc .+ (2 * width) .+ (width/2)
    # for i =1:4
    #     plt.text(x = pos[i]-0.6 , y = record_times["separator"][i]+0.5, s = record_times["separator"][i], size = 7)
    # end
    # title("Influence of instance size on EPB BOB\\&C", fontsize=12)
    # savefig(dir * "/EPBBCperformTime.png")
    # plt.close()

    # # ---------------
    # width = 0.4
    # plt.bar(loc .- width/2, record_nodes["total"], width, label="total", color="darkorange")
    # plt.bar(loc .+ width/2, record_nodes["pruned"], width, label="pruned", color="forestgreen")
    # plt.xticks(loc, labels)
    # plt.xlabel("Number of variables")
    # plt.ylabel("Number of explored nodes", fontsize=14)
    # plt.legend(labelNode)

    # pos = loc .+ (width)/2# .+ (width/2)
    # for i =1:4
    #     plt.text(x = pos[i]-0.6 , y = record_nodes["total"][i]+0.5, s = record_nodes["total"][i], size = 7)
    # end

    # pos = loc .+ (3 *width)/2# .+ (width/2)
    # for i =1:4
    #     plt.text(x = pos[i]-0.6 , y = record_nodes["pruned"][i]+0.5, s = record_nodes["pruned"][i], size = 7)
    # end
    # title("Influence of instance size on EPB BOB\\&C tree", fontsize=12)
    # savefig(dir * "/EPBBCperformNodes.png")
    # plt.close()


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
    methods = ["epsilon", "bb", "bc", "bc_rootRelax", "bb_EPB", "bc_EPB", "bc_rootRelaxEPB", "bc_rootRelaxCP", "bc_rootRelaxCPEPB"] ; record_n = []
    record_times = Dict(k => [] for k in methods) ; record_nodes = Dict(k => [] for k in methods[2:end])

    # ∀ filder_n
    for folder_n in readdir(work_dir * "/epsilon") 
        count = 0
        avg_n = 0 ; avg_m = 0
        avgT = Dict(k => 0.0 for k in methods) ; avgY = Dict(k => 0.0 for k in methods)

        # ∀ file in dicho
        for file in readdir(work_dir * "/epsilon/" * string(folder_n) * "/")
            if split(file, ".")[end] == "png"
                continue
            end

            print(fout, file * " & ")
            times = [] ; pts = []

            # write dichotomy result 
            include(work_dir * "/epsilon/" * string(folder_n) * "/" * file)
            print(fout, string(vars) * " & " * string(constr) * " & ")

            count += 1
            avg_n += vars ; avg_m += constr

            # ∀ method 
            for m in methods
                if isfile(work_dir * "/" * m * "/" * string(folder_n) * "/" * file)
                    include(work_dir * "/" * m * "/" * string(folder_n) * "/" * file)
                    push!(times, total_times_used); push!(pts, size_Y_N)
    
                    avgT[m] += total_times_used ; avgY[m] += size_Y_N
                else
                    push!(times, -1); push!(pts, -1)
                end
            end

            # ------------------
            for i=1:length(methods)-1
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

        avg_n = round(Int, avg_n/count) ; avg_m = round(Int, avg_m/count)
        for m in methods
            avgT[m] = round(avgT[m]/count, digits = 2); avgY[m] = round(avgY[m]/count, digits = 2) 
        end

        print(fout, "\\cline{1-21} \\textbf{avg} & \\textbf{" * string(avg_n) * "} & \\textbf{" * string(avg_m) )

        for m in methods
            print(fout, "} & \\textbf{" * string(avgT[m]) * "} & \\textbf{" * string(avgY[m]))
        end

        println(fout, "} \\\\ \\cline{1-21}")
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
    record_times = Dict(k => [] for k in methods) ; record_nodes = Dict(k => [] for k in methods[2:end])

    # ∀ filder_n
    for folder_n in readdir(work_dir * "/epsilon") 
        count = 0
        avg_n = 0 ; avg_m = 0
        avgT = Dict(k => 0.0 for k in methods) ; avgY = Dict(k => 0.0 for k in methods) ; avgTO = Dict(k => 0 for k in methods)
        countY_N = 0

        # ∀ file in dicho
        for file in readdir(work_dir * "/epsilon/" * string(folder_n) * "/")
            if split(file, ".")[end] == "png"
                continue
            end

            print(fout, file * " & ")
            times = [] ; pts = []

            # write dichotomy result 
            include(work_dir * "/epsilon/" * string(folder_n) * "/" * file)
            print(fout, string(vars) * " & " * string(constr) * " & ")

            count += 1
            avg_n += vars ; avg_m += constr
            countY_N += size_Y_N

            # ∀ method 
            for m in methods
                if isfile(work_dir * "/" * m * "/" * string(folder_n) * "/" * file)
                    include(work_dir * "/" * m * "/" * string(folder_n) * "/" * file)
                    push!(times, total_times_used); avgT[m] += total_times_used
                    total_times_used > 3600.0 ? avgTO[m] += 1 : nothing

                    
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
                    print(fout,  string(times[i]) * " & ")
                else
                    print(fout, string(times[i]) * " & ")
                end
    
                if i==1 continue end 
                if pts[i] == -1
                    print(fout, " - & ")
                else
                    print(fout, string(pts[i]) * " & ")
                end

            end

            if times[end] == minimum(filter(x -> x > 0 ,times))
                print(fout,  string(times[end]) * " & ")
            else
                print(fout, string(times[end]) * " & ") 
            end
            
            println(fout, string(pts[end]) * " & " * string(size_Y_N) * " \\\\")
    
        end

        avg_n = round(Int, avg_n/count) ; avg_m = round(Int, avg_m/count)
        for m in methods
            avgT[m] = round(avgT[m]/count, digits = 2); avgY[m] = round(avgY[m]/count, digits = 2) 
        end

        print(fout, "\\cline{1-9} avg & " * string(avg_n) * " & " * string(avg_m) )

        print(fout, " & " * string(avgT[methods[1]]) )
        for m in methods[2:end]
            print(fout, " & " * string(avgT[m]) * "& " * string(avgY[m]))
        end

        println(fout, " & " * string(round(countY_N/count, digits = 2))* "\\\\ \\cline{1-9}")
        println("n = $folder_n , count = $count TO = $avgTO")

    end

    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout, latex)
    println(fout, "\\caption{Comparison of the B\\&B algorithms performances for instances $instances .}")
    println(fout, "\\label{tab:table_EPSILONvsBBvsEPBBB_$instances }")
    println(fout, "\\end{table}")
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
    methods = ["bc", "bc_EPB"] ; record_n = []
    record_sp = Dict(k => [] for k in methods) ; record_mp = Dict(k => [] for k in methods)

    # ∀ filder_n
    for folder_n in readdir(work_dir * "/epsilon") 
        count = 0
        avg_n = 0 ; avg_m = 0
        avgSP = Dict(k => 0.0 for k in methods) ; avgMP = Dict(k => 0.0 for k in methods)

        # ∀ file in dicho
        for file in readdir(work_dir * "/epsilon/" * string(folder_n) * "/")
            if split(file, ".")[end] == "png"
                continue
            end

            print(fout, file * " & ")
            sp = [] ; mp = []

            # write dichotomy result 
            include(work_dir * "/epsilon/" * string(folder_n) * "/" * file)
            print(fout, string(vars) * " & " * string(constr) * " & ")

            count += 1
            avg_n += vars ; avg_m += constr

            # ∀ method 
            for m in methods
                if isfile(work_dir * "/" * m * "/" * string(folder_n) * "/" * file)
                    include(work_dir * "/" * m * "/" * string(folder_n) * "/" * file)
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

        avg_n = round(Int, avg_n/count) ; avg_m = round(Int, avg_m/count)
        for m in methods
            avgSP[m] = round(avgSP[m]/count, digits = 2); avgMP[m] = round(avgMP[m]/count, digits = 2) 
        end

        print(fout, "\\cline{1-7} \\textbf{avg} & \\textbf{" * string(avg_n) * "} & \\textbf{" * string(avg_m) )

        for m in methods
            print(fout, "} & \\textbf{" * string(avgSP[m]) * "} & \\textbf{" * string(avgMP[m]))
        end

        println(fout, "} " * "\\\\ \\cline{1-7}")
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
    methods = ["bb", "bc", "bb_EPB", "bc_EPB"] ; record_n = []
    record_times = Dict(k => [] for k in methods) ; record_nodes = Dict(k => [] for k in methods)

    # ∀ filder_n
    for folder_n in readdir(work_dir * "/epsilon") 
        count = 0
        avg_n = 0 ; avg_m = 0
        avgT = Dict(k => 0.0 for k in methods) ; avgY = Dict(k => 0.0 for k in methods) ;  avgTO = Dict(k => 0 for k in methods)

        # ∀ file in dicho
        for file in readdir(work_dir * "/epsilon/" * string(folder_n) * "/")
            if split(file, ".")[end] == "png"
                continue
            end

            print(fout, file * " & ")
            times = [] ; pts = []

            # write dichotomy result 
            include(work_dir * "/epsilon/" * string(folder_n) * "/" * file)
            print(fout, string(vars) * " & " * string(constr) * " & ")

            count += 1
            avg_n += vars ; avg_m += constr

            # ∀ method 
            for m in methods
                if isfile(work_dir * "/" * m * "/" * string(folder_n) * "/" * file)
                    include(work_dir * "/" * m * "/" * string(folder_n) * "/" * file)
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
            else
                print(fout, string(times[end]) * " & ") 
            end
            
            println(fout, string(pts[end]) * " \\\\")
    
        end

        avg_n = round(Int, avg_n/count) ; avg_m = round(Int, avg_m/count)
        for m in methods
            avgT[m] = round(avgT[m]/count, digits = 2); avgY[m] = round(avgY[m]/count, digits = 2) 
        end

        print(fout, "\\cline{1-11} \\textbf{avg} & \\textbf{" * string(avg_n) * "} & \\textbf{" * string(avg_m) )

        for m in methods
            print(fout, "} & \\textbf{" * string(avgT[m]) * "} & \\textbf{" * string(avgY[m]))
        end

        println(fout, "} \\\\ \\cline{1-11}")
        println("n = $folder_n , count = $count TO = $avgTO")

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
    methods = ["bc", "bc_rootRelax", "bc_rootRelaxCP", "bc_EPB", "bc_rootRelaxEPB", "bc_rootRelaxCPEPB"] ; record_n = []
    record_times = Dict(k => [] for k in methods) ; record_nodes = Dict(k => [] for k in methods)

    # ∀ filder_n
    for folder_n in readdir(work_dir * "/epsilon") 
        count = 0
        avg_n = 0 ; avg_m = 0
        avgT = Dict(k => 0.0 for k in methods) ; avgY = Dict(k => 0.0 for k in methods) ;  avgTO = Dict(k => 0 for k in methods)

        # ∀ file in dicho
        for file in readdir(work_dir * "/epsilon/" * string(folder_n) * "/")
            if split(file, ".")[end] == "png"
                continue
            end

            print(fout, file * " & ")
            times = [] ; pts = []

            # write dichotomy result 
            include(work_dir * "/epsilon/" * string(folder_n) * "/" * file)
            print(fout, string(vars) * " & " * string(constr) * " & ")

            count += 1
            avg_n += vars ; avg_m += constr

            # ∀ method 
            for m in methods
                if isfile(work_dir * "/" * m * "/" * string(folder_n) * "/" * file)
                    include(work_dir * "/" * m * "/" * string(folder_n) * "/" * file)
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
            else
                print(fout, string(times[end]) * " & ") 
            end
            
            println(fout, string(pts[end]) * " \\\\")
    
        end

        avg_n = round(Int, avg_n/count) ; avg_m = round(Int, avg_m/count)
        for m in methods
            avgT[m] = round(avgT[m]/count, digits = 2); avgY[m] = round(avgY[m]/count, digits = 2) 
        end

        print(fout, "\\cline{1-15} \\textbf{avg} & \\textbf{" * string(avg_n) * "} & \\textbf{" * string(avg_m) )

        for m in methods
            print(fout, "} & \\textbf{" * string(avgT[m]) * "} & \\textbf{" * string(avgY[m]))
        end

        println(fout, "} \\\\ \\cline{1-15}")
        println("n = $folder_n , count = $count TO = $avgTO")

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



function comparisons_tri(instances::String)
    work_dir = "../../results/" * instances
    @assert isdir(work_dir) "This directory doesn't exist $work_dir !"

    fout = open(work_dir * "/comparisonTable2.tex", "w")

    latex = raw"""\begin{table}[!ht]
    \centering
    \resizebox{\columnwidth}{!}{%
    \begin{tabular}{lcccccccc}
    \toprule
    \textbf{Instance} & \textbf{n} & \textbf{m} & \textbf{$|\mathcal{Y}_N|$} & \textbf{$\epsilon$-constraint}  & \multicolumn{2}{c}{\textbf{B\&B}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}}
    \\
    \cmidrule(r){6-7} \cmidrule(r){8-9} 
    ~ & ~ & ~ & ~ & ~ & \textbf{Time(s)} &\textbf{Nodes} & \textbf{Time(s)} &\textbf{Nodes} \\
    \midrule
    """
    println(fout, latex)
    methods = ["epsilon", "bb", "bc_rootRelaxEPB"] 

    # ∀ filder_n
    for folder_n in readdir(work_dir * "/epsilon") 
        count = 0
        avg_n = 0 ; avg_m = 0
        avgT = Dict(k => 0.0 for k in methods) ; avgY = Dict(k => 0.0 for k in methods)

        # ∀ file in dicho
        for file in readdir(work_dir * "/epsilon/" * string(folder_n) * "/")
            if split(file, ".")[end] == "png"
                continue
            end

            print(fout, file * " & ")
            times = [] ; pts = []

            # write dichotomy result 
            include(work_dir * "/epsilon/" * string(folder_n) * "/" * file)
            print(fout, string(vars) * " & " * string(constr) * " & ")

            count += 1
            avg_n += vars ; avg_m += constr

            # ∀ method 
            m = methods[1]
            if isfile(work_dir * "/" * m * "/" * string(folder_n) * "/" * file)
                include(work_dir * "/" * m * "/" * string(folder_n) * "/" * file)
                push!(times, total_times_used); push!(pts, size_Y_N)

                avgT[m] += total_times_used ; avgY[m] += size_Y_N
            else
                push!(times, -1); push!(pts, -1)
            end
            for m in methods[2:end]
                if isfile(work_dir * "/" * m * "/" * string(folder_n) * "/" * file)
                    include(work_dir * "/" * m * "/" * string(folder_n) * "/" * file)
                    push!(times, total_times_used); push!(pts, total_nodes)
    
                    avgT[m] += total_times_used ; avgY[m] += total_nodes
                else
                    push!(times, -1); push!(pts, -1)
                end
            end

            # ------------------
            if pts[1] == -1
                print(fout, " - & ")
            else
                print(fout, string(pts[1]) * " & ")
            end
            if times[1] == -1
                print(fout, " - & ")
            else
                print(fout, string(times[1]) * " & ")
            end

            for i=2:length(methods)-1
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

        print(fout, "\\hline \\textbf{avg} & " * string(avg_n) * " & " * string(avg_m) )

        m = methods[1]
        print(fout, " & " * string(avgY[m]) * " & " * string(avgT[m]))
        for m in methods[2:end]
            print(fout, " & " * string(avgT[m]) * " & " * string(avgY[m]))
        end

        println(fout, " " * "\\\\ \\hline")
    end

    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout, latex)
    println(fout, "\\caption{ instances $instances .}")
    println(fout, "\\label{tab:table_compareBB_$instances }")
    println(fout, "\\end{table}")
    close(fout)

end


function comparisonsTree(instances::String)
    work_dir = "../../results/" * instances
    @assert isdir(work_dir) "This directory doesn't exist $work_dir !"

    fout = open(work_dir * "/comparisonTreeTable.tex", "w")

    latex = raw"""\begin{sidewaystable}[!ht]
    \centering
    \resizebox{\columnwidth}{!}{%
    \begin{tabular}{lcccccccccccccccccc}
    \toprule
    \textbf{Instance} & \textbf{n} & \textbf{m}  & \multicolumn{2}{c}{\textbf{B\&B}} & \multicolumn{2}{c}{\textbf{B\&C(LP+CP)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{EPB B\&B}} & \multicolumn{2}{c}{\textbf{EPB B\&C(LP+CP)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex)}} & \multicolumn{2}{c}{\textbf{B\&C(cplex+CP)}} & \multicolumn{2}{c}{\textbf{EPB B\&C(cplex+CP)}} \\
    
    \cmidrule(r){4-5} \cmidrule(r){6-7} \cmidrule(r){8-9} \cmidrule(r){10-11} \cmidrule(r){12-13} \cmidrule(r){14-15} \cmidrule(r){16-17} \cmidrule(r){18-19}
    ~ & ~ & ~ & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes} & \textbf{Time(s)} & \textbf{Nodes}  \\
    \midrule
    """
    println(fout, latex)
    methods = ["bb", "bc", "bc_rootRelax", "bb_EPB", "bc_EPB", "bc_rootRelaxEPB", "bc_rootRelaxCP", "bc_rootRelaxCPEPB"] ; record_n = []
    record_times = Dict(k => [] for k in methods) ; record_nodes = Dict(k => [] for k in methods)

    # ∀ filder_n
    for folder_n in readdir(work_dir * "/epsilon") 
        count = 0
        avg_n = 0 ; avg_m = 0
        avgT = Dict(k => 0.0 for k in methods) ; avgY = Dict(k => 0.0 for k in methods)

        # ∀ file in dicho
        for file in readdir(work_dir * "/epsilon/" * string(folder_n) * "/")
            if split(file, ".")[end] == "png"
                continue
            end

            print(fout, file * " & ")
            times = [] ; pts = []

            # write dichotomy result 
            include(work_dir * "/epsilon/" * string(folder_n) * "/" * file)
            print(fout, string(vars) * " & " * string(constr) * " & ")

            count += 1
            avg_n += vars ; avg_m += constr

            # ∀ method 
            for m in methods
                if isfile(work_dir * "/" * m * "/" * string(folder_n) * "/" * file)
                    include(work_dir * "/" * m * "/" * string(folder_n) * "/" * file)
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

        avg_n = round(Int, avg_n/count) ; avg_m = round(Int, avg_m/count)
        for m in methods
            avgT[m] = round(avgT[m]/count, digits = 2); avgY[m] = round(avgY[m]/count, digits = 2) 
        end

        print(fout, "\\cline{1-19} \\textbf{avg} & \\textbf{" * string(avg_n) * "} & \\textbf{" * string(avg_m) )

        for m in methods
            print(fout, "} & \\textbf{" * string(avgT[m]) * "} & \\textbf{" * string(avgY[m]))
        end

        println(fout, "} \\\\ \\cline{1-19}")
    end

    latex = raw"""\bottomrule
    \end{tabular}
    }%"""
    println(fout, latex)
    println(fout, "\\caption{Comparison of the different algorithms B\\&B tree for instances $instances .}")
    println(fout, "\\label{tab:table_compare_tree_$instances }")
    println(fout, "\\end{sidewaystable}")
    close(fout)

end




# -------------------------------------------------
# comparisons("momhMKPstu/MOBKP/set3")
# comparisonsTree("momhMKPstu/MOBKP/set3")
# comparisons_eps_BB_EPB("momhMKPstu/MOBKP/set3")


comparisons5("momhMKPstu/MOBKP/set3")
# comparisons4("momhMKPstu/MOBKP/set3")
# comparisonsCP("momhMKPstu/MOBKP/set3")
# comparisonThreeMethods("momhMKPstu/MOBKP/set3")
# comparisons_tri("momhMKPstu/MOBKP/set3")
# -------------------------------------------------




# detailedMOBB_perform("momhMKPstu/MOBKP/set3")

# comparisonThreeMethods("momhMKPstu/MOBKP/set3")

# MOBC_perform("momhMKPstu/MOBKP/set3")

