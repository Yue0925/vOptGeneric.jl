# Bi-objective unidimensionnal 01 knapsack problem (biukp)
#
# Example 2 with 20 variables


# ---- Packages to use
using JuMP, CPLEX, LinearAlgebra

include("../../../src/vOptGeneric.jl")
include("../../../src/BO01BB/displayGraphic.jl")
using .vOptGeneric


function writeResults(vars::Int64, constr::Int64, fname::String, outputName::String, method, Y_N, X_E; total_time=nothing, infos=nothing)

    fout = open(outputName, "w")
    println(fout, "vars = $vars ; constr = $constr ")
  
    if method == :bb || method == :bc
        println(fout, infos)
    else
      println(fout, "total_times_used = $total_time")
    end
    println(fout, "size_Y_N = ", length(Y_N))
    println(fout, "Y_N = ", Y_N)
    println(fout)
    println(fout, "size_X_E = ", length(X_E))
  
    close(fout)
  
    displayGraphics(fname,Y_N, outputName)
end


function BOUKP(method, fname; step=0.5)

    # ---- Values of the instance to solve
    p1 = [77,94,71,63,96,82,85,75,72,91,99,63,84,87,79,94,90,60,69,62]
    p2 = [65,90,90,77,95,84,70,94,66,92,74,97,60,60,65,97,93,60,69,74]
    w  = [80,87,68,72,66,77,99,85,70,93,98,72,100,89,67,86,91,79,71,99]
    c  = 900
    size = length(p1)


    # ---- setting the model
    m = vModel( CPLEX.Optimizer ) ; JuMP.set_silent(m)
    @variable( m, x[1:size], Bin )
    @addobjective( m, Max, dot(x, p1) )
    @addobjective( m, Max, dot(x, p2) )
    @constraint( m, dot(x, w) <= c )


    if method == :bb
        infos = vSolve( m, method=:bb, verbose=false )
    elseif method == :bc 
        infos = vSolve( m, method=:bc, verbose=false )
    elseif method == :dicho 
        start = time()
        vSolve( m, method=:dicho, verbose=false )
        total_time = round(time() - start, digits = 2)
    elseif method==:epsilon 
        start = time()
        vSolve( m, method=:epsilon, step=step, verbose=false )
        total_time = round(time() - start, digits = 2)
    else
        @error "Unknown method parameter $(method) !"
    end

    # ---- Querying the results
    Y_N = getY_N( m )

    X_E = getX_E( m )


    (method == :bb || method == :bc) ? 
        writeResults(size, 1, "UKnapsackExample2", fname, method, Y_N, X_E; infos) : 
        writeResults(size, 1, "UKnapsackExample2", fname, method, Y_N, X_E; total_time)

end


function main()
    folder = "../../results/smallExamples"
    for method in [:bc] #  :dicho, :bb
        result_dir = method≠:bb ? folder * "/" * string(method) : folder * "/" * string(method) * "/default"
            if !isdir(result_dir)
                mkdir(result_dir)
            end
            fname = result_dir * "/" * "UKnapsackExample2"

            BOUKP(method, fname) 
    end

    # for step in ["0.1", "0.5", "1"]
    #     run_epsilon_ctr(step)
    # end
end


function run_epsilon_ctr(epsilon::String)
    step = parse(Float64, epsilon)
    folder = "../../results/smallExamples"
    method = :epsilon
    result_dir = folder * "/" * string(method) * "/" * string(method) * "_" * string(step) * "/"
    if !isdir(result_dir)
            mkdir(result_dir)
    end
    fname = result_dir * "UKnapsackExample2"

    BOUKP(method, fname; step=step)  
end

# run_epsilon_ctr(ARGS[1])

main()

# vars = 20 ; constr = 1 
#  # informations of B&B algorithm : 
# total_times_used = 24.85 
# total_nodes = 5245 
# pruned_nodes = 2623 
# GAP = 0.0 
# relaxation_time = 3.5 
# test_dominance_time = 0.05 
# update_incumbent_time = 0.0 
# tree_size = 12.645 


#  # ----------- info about cuts : 
# ite_total = 21267 
# cuts_applied = 8821 
# sp_cuts = 2173 
# mp_cuts = 6648 
# cuts_total = 8821 
# times_calling_dicho = 18.16 
# times_calling_separators = 0.25 
# times_oper_cutPool = 0.65 
# times_total_for_cuts = 20.23 
# times_add_retrieve_cuts = 0.54 



# size_Y_N = 11
# Y_N = [[918.0, 984.0], [924.0, 975.0], [927.0, 972.0], [934.0, 971.0], [935.0, 947.0], [940.0, 943.0], [943.0, 940.0], [948.0, 939.0], [949.0, 915.0], [952.0, 909.0], [955.0, 906.0]]

# size_X_E = 11
