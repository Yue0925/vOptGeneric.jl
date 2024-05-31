# ---- Packages to use
using JuMP, CPLEX

include("../../../src/vOptGeneric.jl")
include("parserAP.jl")
using .vOptGeneric


function vOptBOAP(n :: Int64, C1 :: Matrix{Int64}, C2 :: Matrix{Int64},
                method, fname, outputName , limit; step = 0.01
        )

    # ---- setting the model
    model = vModel( CPLEX.Optimizer )
    JuMP.set_silent(model)

    @variable(model, x[1:n, 1:n], Bin)
    @constraint(model, [i = 1:n], sum(x[i, j] for j in 1:n) == 1)
    @constraint(model, [j = 1:n], sum(x[i, j] for i in 1:n) == 1)
    @addobjective(model, Min, sum(C1[i, j] * x[i, j] for i in 1:n for j in 1:n))
    @addobjective(model, Min, sum(C2[i, j] * x[i, j] for i in 1:n for j in 1:n))

    exhaustive = !(limit > 0) 
    if limit == 0 limit = 2^20 end 
    
    if method == :bb
        infos = vSolve( model, method=:bb, verbose=false )
        println(infos)
    elseif method == :bc 
        infos = vSolve( model, method=:bc, verbose=false )
        println(infos)
    elseif method == :dicho 
        start = time()
        vSolve( model, method=:dicho, verbose=false )
        total_time = round(time() - start, digits = 2)
    elseif method==:epsilon 
        start = time()
        vSolve( model, method=:epsilon, step=step, verbose=false )
        total_time = round(time() - start, digits = 2)

    elseif method == :bc_rootRelax 
        infos = vSolve( model, method=:bc_rootRelax, verbose=false, λ_limit = limit )
        println(infos)
    elseif method == :bc_rootRelaxCP 
        infos = vSolve( model, method=:bc_rootRelaxCP, verbose=false, λ_limit = limit )
        println(infos)
    elseif method == :bb_EPB 
        infos = vSolve( model, method=:bb_EPB, verbose=false )
        println(infos)
    elseif method == :bc_rootRelaxEPB 
        infos = vSolve( model, method=:bc_rootRelaxEPB, verbose=false , λ_limit = limit)
        println(infos)
    elseif method == :bc_rootRelaxCPEPB 
        infos = vSolve( model, method=:bc_rootRelaxCPEPB, verbose=false , λ_limit = limit)
        println(infos)
    elseif method == :bc_EPB 
        infos = vSolve( model, method=:bc_EPB, verbose=false )
        println(infos)

    else
        @error "Unknown method parameter $(method) !"
    end

    # ---- Querying the results
    Y_N = getY_N( model )
    println("length Y_N = ", length(Y_N))

    X_E = getX_E( model )
    # println("length X_E = ", length(X_E))



    (method != :dicho && method != :epsilon) ? writeResults(n*n, 2*n, fname, outputName, method, Y_N, X_E; infos) :
        writeResults(n*n, 2*n, fname, outputName, method, Y_N, X_E; total_time)

end

function solve(fname::String, method)
    ap = readAP(fname)
    n = ap.n

    println("\n -----------------------------")
    println(" solving mono $(fname) ... ")
    println(" -----------------------------")

    model = Model( CPLEX.Optimizer ) ; JuMP.set_silent(model)

    @variable(model, x[1:n, 1:n], Bin)
    @constraint(model, [i = 1:n], sum(x[i, j] for j in 1:n) == 1)
    @constraint(model, [j = 1:n], sum(x[i, j] for i in 1:n) == 1)
    @objective(model, Min, sum(ap.C1[i, j] * x[i, j] for i in 1:n for j in 1:n))

    optimize!(model) ; solved_time = round(solve_time(model), digits = 2)

    println("solved time $(solved_time)" )

    status = termination_status(model)
    if status != MOI.OPTIMAL
        @info "mono instance is not feasible"
        return 
    end

    for limit in [0, 2, 3, 6, 8] # 

        println("\n -----------------------------")
        println(" solving $(fname) by $method λ_limit = $limit ... ")
        println(" -----------------------------")

        folder = "../../results/MOAP/"
        if !isdir(folder)
            mkdir(folder)
        end
        folder = "../../results/MOAP/AP/"
        if !isdir(folder)
            mkdir(folder)
        end

        result_folder = folder * "lambda_limit/"
        if !isdir(result_folder)
            mkdir(result_folder)
        end

        result_folder = folder * "lambda_limit/" * string(limit) * "/"
        if !isdir(result_folder)
            mkdir(result_folder)
        end

        result_folder = folder * "lambda_limit/" * string(limit) * "/" * string(method)
        if !isdir(result_folder)
            mkdir(result_folder)
        end

        outputName = result_folder * "/" * ap.name * ".dat"
        if isfile(outputName) continue end #TODO : ignore existed file  

        vOptBOAP(ap.n, ap.C1, ap.C2, Symbol(method), ap.name, outputName, limit)
    end

end

solve(ARGS[1], ARGS[2])

# solve("./AP/AP_p-3_n-5_ins-1.dat", :bb)
