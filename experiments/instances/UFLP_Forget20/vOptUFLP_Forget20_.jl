using JuMP, CPLEX

include("../../../src/vOptGeneric.jl")
include("parserForget20.jl")
using .vOptGeneric


function vOptUFLP(forget::Forget,
    method, fname, outputName , limit; step = 0.01
)

    exhaustive = !(limit > 0) 
    if limit == 0 limit = 2^20 end 

    # ---- setting the model
    model = vModel( CPLEX.Optimizer )
    JuMP.set_silent(model)

    n = forget.n ; m = length(forget.m[0])-1 + length(forget.m[1])-1 + length(forget.m[2])-1


    @variable(model, x[1:n], Bin)

    if length(forget.m[0]) > 1
        @constraint(model, [i = 1:size(forget.A0, 1)], forget.A0[i, :]'* x >= forget.b0[i])
    end

    if length(forget.m[1]) > 1
        @constraint(model, [i = 1:size(forget.A1, 1)], forget.A1[i, :]'* x <= forget.b1[i])
    end

    if length(forget.m[2]) > 1
        @constraint(model,  [i = 1:size(forget.A2, 1)], forget.A2[i, :]'* x == forget.b2[i])
    end

    @addobjective(model, Min, forget.c1'* x)
    @addobjective(model, Min, forget.c2'* x)

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
        infos = vSolve( model, method=:bc_rootRelax, verbose=false  , λ_limit = limit )
        println(infos)
    elseif method == :bc_rootRelaxCP 
        infos = vSolve( model, method=:bc_rootRelaxCP, verbose=false  , λ_limit = limit )
        println(infos)
    elseif method == :bb_EPB 
        infos = vSolve( model, method=:bb_EPB, verbose=false )
        println(infos)
    elseif method == :bc_rootRelaxEPB 
        infos = vSolve( model, method=:bc_rootRelaxEPB, verbose=false  , λ_limit = limit)
        println(infos)
    elseif method == :bc_rootRelaxCPEPB 
        infos = vSolve( model, method=:bc_rootRelaxCPEPB, verbose=false  , λ_limit = limit)
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



    (method != :dicho && method != :epsilon) ? writeResults(n, m, fname, outputName, method, Y_N, X_E; infos) :
                    writeResults(n, m, fname, outputName, method, Y_N, X_E; total_time)

end




function solve(fname::String, method)
    forget = readForget(fname)
    n = forget.n ; m = length(forget.m[0])-1 + length(forget.m[1])-1 + length(forget.m[2])-1

    println("\n -----------------------------")
    println(" solving mono $(fname) with n=$n m=$m ... ")
    println(" -----------------------------")

    model = Model( CPLEX.Optimizer ) ; JuMP.set_silent(model)

    @variable(model, x[1:n], Bin)

    if length(forget.m[0]) > 1
        @constraint(model, [i = 1:size(forget.A0, 1)], forget.A0[i, :]'* x >= forget.b0[i])
    end

    if length(forget.m[1]) > 1
        @constraint(model, [i = 1:size(forget.A1, 1)], forget.A1[i, :]'* x <= forget.b1[i])
    end

    if length(forget.m[2]) > 1
        @constraint(model,  [i = 1:size(forget.A2, 1)], forget.A2[i, :]'* x == forget.b2[i])
    end

    @objective(model, Min, forget.c1'* x)

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

        folder = "../../results/UFLP_Forget"
        if !isdir(folder)
            mkdir(folder)
        end

        result_folder = folder * "/lambda_limit/"
        if !isdir(result_folder)
            mkdir(result_folder)
        end

        result_folder = folder * "/lambda_limit/" * string(limit) * "/"
        if !isdir(result_folder)
            mkdir(result_folder)
        end

        result_folder = folder * "/lambda_limit/" * string(limit) * "/" * string(method)
        if !isdir(result_folder)
            mkdir(result_folder)
        end


        outputName = result_folder * "/" * forget.name * ".dat"
        if isfile(outputName) continue end #TODO : ignore existed file  
    
        vOptUFLP(forget, Symbol(method), forget.name, outputName, limit)

    end

end

solve(ARGS[1], ARGS[2])
