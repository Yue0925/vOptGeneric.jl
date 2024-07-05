using JuMP, CPLEX

include("../../../src/vOptGeneric.jl")
include("parserKP_Forget20.jl")
using .vOptGeneric



function vOptKP(kp::ForgetKP,
    method, fname, outputName , heuristic, mixed, mixed2 ; step = 0.01
)

    # ---- setting the model
    model = vModel( CPLEX.Optimizer )
    JuMP.set_silent(model)

    n = kp.n 

    @variable(model, x[1:n], Bin)
    @constraint(model, kp.A'* x <= kp.b)

    @addobjective(model, Max, kp.c1'* x)
    @addobjective(model, Max, kp.c2'* x)

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
        infos = vSolve( model, method=:bc_rootRelax, verbose=false , heuristic = heuristic, mixed = mixed, mixed2 = mixed2  )
        println(infos)
    elseif method == :bc_rootRelaxCP 
        infos = vSolve( model, method=:bc_rootRelaxCP, verbose=false )
        println(infos)
    elseif method == :bb_EPB 
        infos = vSolve( model, method=:bb_EPB, verbose=false )
        println(infos)
    elseif method == :bc_rootRelaxEPB 
        infos = vSolve( model, method=:bc_rootRelaxEPB, verbose=false )
        println(infos)
    elseif method == :bc_rootRelaxCPEPB 
        infos = vSolve( model, method=:bc_rootRelaxCPEPB, verbose=false )
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



    (method != :dicho && method != :epsilon) ? writeResults(n, 1, fname, outputName, method, Y_N, X_E; infos) :
                    writeResults(n, 1, fname, outputName, method, Y_N, X_E; total_time)

end


function solve(fname::String, method)
    kp = readForgetKP(fname)
    n = kp.n

    println("\n -----------------------------")
    println(" solving mono $(fname) with n=$n... ")
    println(" -----------------------------")

    model = Model( CPLEX.Optimizer ) ; JuMP.set_silent(model)

    @variable(model, x[1:n], Bin)
    @constraint(model, kp.A'* x <= kp.b)

    @objective(model, Max, kp.c1'* x)

    optimize!(model) ; solved_time = round(solve_time(model), digits = 2)

    println("solved time $(solved_time)" )

    status = termination_status(model)
    if status != MOI.OPTIMAL
        @info "mono instance is not feasible"
        return 
    end



    println("\n -----------------------------")
    println(" solving $(fname) by $method  ... ")
    println(" -----------------------------")

    folder = "../../results/KP_Forget"
    if !isdir(folder)
        mkdir(folder)
    end

    mixed = false

    heuristic = false

    mixed2 = true

    folder = "../../results/KP_Forget/mixed2"
    if !isdir(folder)
        mkdir(folder)
    end

    result_folder = folder * "/" * string(method)
    if !isdir(result_folder)
        mkdir(result_folder)
    end

    outputName = result_folder * "/" * kp.name * ".dat"
    if isfile(outputName) return end #TODO : ignore existed file  

    vOptKP(kp, Symbol(method), kp.name, outputName, heuristic, mixed,  mixed2)
end


solve(ARGS[1], ARGS[2])




