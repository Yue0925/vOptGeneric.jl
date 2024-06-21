using JuMP, CPLEX

include("../../../src/vOptGeneric.jl")
include("parserForget20.jl")
using .vOptGeneric



function vOptUFLP(forget::Forget,
    method, fname, outputName, heuristic, mixed, mixed2 ; step = 0.01
)

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
        infos = vSolve( model, method=:bb, verbose=false , heuristic = heuristic)
        println(infos)
    elseif method == :bc 
        infos = vSolve( model, method=:bc, verbose=false , heuristic = heuristic)
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
        infos = vSolve( model, method=:bc_rootRelax, verbose=false , heuristic = heuristic, mixed = mixed, mixed2 = mixed2 )
        println(infos)
    elseif method == :bc_rootRelaxCP 
        infos = vSolve( model, method=:bc_rootRelaxCP, verbose=false, heuristic = heuristic , mixed = mixed, mixed2 = mixed2 )
        println(infos)
    elseif method == :bb_EPB 
        infos = vSolve( model, method=:bb_EPB, verbose=false , heuristic = heuristic)
        println(infos)
    elseif method == :bc_rootRelaxEPB 
        infos = vSolve( model, method=:bc_rootRelaxEPB, verbose=false , heuristic = heuristic, mixed = mixed, mixed2 = mixed2 )
        println(infos)
    elseif method == :bc_rootRelaxCPEPB 
        infos = vSolve( model, method=:bc_rootRelaxCPEPB, verbose=false, heuristic = heuristic, mixed = mixed, mixed2 = mixed2 )
        println(infos)
    elseif method == :bc_EPB 
        infos = vSolve( model, method=:bc_EPB, verbose=false , heuristic = heuristic)
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



    println("\n -----------------------------")
    println(" solving $(fname) by $method  ... ")
    println(" -----------------------------")

    folder = "../../results/UFLP_Forget"
    if !isdir(folder)
        mkdir(folder)
    end



    mixed2 = false ; mixed = false

    heuristic = false

    mixed2 = true

    folder = "../../results/UFLP_Forget/mixed2"
    if !isdir(folder)
        mkdir(folder)
    end

    result_folder = folder * "/" * string(method)
    if !isdir(result_folder)
        mkdir(result_folder)
    end

    outputName = result_folder * "/" * forget.name * ".dat"
    # if isfile(outputName) return end #TODO : ignore existed file  

    vOptUFLP(forget, Symbol(method), forget.name, outputName, heuristic, mixed,  mixed2)


    heuristic = true


    folder = "../../results/UFLP_Forget/mixed2heuristic"
    if !isdir(folder)
        mkdir(folder)
    end

    result_folder = folder * "/" * string(method)
    if !isdir(result_folder)
        mkdir(result_folder)
    end


    outputName = result_folder * "/" * forget.name * ".dat"
    # if isfile(outputName) return end #TODO : ignore existed file  

    vOptUFLP(forget, Symbol(method), forget.name, outputName, heuristic, mixed,  mixed2)



end


solve(ARGS[1], ARGS[2])

# solve("./instance/Forget20-UFLP_5_3_1-1000_1-100_spheredown_1_1.raw", :bb)

