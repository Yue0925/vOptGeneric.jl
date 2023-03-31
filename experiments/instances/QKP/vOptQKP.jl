
using JuMP, CPLEX

include("../../../src/vOptGeneric.jl")
using .vOptGeneric 



function writeResults(vars::Int64, constr::Int64, fname::String, outputName::String, method, Y_N, X_E; total_time=nothing, infos=nothing)

    fout = open(outputName, "w")
    println(fout, "vars = $vars ; constr = $constr ")
  
    if method != :dicho && method != :epsilon
        println(fout, infos)
    else
      println(fout, "total_times_used = $total_time")
    end

    println(fout, "size_Y_N = ", length(Y_N))
    println(fout, "Y_N = ", Y_N)
    println(fout)
    println(fout, "size_X_E = ", length(X_E))
  
    close(fout)
  
    # displayGraphics(fname,Y_N, outputName)
end



function vopt_solve(method, outputName; step=0.5) 
    # ---- setting the model
    model = vModel( CPLEX.Optimizer ) ; JuMP.set_silent(model)

    # todo : binary for epsilon 
    @variable(model, x[i=1:n, j=1:i], Bin )
    @variable(model, y[1:n], Bin)

    @addobjective(model, Max, sum([Q1[i, j]*x[i, j] for i=1:n for j=1:i]))
    @addobjective(model, Max, sum([Q2[i, j]*x[i, j] for i=1:n for j=1:i]))

    @constraint(model, y'*w ≤ W)
    @constraint(model, [i=1:n, j=1:i], x[i, j] ≥ y[i] + y[j] -1 )
    @constraint(model, [i=1:n, j=1:i], x[i, j] ≤ y[i])
    @constraint(model, [i=1:n, j=1:i], x[i, j] ≤ y[j])
    @constraint(model, sum(y) == k)


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
        println(" total_time = $total_time ")

    elseif method==:epsilon 
        start = time()
        vSolve( model, method=:epsilon, step=step, verbose=false )
        total_time = round(time() - start, digits = 2)
        println(" total_time = $total_time ")

    elseif method == :bc_rootRelax 
        infos = vSolve( model, method=:bc_rootRelax, verbose=false )
        println(infos)
    elseif method == :bc_rootRelaxCP  
        infos = vSolve( model, method=:bc_rootRelaxCP , verbose=false )
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


    (method != :dicho && method != :epsilon) ? writeResults(round(Int64, n*(n+1)/2+n), round(Int64,n*(n+1)*3/2+1), inst_name, outputName, method, Y_N, X_E; infos) :
        writeResults(round(Int64, n*(n+1)/2+n), round(Int64,n*(n+1)*3/2+1), inst_name, outputName, method, Y_N, X_E; total_time)

end



function solve(fname::String, method::String)
    include(fname)

    folder = "../../results/QKP"
    if !isdir(folder)
        mkdir(folder)
    end
    result_folder = folder * "/" * string(method)
    if !isdir(result_folder)
        mkdir(result_folder)
    end

    println("\n -----------------------------")
    println(" solving mono $(inst_name) ... ")
    println(" -----------------------------")

    model = Model(CPLEX.Optimizer) ; JuMP.set_silent(model)
    @variable(model, x[i=1:n, j=1:i], Bin )
    @variable(model, y[1:n], Bin)

    @objective(model, Max, sum([Q1[i, j]*x[i, j] for i=1:n for j=1:i]))

    @constraint(model, y'*w ≤ W)
    @constraint(model, [i=1:n, j=1:i], x[i, j] ≥ y[i] + y[j] -1 )
    @constraint(model, [i=1:n, j=1:i], x[i, j] ≤ y[i])
    @constraint(model, [i=1:n, j=1:i], x[i, j] ≤ y[j])
    @constraint(model, sum(y) == k)

    # optimize
    optimize!(model) ; solved_time = round(solve_time(model), digits = 2)
    println(" n = $(n*(n+1)/2 + n)")
    println("solved time $(solved_time)" )

    status = termination_status(model)
    if status != MOI.OPTIMAL
        @info "mono instance is not feasible"
        return 
    end

    println("\n -----------------------------")
    println(" solving $(inst_name) by $method  ... ")
    println(" -----------------------------")
    # solve bo-pb 
    outputName = result_folder * "/" * inst_name
    # if isfile(outputName) return end
    vopt_solve(Symbol(method), outputName)

end


solve(ARGS[1], ARGS[2])
