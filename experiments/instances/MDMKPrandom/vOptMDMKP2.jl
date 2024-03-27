
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
    # println(fout, "size_X_E = ", length(X_E))
  
    close(fout)
  
    # displayGraphics(fname,Y_N, outputName)
end



function vopt_solve(method, outputName, limit; step=0.5) 
    # ---- setting the model
    model = vModel( CPLEX.Optimizer ) ; JuMP.set_silent(model)

    @variable(model, x[1:n], Bin )
    @addobjective(model, Max, x'* c)
    @addobjective(model, Max, x'* c2)

    @constraint(model, [i in 1:m], A[i, :]'*x ≤ b[i])
    @constraint(model, [i in 1+m:m+q], A[i, :]'*x ≥ b[i])

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
        println(" total_time = $total_time ")

    elseif method==:epsilon 
        start = time()
        vSolve( model, method=:epsilon, step=step, verbose=false )
        total_time = round(time() - start, digits = 2)
        println(" total_time = $total_time ")

    elseif method == :bc_rootRelax 
        infos = vSolve( model, method=:bc_rootRelax, verbose=false , LBSexhaustive = exhaustive, λ_limit = limit)
        println(infos)
    elseif method == :bc_rootRelaxCP  
        infos = vSolve( model, method=:bc_rootRelaxCP , verbose=false , LBSexhaustive = exhaustive,  λ_limit = limit )
        println(infos)    
    elseif method == :bb_EPB 
        infos = vSolve( model, method=:bb_EPB, verbose=false )
        println(infos)
    elseif method == :bc_rootRelaxEPB
        infos = vSolve( model, method=:bc_rootRelaxEPB, verbose=false , LBSexhaustive = exhaustive,  λ_limit = limit)
        println(infos)
    elseif method == :bc_rootRelaxCPEPB
        infos = vSolve( model, method=:bc_rootRelaxCPEPB, verbose=false , LBSexhaustive = exhaustive, λ_limit = limit)
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


    (method != :dicho && method != :epsilon) ? writeResults(n, m, inst_name, outputName, method, Y_N, X_E; infos) :
        writeResults(n, m, inst_name, outputName, method, Y_N, X_E; total_time)

end



function solve(fname::String, method::String)
    include(fname)

    folder = "../../results/MDMKPrandom"
    if !isdir(folder)
        mkdir(folder)
    end

    # todo : 
    if !(n == 10 || n==20 || n== 30 || n== 40 ) return end

    println("\n -----------------------------")
    println(" solving mono $(name) ... ")
    println(" -----------------------------")
        
    model = Model(CPLEX.Optimizer) ; JuMP.set_silent(model)
    @variable(model, x[1:n], Bin )
    @objective(model, Max, x'* c)
    @constraint(model, [i in 1:m], A[i, :]'*x ≤ b[i])
    @constraint(model, [i in 1+m:m+q], A[i, :]'*x ≥ b[i])


    # optimize
    optimize!(model) ; solved_time = round(solve_time(model), digits = 2)
    println(" n = $(n) , m = $(m), q = $(q)")
    println("solved time $(solved_time)" )

    status = termination_status(model)
    if status != MOI.OPTIMAL
        @info "mono instance is not feasible"
        return 
    end

    for limit in [2, 3, 6, 8] # 

        println("\n -----------------------------")
        println(" solving $(name) by $method λ_limit = $limit ... ")
        println(" -----------------------------")

        result_folder = folder * "/" * "lambda_limit"
        if !isdir(result_folder)
            mkdir(result_folder)
        end

        result_folder = folder * "/" * "lambda_limit" * "/" * string(limit)
        if !isdir(result_folder)
            mkdir(result_folder)
        end

        result_folder = folder * "/" * "lambda_limit" * "/" * string(limit) * "/" * string(method)
        if !isdir(result_folder)
            mkdir(result_folder)
        end

        # solve bo-pb 
        outputName = result_folder * "/" * inst_name

        # todo : if the output file already exists 
        if isfile(outputName)
            return
        end
        
        # if isfile(outputName) return end
        vopt_solve(Symbol(method), outputName, limit)

    end


end


solve(ARGS[1], ARGS[2])
