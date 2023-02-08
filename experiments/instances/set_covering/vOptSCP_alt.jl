
include("parserSC.jl")

using JuMP, CPLEX


include("../../../src/vOptGeneric.jl")
using .vOptGeneric 


"""
Coefficients generator proposed by Pedersen et al, 2008.
"""
function generateC2(c1::Vector{Int64})::Vector{Int64}
    mini = minimum(c1) ; maxi = maximum(c1)
    c2 = zeros(Int64, length(c1))

    for i in 1:length(c1)
        if c1[i] < (maxi-mini)/2
            c2[i] = rand((maxi-c1[i]+mini):(maxi))
        elseif c1[i] >= (maxi-mini)/2
            c2[i] = rand((mini):(maxi-c1[i]+mini))
        else
            error("Coefficient error c1[$i] = $(c1[i]), maximum = $maxi, minimum = $mini !")
        end
    end
    return c2
end


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



function vopt_solve(inst::SCP, method, outputName; step=0.5) # fname, outputName
    # ---- setting the model
    model = vModel( CPLEX.Optimizer ) ; JuMP.set_silent(model)

    @variable(model, x[1:inst.n], Bin)
    @constraint(model, [i in 1:inst.m], sum(x[j] for j in inst.cover[i]) >= 1)
    
    include("./objective/" * inst.name)
    @addobjective(model, Min, x'* c1)
    @addobjective(model, Min, x'* c2)

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


    (method != :dicho && method != :epsilon) ? writeResults(inst.n, inst.m, inst.name, outputName, method, Y_N, X_E; infos) :
        writeResults(inst.n, inst.m, inst.name, outputName, method, Y_N, X_E; total_time)

end




function solve(fname::String, method::String)
    inst = readInstances(fname)

    folder = "../../results/BOSCP"
    if !isdir(folder)
        mkdir(folder)
    end
    result_folder = folder * "/" * string(method)
    if !isdir(result_folder)
        mkdir(result_folder)
    end

    fileName = "./objective/" * inst.name
    if !isfile(fileName) return end

    println("\n -----------------------------")
    println(" solving mono $(inst.name) ... ")
    println(" -----------------------------")

    model = Model(CPLEX.Optimizer) ; JuMP.set_silent(model)
    @variable(model, x[1:inst.n], Bin )
    include("./objective/" * inst.name)
    @objective(model, Min, x'* c1)
    @constraint(model, [i in 1:inst.m], sum(x[j] for j in inst.cover[i]) >= 1)

    relax_integrality(model)
    # optimize
    optimize!(model) ; solved_time = round(solve_time(model), digits = 2)
    println(" n = $(inst.n) , m = $(inst.m)")
    println("solved time $(solved_time)" )

    println("\n -----------------------------")
    println(" solving $(inst.name) by $method  ... ")
    println(" -----------------------------")
    # solve bo-pb 
    outputName = result_folder * "/" * inst.name
    if isfile(outputName) return end
    vopt_solve(inst, Symbol(method), outputName)

end


solve(ARGS[1], ARGS[2])
