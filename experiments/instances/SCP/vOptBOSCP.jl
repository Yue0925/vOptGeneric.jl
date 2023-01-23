const verbose = false
const graphic = true

using JuMP, CPLEX
include("../../../src/vOptGeneric.jl")
using .vOptGeneric 

include("parserBOSCP.jl")



function solveBOSCP(  inst::BOSCP, method, fname, outputName; step=0.5)
    println("n = $(inst.n) m = $(inst.m) ")
    # ---- setting the model
    model = vModel( CPLEX.Optimizer )
    JuMP.set_silent(model)

    @variable(model, x[1:inst.n], Bin)
    @constraint(model, [i=1:inst.m], sum(x[j] for j in inst.Cover[i]) >= 1)
    @addobjective(model, Min, sum(inst.C1[j] * x[j] for j=1:inst.n))
    @addobjective(model, Min, sum(inst.C2[j] * x[j] for j=1:inst.n))

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
    println("length X_E = ", length(X_E))


    (method != :dicho && method != :epsilon) ? writeResults(inst.n, inst.m, fname, outputName, method, Y_N, X_E; infos) :
        writeResults(inst.n, inst.m, fname, outputName, method, Y_N, X_E; total_time)

end


function solve(fname::String, method)

    # load a numerical instance of 2SPA ----------------------------------------
    inst = readingBOSCP(fname)
    if inst.n >= 1000 return end

    folder = "../../results/SCP/BOSCP"
    if !isdir(folder)
        mkdir(folder)
    end
    result_folder = folder * "/" * string(method)
    if !isdir(result_folder)
        mkdir(result_folder)
    end
    inst_name = split(fname, "/")[end]

    outputName = result_folder * "/" * inst_name
    solveBOSCP(inst, method, string(split(inst_name, ".")[1]), outputName)
end

solve(ARGS[1], :dicho)
solve(ARGS[1], :epsilon)
solve(ARGS[1], :bb)
solve(ARGS[1], :bc)

solve(ARGS[1], :bc_rootRelax)
solve(ARGS[1], :bc_EPB)
solve(ARGS[1], :bc_rootRelaxCPEPB)
solve(ARGS[1], :bb_EPB)
solve(ARGS[1], :bc_rootRelaxCP)
solve(ARGS[1], :bc_rootRelaxEPB)

