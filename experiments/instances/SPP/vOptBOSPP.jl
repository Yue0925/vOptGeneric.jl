
const verbose = false
const graphic = true

using JuMP, CPLEX
include("../../../src/vOptGeneric.jl")
using .vOptGeneric 

include("parserBOSPP.jl")


function solve(fname::String, method ; step=0.5)
    if !isfile(fname)
        @error "This file doesn't exist ! $fname"
    end

    result_dir = "../../results/SPP/" * split(fname, "/")[2]
    if !isdir(result_dir)
        mkdir(result_dir)
    end

    folder = result_dir * "/" * string(method)
    if !isdir(folder)
        mkdir(folder)
    end

    outputName = folder * "/" * split(fname, "/")[end]
    # if isfile(outputName) return end #TODO : ignore existed file  
    inst = readingBOSPP(fname)
    if inst.n >= 3000 return end
    

    # ---- setting the model
    println("Building...")
    bospp = vModel( CPLEX.Optimizer ) ; JuMP.set_silent(bospp)
    @variable( bospp, x[1:inst.n], Bin )
    @addobjective( bospp, Max, sum(inst.C1[j] * x[j] for j=1:inst.n) )
    @addobjective( bospp, Max, sum(inst.C2[j] * x[j] for j=1:inst.n) )

    @constraint( bospp, cte[i=1:inst.m], sum(x[j] for j in inst.Cover[i]) <= 1)

    println("Solving...")



    if method == :bb
        infos = vSolve( bospp, method=:bb, verbose=false )
        println(infos)
    elseif method == :bc 
        infos = vSolve( bospp, method=:bc, verbose=false )
        println(infos)
    elseif method == :dicho 
        start = time()
        vSolve( bospp, method=:dicho, verbose=false)
        total_time = round(time() - start, digits = 2)
    elseif method==:epsilon 
        start = time()
        vSolve( bospp, method=:epsilon, step=0.5, verbose=false)
        total_time = round(time() - start, digits = 2)

    elseif method == :bc_rootRelax 
        infos = vSolve( bospp, method=:bc_rootRelax, verbose=false )
        println(infos)
    elseif method == :bc_rootRelaxCP  
        infos = vSolve( bospp, method=:bc_rootRelaxCP , verbose=false )
        println(infos)    
    elseif method == :bb_EPB 
        infos = vSolve( bospp, method=:bb_EPB, verbose=false )
        println(infos)
    elseif method == :bc_rootRelaxEPB
        infos = vSolve( bospp, method=:bc_rootRelaxEPB, verbose=false )
        println(infos)
    elseif method == :bc_rootRelaxCPEPB
        infos = vSolve( bospp, method=:bc_rootRelaxCPEPB, verbose=false )
        println(infos)
    elseif method == :bc_EPB
        infos = vSolve( bospp, method=:bc_EPB, verbose=false )
        println(infos)
    else
        @error "Unknown method parameter $(method) !"
    end


    # --------------------------------------------


    # ---- Querying the results
    println("Querying...")
    Y_N = getY_N( bospp )
    println("length Y_N = ", length(Y_N))

    X_E = getX_E( bospp )
    println("length X_E = ", length(X_E))

    # ---- Writing the results
    (method != :dicho && method != :epsilon) ? writeResults(inst.n, inst.m, fname, outputName, method, Y_N, X_E; infos) :
        writeResults(inst.n, inst.m, fname, outputName, method, Y_N, X_E; total_time)

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

