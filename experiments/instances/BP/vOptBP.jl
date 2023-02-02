
include("parserBP.jl")

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




function vopt_solve(inst::BP, method, outputName; step=0.5) # fname, outputName
    # ---- setting the model
    model = vModel( CPLEX.Optimizer ) ; JuMP.set_silent(model)

    @variable(model, x[1:inst.n], Bin)
    @constraint(model, [i in 1:inst.m], x'* inst.A[i, :] ≤ inst.b[i])

    include("./objective/" * inst.name)
    inst.objSense == MIN_SENSE ? @addobjective(model, Min, x'* c1) : @addobjective(model, Max, x'* c1)
    inst.objSense == MIN_SENSE ? @addobjective(model, Min, x'* c2) : @addobjective(model, Max, x'* c2)

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
    inst = readModel(fname)

    folder = "../../results/BP"
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
    println(" solving $(inst.name) by $method  ... ")
    println(" -----------------------------")

    # solve bo-pb 
    outputName = result_folder * "/" * inst.name
    # todo avoid double computation 
    # if isfile(outputName) return end
    vopt_solve(inst, Symbol(method), outputName)
end


solve(ARGS[1], ARGS[2])
