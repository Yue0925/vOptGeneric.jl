# =============================================================================
"""
    0/1 uncapacitated facility location problem with two objectives (biUFLP)
    uncapacitatedFacilityLocation.jl
    September 2022
"""

const verbose = false
const graphic = true
# =============================================================================
using JuMP, CPLEX, Printf
include("../../../src/vOptGeneric.jl")
using .vOptGeneric 

include("../../../src/BO01BB/displayGraphic.jl")



# =============================================================================
# structure of a 2UFLP instance

mutable struct Instance
    fname :: String        # name of the file
    nI    :: Int64         # number of users
    nJ    :: Int64         # number of services
    c1    :: Array{Int,2}  # assignment costs users-services for objective 1
    c2    :: Array{Int,2}  # assignment costs users-services for objective 2
    r1    :: Array{Int,1}  # running costs for objective 1
    r2    :: Array{Int,1}  # running costs for objective 2
end


# =============================================================================
"""
    load2UFLP(fdirectory::String, fname::String)
    Load an instance of bi-objective UFLP
      ↓  fdirectory : path to the data file
      ↓  fname      : file name of the data
      ↑             : data of an instance
"""

function load2UFLP(fname::String)

    f=open(fname)
 
    # read the number of users (nI)
    nI = parse(Int, readline(f) )
    # read the number of services (nJ) 
    nJ = parse(Int, readline(f) )
    # read the following line (separator -> line without information)
    useless = readline(f)

    # assignment costs users-services (matrix nI x nJ) ------------------------
    c1 = Array{Int64,2}(undef,nI,nJ)
    c2 = Array{Int64,2}(undef,nI,nJ)

    # objective 1
    for i=1:nI
        c1[i,:] = parse.(Int64, split(readline(f)) )
    end
    # read the following line (separator -> line without information)
    useless = readline(f)

    # objective 2
    for i=1:nI
        c2[i,:] = parse.(Int64, split(readline(f)) )
    end
    # read the following line (separator -> line without information)
    useless = readline(f)

    # running costs of services (vector nJ) -----------------------------------
    r1 = Array{Int64,1}(undef,nJ)
    r2 = Array{Int64,1}(undef,nJ)

    # objective 1
    r1[:] = parse.(Int64, split(readline(f)) )
    # read the following line (separator -> line without information)
    useless = readline(f)

    # objective 2
    r2[:] = parse.(Int64, split(readline(f)) )
   
    close(f)

    return Instance(fname,nI,nJ,c1,c2,r1,r2)

end


# ==============================================================================
"""
    createModel2UFLP(solver::DataType, data::Instance)
    Create the vOptGeneric model of 2UFLP 
"""
function createModel2UFLP(solver::DataType, data::Instance)
    model = vModel( solver )
    
    @variable(model, x[1:data.nI,1:data.nJ], Bin)
    @variable(model, s[1:data.nJ], Bin)

    @addobjective( model, Min, sum(data.c1[i,j]*x[i,j] for i in 1:data.nI, j in 1:data.nJ) + sum(data.r1[j]*s[j] for j in 1:data.nJ) )
    @addobjective( model, Min, sum(data.c2[i,j]*x[i,j] for i in 1:data.nI, j in 1:data.nJ) + sum(data.r2[j]*s[j] for j in 1:data.nJ) )

    @constraint( model, [i=1:data.nI], sum(x[i,j] for j in 1:data.nJ) == 1 )
    @constraint( model, [i=1:data.nI,j=1:data.nJ], x[i,j] <= s[j] )    

    return model
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
  
    displayGraphics(fname,Y_N, outputName)
end



# ==============================================================================
function solve(fname::String, method ;step=0.5)

    data  = load2UFLP(fname)

    # -------------------------------------------------------------------------
    # Display information about the instance to solve

    println("\nfilename      : $(data.fname)") 
    println("nI (users)    : $(data.nI)") 
    println("nJ (services) : $(data.nJ)\n") 

    folder = "../../results/UFLP"
    if !isdir(folder)
        mkdir(folder)
    end
    result_folder = folder * "/" * string(method)
    if !isdir(result_folder)
        mkdir(result_folder)
    end
    inst_name = split(fname, "/")[end]

    outputName = result_folder * "/" * inst_name ; fname = string(split(inst_name, ".")[1])

    # -------------------------------------------------------------------------
    solver = CPLEX.Optimizer    
    mod2UFLP = createModel2UFLP(solver, data)
    set_silent(mod2UFLP)

    n = data.nI * data.nJ + data.nJ
    m = data.nI + data.nI * data.nJ
    println("n = $n m = $m ")


    if method == :bb
        infos = vSolve( mod2UFLP, method=:bb, verbose=false )
        println(infos)
    elseif method == :bc 
        infos = vSolve( mod2UFLP, method=:bc, verbose=false )
        println(infos)
    elseif method == :dicho 
        start = time()
        vSolve( mod2UFLP, method=:dicho, verbose=false )
        total_time = round(time() - start, digits = 2)
    elseif method==:epsilon 
        start = time()
        vSolve( mod2UFLP, method=:epsilon, step=step, verbose=false )
        total_time = round(time() - start, digits = 2)

    elseif method == :bc_rootRelax 
        infos = vSolve( mod2UFLP, method=:bc_rootRelax, verbose=false )
        println(infos)
    elseif method == :bc_rootRelaxCP  
        infos = vSolve( mod2UFLP, method=:bc_rootRelaxCP , verbose=false )
        println(infos)    
    elseif method == :bb_EPB 
        infos = vSolve( mod2UFLP, method=:bb_EPB, verbose=false )
        println(infos)
    elseif method == :bc_rootRelaxEPB
        infos = vSolve( mod2UFLP, method=:bc_rootRelaxEPB, verbose=false )
        println(infos)
    elseif method == :bc_rootRelaxCPEPB
        infos = vSolve( mod2UFLP, method=:bc_rootRelaxCPEPB, verbose=false )
        println(infos)
    elseif method == :bc_EPB
        infos = vSolve( mod2UFLP, method=:bc_EPB, verbose=false )
        println(infos)
    else
        @error "Unknown method parameter $(method) !"
    end

    # -------------------------------------------------------------------------
    # ---- Querying the results
    Y_N = getY_N( mod2UFLP )
    println("length Y_N = ", length(Y_N))

    X_E = getX_E( mod2UFLP )
    println("length X_E = ", length(X_E))


    (method != :dicho && method != :epsilon) ? writeResults(n, m, fname, outputName, method, Y_N, X_E; infos) :
        writeResults(n, m, fname, outputName, method, Y_N, X_E; total_time)
    return nothing
end

# ==============================================================================




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
