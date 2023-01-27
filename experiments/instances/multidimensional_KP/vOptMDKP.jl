# ----------------------------------------------------------
# Multi-dimensional knapsack problem 
# Max  sum{j=1,...,n} p(j)x(j)
# st   sum{j=1,...,n} r(i,j)x(j) <= b(i)       i=1,...,m
#                     x(j)=0 or 1
# ----------------------------------------------------------

include("parserMDKP.jl")

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



function vopt_solve(inst::MDKP, method; step=0.5) # fname, outputName
    # ---- setting the model
    model = vModel( CPLEX.Optimizer ) ; JuMP.set_silent(model)

    @variable(model, x[1:inst.n], Bin)
    @constraint(model, [i in 1:inst.m], x'* inst.A[i, :] ≤ inst.b[i])
    @addobjective(model, Max, x'* inst.c)
    #todo : construct C2 
    c2 = generateC2(inst.c)
    @addobjective(model, Max, x'* c2)

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
        println("epsilon_ctr total_time = $total_time ")

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


    # (method != :dicho && method != :epsilon) ? writeResults(inst.n, inst.m, fname, outputName, method, Y_N, X_E; infos) :
    #     writeResults(inst.n, inst.m, fname, outputName, method, Y_N, X_E; total_time)

end




function solve(fname::String)
    instances = readMDKP(fname) 

    for inst in instances
        if inst.n >= 1000 continue end 

        println("\n -----------------------------")
        println(" solving mono $(inst.name) ... ")
        println(" -----------------------------")
 
        model = Model(CPLEX.Optimizer) ; JuMP.set_silent(model)
        @variable(model, x[1:inst.n], Bin )
        @objective(model, Max, x'* inst.c)

        @constraint(model, [i in 1:inst.m], x'* inst.A[i, :] ≤ inst.b[i])

        # optimize
        optimize!(model)  ; solved_time = round(solve_time(model), digits = 2)
        println(" n = $(inst.n) , m = $(inst.m)")
        println("solved time $(solved_time)" )

        if solved_time <= 500.0
            vopt_solve(inst, :epsilon)
        end 
    end
end

solve(ARGS[1])