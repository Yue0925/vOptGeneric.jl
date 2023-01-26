# ----------------------------------------------------------
# Multi-dimensional knapsack problem 
# Max  sum{j=1,...,n} p(j)x(j)
# st   sum{j=1,...,n} r(i,j)x(j) <= b(i)       i=1,...,m
#                     x(j)=0 or 1
# ----------------------------------------------------------

include("parserMDKP.jl")

using JuMP, CPLEX



function solve(fname::String)
    instances = readMDKP(fname) 

    for inst in instances
        println("\n -----------------------------")
        println(" building $(inst.name) ... ")
        println(" -----------------------------")
 
        model = Model(CPLEX.Optimizer) ; JuMP.set_silent(model)
        @variable(model, x[1:inst.n], Bin )
        @objective(model, Max, x'* inst.c)

        @constraint(model, [i in 1:inst.m], x'* inst.A[i, :] ≤ inst.b[i])

        # optimize
        optimize!(model) 
        println(" n = $(inst.n) , m = $(inst.m)")
        println("solved time $(round(solve_time(model), digits = 2))" )
    end
end

solve(ARGS[1])