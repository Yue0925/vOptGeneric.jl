# ----------------------------------------------------------------------------
# Multi-demande multi-dimensional Knapsack Problem 
#
# (MDMKP) max  sum{j=1,...,n} c(j) x(j) 
#         s.t. sum{j=1,...,n} a(i,j) x(j) <= b(i) i=1,...,m
#              sum{j=1,...,n} a(i,j) x(j) >= b(i) i=m+1,...,m+q
#              x(j) = 0 or 1                      j=1,...,n
#  ---------------------------------------------------------------------------

include("parserMDMDKP.jl")

using JuMP, CPLEX



function solve(fname::String)
    instances = readInstances(fname)

    for inst in instances
        println("\n -----------------------------")
        println(" building $(inst.name) ... ")
        println(" -----------------------------")

        model = Model(CPLEX.Optimizer) ; JuMP.set_silent(model)
        @variable(model, x[1:inst.n], Bin )
        @objective(model, Max, x'* inst.c)

        @constraint(model, [i in 1:inst.m], x'* inst.A_inf[i, :] ≤ inst.b_inf[i])

        @constraint(model, [i in 1:inst.m], x'* inst.A_sup[i, :] ≥ inst.b_sup[i])

        # optimize
        optimize!(model) 
        println(" n = $(inst.n) , m = $(inst.m * 2)")
        println("solved time $(round(solve_time(model), digits = 2))" )
    end
end

solve(ARGS[1])