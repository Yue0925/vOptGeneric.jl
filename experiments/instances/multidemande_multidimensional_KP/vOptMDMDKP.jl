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



function solve(fname::String)
    instances = readInstances(fname)

    for inst in instances
        if inst.n >= 1000 continue end 

        println("\n -----------------------------")
        println(" solving mono $(inst.name) ... ")
        println(" -----------------------------")

        # model = Model(CPLEX.Optimizer) ; JuMP.set_silent(model)
        # @variable(model, x[1:inst.n], Bin )
        # @objective(model, Max, x'* inst.c)

        # @constraint(model, [i in 1:inst.m], x'* inst.A_inf[i, :] ≤ inst.b_inf[i])

        # @constraint(model, [i in 1:inst.m], x'* inst.A_sup[i, :] ≥ inst.b_sup[i])

        # # optimize
        # optimize!(model) ; solved_time = round(solve_time(model), digits = 2)
        # println(" n = $(inst.n) , m = $(inst.m * 2)")
        # println("solved time $(solved_time)" )

        c2 = generateC2(inst.c)
        folder = "./objective"
        if !isdir(folder)
            mkdir(folder)
        end

        outputName = folder * "/" * inst.name
        fout = open(outputName, "w")
        println(fout, "c2 = $c2 ")
        close(fout)

    end
end

solve(ARGS[1])