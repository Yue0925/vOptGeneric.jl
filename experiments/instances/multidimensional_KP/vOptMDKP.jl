# ----------------------------------------------------------
# Multi-dimensional knapsack problem 
# Max  sum{j=1,...,n} p(j)x(j)
# st   sum{j=1,...,n} r(i,j)x(j) <= b(i)       i=1,...,m
#                     x(j)=0 or 1
# ----------------------------------------------------------

include("parserMDKP.jl")

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
    instances = readMDKP(fname) 

    for inst in instances
        if inst.n >= 1000 continue end 

        println("\n -----------------------------")
        println(" solving mono $(inst.name) ... ")
        println(" -----------------------------")
 
        # model = Model(CPLEX.Optimizer) ; JuMP.set_silent(model)
        # @variable(model, x[1:inst.n], Bin )
        # @objective(model, Max, x'* inst.c)

        # @constraint(model, [i in 1:inst.m], x'* inst.A[i, :] ≤ inst.b[i])

        # # optimize
        # optimize!(model)  ; solved_time = round(solve_time(model), digits = 2)
        # println(" n = $(inst.n) , m = $(inst.m)")
        # println("solved time $(solved_time)" )

        # if solved_time <= 300.0
            c2 = generateC2(inst.c)
            folder = "./objective"
            if !isdir(folder)
                mkdir(folder)
            end

            outputName = folder * "/" * inst.name
            fout = open(outputName, "w")
            println(fout, "c2 = $c2 ")
            close(fout)
            
        #     vopt_solve(inst, :epsilon)
        # end 
    end
end

solve(ARGS[1])