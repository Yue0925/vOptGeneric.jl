
include("parserSC.jl")

using JuMP, CPLEX


# include("../../../src/vOptGeneric.jl")
# using .vOptGeneric 


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
    inst = readInstances(fname)

    if inst.n >= 1000 return end 

    println("\n -----------------------------")
    println(" solving mono $(inst.name) ... ")
    println(" -----------------------------")

    # model = Model(CPLEX.Optimizer) ; JuMP.set_silent(model)
    # @variable(model, x[1:inst.n], Bin )
    # @objective(model, Min, x'* inst.c)
    # @constraint(model, [i in 1:inst.m], sum(x[j] for j in inst.cover[i]) >= 1)

    # # optimize
    # optimize!(model) ; solved_time = round(solve_time(model), digits = 2)
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
    # end 

end

solve(ARGS[1])
