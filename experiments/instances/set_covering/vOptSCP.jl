
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

    folder = "./objective"
    if !isdir(folder)
        mkdir(folder)
    end
    outputName = folder * "/" * inst.name
    fout = open(outputName, "w")
    c1 = [rand(1:inst.n) for _ in 1:inst.n]
    c2 = generateC2(c1)
    println(fout, "c1 = $(c1) ")
    println(fout, "c2 = $c2 ")
    close(fout)

end

solve(ARGS[1])
