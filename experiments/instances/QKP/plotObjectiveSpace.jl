using LinearAlgebra
using PyPlot

# using JuMP, CPLEX

# include("../../../src/vOptGeneric.jl")
# using .vOptGeneric 

function is_SDP(fname::String)
    include(fname)
    println("Q1 is SDP ? ", isposdef(Q1))
    println("Q2 is SDP ? ", isposdef(Q2))
end


function write_ieq(fname::String)
    include(fname)

    f_ieq = "./porta/" * inst_name * ".ieq"
    fout = open(f_ieq, "w")

    println(fout, "DIM = ", n) ; println(fout)

    println(fout, "LOWER_BOUNDS")

    for i in 1:n
        print(fout, "0 ")
    end
    println(fout)

    println(fout, "UPPER_BOUNDS")
    for i in 1:n
        print(fout, "1 ")
    end
    println(fout)


    println(fout)
    println(fout, "INEQUALITIES_SECTION")

    # wx <= W
    print(fout, "$(w[1])x1 ")
    for i in 2:n
        print(fout, "+ $(w[i])x$i ")
    end
    print(fout, "<= ")

    println(fout, W)

    println(fout)

    # k-cluster
    print(fout, "x1 ")
    for i in 2:n
        print(fout, "+ x$i ")
    end
    print(fout, "== ")

    println(fout, k)

    println(fout)

    println(fout, "END")

    close(fout)

end

# nonconvex QP continous relaxation using cplex 
function polyhedron_Y_supp(fname::String)
    include(fname)

    # x = zeros(n)
    dom = [i for i in 0:0.1:1]
    Y_1 = []
    Y_2 = []

    for x1 in dom
        for x2 in dom
            for x3 in dom
                for x4 in dom
                    for x5 in dom
                        for x6 in dom
                            for x7 in dom
                                x = [x1, x2, x3, x4, x5, x6, x7]
                                if sum(x) == k && sum(x.*w) <= W 
                                    push!(Y_1, dot(x, Q1*x))
                                    push!(Y_2, dot(x, Q2*x))
                                end  
                            end
                        end
                    end
                end
            end
        end
    end

    f = "./continousY/" * inst_name
    fout = open(f, "w")
    println(fout, "Y_1 = ", Y_1)
    println(fout, "Y_2 = ", Y_2)

    close(fout)
    plot_polyY(f)
end

function feasible_set_Y(fname::String)
    include(fname)

    feasibleX = Vector{Vector{Int64}}()
    feasibleY = Vector{Vector{Int64}}()

    filename = "./porta/" * inst_name * ".poi"

    f = open(filename)
    lines = readlines(f)
    close(f)

    for line in lines
        if length(line) == 0 || line[1] != '(' continue end 
        l = parse.(Int64, split(split(line, ")")[end], " ")[2:end-1])

        @assert length(l) == n

        push!(feasibleX, l)
        push!(feasibleY, [dot(l, Q1*l), dot(l, Q2*l)])
    end

    # println("feasibleX : ", feasibleX)
    # println("feasibleY : ", feasibleY)

    f_Y = "./setY/" * inst_name
    fout = open(f_Y, "w")
    println(fout, "feasibleY = ", feasibleY)

    close(fout)
end

function plot_Y(fname::String)
    include(fname)

    # scatter points Y 
    f_Y = "./setY/" * inst_name
    include(f_Y)

    Y_1=[];Y_2=[]
    for i in 1:length(feasibleY)
        push!(Y_2, feasibleY[i][2])
        push!(Y_1, feasibleY[i][1])
    end

    figure(inst_name) #  figsize=(6.5,5)
    xlabel(L"z^1(x)")
    ylabel(L"z^2(x)")
    PyPlot.title("BO-kQKP")
    scatter(Y_1, Y_2, color="black", marker="o", label = L"y \in \mathcal{Y}")

    legend(bbox_to_anchor=[1,1], loc=0, borderaxespad=0, fontsize = "x-small")
    savefig("./setY/" * inst_name * ".png")
    PyPlot.close()

end



function plot_polyY(fname::String)
    include(fname)

    # scatter points Y 
    inst_name = split(fname, "/")[end]


    figure(inst_name) #  figsize=(6.5,5)
    xlabel(L"z^1(x)")
    ylabel(L"z^2(x)")
    PyPlot.title("BO-kQKP")
    scatter(Y_1, Y_2, color="blue", marker="+", label = L"\tilde{\mathcal{Y}}")

    legend(bbox_to_anchor=[1,1], loc=0, borderaxespad=0, fontsize = "x-small")
    savefig("./continousY/" * inst_name * ".png")
    PyPlot.close()

end

# function solve(fname::String, method::String)
#     include(fname)

#     folder = "../../results/QKP"
#     if !isdir(folder)
#         mkdir(folder)
#     end
#     result_folder = folder * "/" * string(method)
#     if !isdir(result_folder)
#         mkdir(result_folder)
#     end

#     # if n%3 != 0 return end # todo : 
#     println("\n -----------------------------")
#     println(" solving mono $(inst_name) ... ")
#     println(" -----------------------------")

#     model = Model(CPLEX.Optimizer) ; JuMP.set_silent(model)
#     @variable(model, x[i=1:n, j=1:i], Bin )
#     @variable(model, y[1:n], Bin)

#     @objective(model, Max, sum([Q1[i, j]*x[i, j] for i=1:n for j=1:i]))

#     @constraint(model, y'*w ≤ W)
#     @constraint(model, [i=1:n, j=1:i], x[i, j] ≥ y[i] + y[j] -1 )
#     @constraint(model, [i=1:n, j=1:i], x[i, j] ≤ y[i])
#     @constraint(model, [i=1:n, j=1:i], x[i, j] ≤ y[j])
#     @constraint(model, sum(y) == k)

#     # optimize
#     optimize!(model) ; solved_time = round(solve_time(model), digits = 2)
#     println(" n = $(n*(n+1)/2 + n)")
#     println("solved time $(solved_time)" )

#     status = termination_status(model)
#     if status != MOI.OPTIMAL
#         @info "mono instance is not feasible"
#         return 
#     end

#     println("\n -----------------------------")
#     println(" solving $(inst_name) by $method  ... ")
#     println(" -----------------------------")
#     # solve bo-pb 
#     outputName = result_folder * "/" * inst_name
#     # if isfile(outputName) return end
#     vopt_solve(Symbol(method), outputName)

# end


# write_ieq(ARGS[1])
# feasible_set_Y(ARGS[1])
# plot_Y(ARGS[1])
# is_SDP(ARGS[1])
polyhedron_Y_supp(ARGS[1])