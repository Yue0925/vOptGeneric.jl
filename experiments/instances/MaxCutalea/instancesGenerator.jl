"""
Max cut instances generation for complete graphs Kₙ with random edge weights (between 500 and 1000). 
"""


using JuMP, CPLEX 


function completeG(n::Int64)
    W = zeros(Int64, n, n) ; mini = 1000
    for i in 1:n-1
        for j in i+1:n
            W[i, j] = rand(500:1000)
            mini = (mini>W[i, j]) ? W[i, j] : mini
        end        
    end
    W_bis = zeros(Int64, n, n) ; maxi = maximum(W)
    len = n*(n-1)/2
    for i in 1:n-1 
        for j in i+i:n
            if W[i, j] < (maxi-mini)/2
                W_bis[i, j] = rand((maxi-W[i, j]+mini):(maxi)) 
            elseif W[i, j] >= (maxi-mini)/2
                W_bis[i, j] = rand((mini):(maxi-W[i, j]+mini))
            else
                error("Coefficient error W[$i, $j] = $(W[i, j]), maximum = $maxi, minimum = $mini !")
            end            
        end
    end

    # --------------------
    name = "MaxCut_$n"
    folder = "./instances/"
    if !isdir(folder)
        mkdir(folder)
    end

    println("\n -----------------------------")
    println(" solving mono $(name) ... ")
    println(" -----------------------------")

    model = Model(CPLEX.Optimizer) ; JuMP.set_silent(model)
    @variable(model, x[i=1:n-1, j=i+1:n], Bin )
    @variable(model, y[1:n], Bin)
    @objective(model, Max, sum([W[i, j]*x[i, j] for i=1:n-1 for j=i+1:n]))
    @constraint(model, [i=1:n-1, j=i+1:n], x[i, j] ≤ y[i] + y[j])
    @constraint(model, [i=1:n-1, j=i+1:n], x[i, j] ≤ 2-(y[i] + y[j]))


    # relax_integrality(model)
    # optimize
    optimize!(model) ; solved_time = round(solve_time(model), digits = 2)
    println(" n = $(n)")
    println("solved time $(solved_time)" )
    status = JuMP.termination_status(model)

    if status != MOI.OPTIMAL
        @info "Not OPT ..."
        return #continue
    end

    fout = open(folder * name, "w")
    println(fout, "inst_name = \"$name\" ")
    println(fout, "n = $(n)")

    println(fout, "W = $W")
    println(fout, "W_bis = $W_bis")
    close(fout)

end



for n in 5:9
    completeG(n)
end