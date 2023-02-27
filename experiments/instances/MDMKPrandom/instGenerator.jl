"""
MDMKP       Max ∑ⱼ cⱼxⱼ
 s.t.          ∑ⱼ aᵢⱼ xⱼ ≤ bᵢ   ∀ i=1...m 
               ∑ⱼ aᵢⱼ xⱼ ≥ bᵢ   ∀ i=m+1...m+q 
               xⱼ ∈ {0,1}

!!! coefficients are not necessarily integers 


n = 10, 20 ... 50
m = 2, 4, ... n/5
q = [1, m/2, m]

"""

using JuMP, CPLEX 

vars = [10, 20, 30, 40, 50]
α_ = [0.25, 0.5, 0.75]


function generateC2(c1::Vector{Float64})::Vector{Float64}
    mini = minimum(c1) ; maxi = maximum(c1)
    c2 = zeros(length(c1))

    for i in 1:length(c1)
        if c1[i] < (maxi-mini)/2
            c2[i] = rand((maxi-c1[i]+mini):(maxi))      # todo : not work for negative coeff for now 
        elseif c1[i] >= (maxi-mini)/2
            c2[i] = rand((mini):(maxi-c1[i]+mini))
        else
            error("Coefficient error c1[$i] = $(c1[i]), maximum = $maxi, minimum = $mini !")
        end
    end
    return c2
end



"""
MDMKP instances with positive cost coefficients. 
"""
function instancesPCC(n::Int64, α::Float64)
    for m in 2:2:round(Int64, n/5)
        for q in round.(Int64, Set([1, m/2, m]))

            A = zeros(Int64, m+q, n)
            for i in 1:m+q for j in 1:n A[i, j] = rand(1:1000) end end 
            
            b = [α* sum(A[i, :]) for i in 1:m+q]

            c = [sum(A[1:m, j])/m - sum(A[m+1:m+q, j])/m for j in 1:n]

            Δ = minimum(c)<0 ? 1+abs(minimum(c)) : 0 

            for j in 1:n
                c[j] = c[j] + Δ + 500*rand(0.0:0.01:1.0)
            end

            name = "MDMKP_n$(n)_m$(m)_q$(q)_α$(α)"
            folder = "./instances/PCC/"
            if !isdir(folder)
                mkdir(folder)
            end

            # --------------------

            println("\n -----------------------------")
            println(" solving mono $(name) ... ")
            println(" -----------------------------")
        
            model = Model(CPLEX.Optimizer) ; JuMP.set_silent(model)
            @variable(model, x[1:n], Bin )
            @objective(model, Max, x'* c)
            @constraint(model, [i in 1:m], A[i, :]'*x ≤ b[i])
            @constraint(model, [i in 1+m:m+q], A[i, :]'*x ≥ b[i])

        
            # relax_integrality(model)
            # optimize
            optimize!(model) ; solved_time = round(solve_time(model), digits = 2)
            println(" n = $(n) , m = $(m), q = $(q)")
            println("solved time $(solved_time)" )
            status = JuMP.termination_status(model)

            if status != MOI.OPTIMAL
                @info "Not OPT ..."
                continue
            end
        
            # --------------------
            c2 = generateC2(c)
            fout = open(folder * name, "w")
            println(fout, "inst_name = \"$name\" ")
            println(fout, "n = $n ") ; println(fout, "m = $m ")
            println(fout, "q = $q ") 
            println(fout, "A = $A") ; println(fout, "b = $b")
            println(fout, "c = $c ") ; println(fout, "c2 = $c2 ")
            close(fout)
        end
    end
end


function instancesMCC(n::Int64, α::Float64)
    for m in 2:2:round(Int64, n/5)
        for q in round.(Int64, Set([1, m/2, m]))

            A = zeros(Int64, m+q, n)
            for i in 1:m+q for j in 1:n A[i, j] = rand(1:1000) end end 
            
            b = [α* sum(A[i, :]) for i in 1:m+q]

            c = [sum(A[1:m, j])/m - sum(A[m+1:m+q, j])/m + 500*rand(0.0:0.01:1.0) for j in 1:n]
            
            name = "MDMKP_n$(n)_m$(m)_q$(q)_α$(α)"
            folder = "./instances/MCC/"
            if !isdir(folder)
                mkdir(folder)
            end

            # --------------------

            println("\n -----------------------------")
            println(" solving mono $(name) ... ")
            println(" -----------------------------")
        
            model = Model(CPLEX.Optimizer) ; JuMP.set_silent(model)
            @variable(model, x[1:n], Bin )
            @objective(model, Max, x'* c)
            @constraint(model, [i in 1:m], A[i, :]'*x ≤ b[i])
            @constraint(model, [i in 1+m:m+q], A[i, :]'*x ≥ b[i])

        
            # relax_integrality(model)
            # optimize
            optimize!(model) ; solved_time = round(solve_time(model), digits = 2)
            println(" n = $(n) , m = $(m), q = $(q)")
            println("solved time $(solved_time)" )
        
            status = JuMP.termination_status(model)

            if status != MOI.OPTIMAL
                @info "Not OPT ..."
                continue
            end

            # --------------------
            c2 = generateC2(c)
            fout = open(folder * name, "w")
            println(fout, "inst_name = \"$name\" ")
            println(fout, "n = $n ") ; println(fout, "m = $m ")
            println(fout, "q = $q ") 
            println(fout, "A = $A") ; println(fout, "b = $b")
            println(fout, "c = $c "); println(fout, "c2 = $c2 ")
            close(fout)
        end
    end
end



for n in vars
    for α in α_

        instancesPCC(n, α)
        # instancesMCC(n, α)
    end
    
end