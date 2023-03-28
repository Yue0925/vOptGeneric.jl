using JuMP, CPLEX 


function colpmeleQ2()
    folder = "./instances/"
    for file in readdir(folder)
        include(folder * file)
        println(file * "...")
        fout = open(folder * file, "a")
        println(fout, "inst_name = \"$file\" ")

        # println("w=$w", w) ; exit(0)

        # mono is feasible 
        println("\n -----------------------------")
        println(" solving mono $(file) ... ")
        println(" -----------------------------")
    
        model = Model(CPLEX.Optimizer) ; JuMP.set_silent(model)
        @variable(model, x[i=1:n, j=1:i], Bin )
        @variable(model, y[1:n], Bin)

        @objective(model, Max, sum([Q1[i, j]*x[i, j] for i=1:n for j=1:i]))

        @constraint(model, y'*w ≤ W)
        @constraint(model, [i=1:n, j=1:i], x[i, j] ≥ y[i] + y[j] -1 )
        @constraint(model, [i=1:n, j=1:i], x[i, j] ≤ y[i])
        @constraint(model, [i=1:n, j=1:i], x[i, j] ≤ y[j])
    
    
        # relax_integrality(model)
        # optimize
        optimize!(model) ; solved_time = round(solve_time(model), digits = 2)
        println(" n = $(n)")
        println("solved time $(solved_time)" )
        status = JuMP.termination_status(model)
    
        if status != MOI.OPTIMAL
            @info "Not OPT ..."
            println(fout, "feasible = false ")
            close(fout)
            continue
        end


        println(fout, "feasible = true ")

        Q2 = zeros(Int64, n, n) ; maxi = maximum(Q1)
        len = n*(n-1)/2 ; mini = minimum(Q1)
        for i in 1:n 
            for j in 1:i
                # if Q1[i, j] == 0
                #     nothing
                if Q1[i, j] < (maxi-mini)/2
                    Q2[i, j] = rand((maxi-Q1[i, j]+mini):(maxi)) 
                elseif Q1[i, j] >= (maxi-mini)/2
                    Q2[i, j] = rand((mini):(maxi-Q1[i, j]+mini))
                else
                    error("Coefficient error Q1[$i, $j] = $(Q1[i, j]), maximum = $maxi, minimum = $mini !")
                end            
            end
        end

        println(fout, "Q2=$Q2")
        close(fout)
    end
end

colpmeleQ2()