# MIT License
# Copyright (c) 2017: Xavier Gandibleux, Anthony Przybylski, Gauthier Soleilhac, and contributors.

using JuMP, CPLEX 

function solve_lexico(m::JuMP.Model, verbose; args...)
    #Retrieve objectives and their senses from vOptData
    vd = getvOptData(m)
    empty!(vd.Y_N) ; empty!(vd.X_E) #; empty!(vd.logObjs)
    objs = vd.objs
    objSenses = vd.objSenses
    nbObj = length(objs)

    #Check that problem is feasible and that no objective is unbounded
    for i = 1:nbObj
        JuMP.set_objective(m, objSenses[i], objs[i])
        JuMP.optimize!(m, ignore_optimize_hook=true)
        status = JuMP.termination_status(m)
        if status != MOI.OPTIMAL
            return status
        end
    end

    #Create the constraints on the objectives
    cstr_rhs = [JuMP.@variable(m) for _ = 1:nbObj]
    cstr_obj = [objSenses[i] == MOI.MAX_SENSE ? JuMP.@constraint(m, objs[i] >= cstr_rhs[i]) : JuMP.@constraint(m, objs[i] <= cstr_rhs[i]) for i=1:nbObj]

    for p in Combinatorics.permutations(1:nbObj, nbObj)
        verbose && println("solving for objectives $p")
        solve_permutation(m, p, cstr_obj, cstr_rhs ; args...)
    end

    JuMP.delete.(m, cstr_obj)
    JuMP.delete.(m, cstr_rhs)
   
    return MOI.OPTIMAL
end

function solve_permutation(m::JuMP.Model, p, cstr_obj, cstr_rhs ; args...)

    vd = getvOptData(m)
    objs = vd.objs
    objSenses = vd.objSenses

    #Set the first objective of the permutation as an objective in the JuMP JuMP.Model
    JuMP.set_objective(m, objSenses[first(p)], objs[first(p)])

    #Solve with that objective
    JuMP.optimize!(m, ignore_optimize_hook=true)
    
    status = JuMP.termination_status(m)
    if status != MOI.OPTIMAL
        return status
    end

    for i = 2:length(p)
        fVal = JuMP.value(objs[p[i-1]]) #get the value for the last objective solved
        slack = objSenses[p[i-1]] == MOI.MAX_SENSE ? -1e-8 : 1e-8
        # JuMP.fix(cstr_rhs[p[i-1]], fVal - objs[p[i-1]].constant + slack)
        JuMP.fix(cstr_rhs[p[i-1]], fVal + slack)
        JuMP.set_objective(m, objSenses[p[i]], objs[p[i]]) #set the i-th objective of the permutation in the JuMP JuMP.Model
        JuMP.optimize!(m, ignore_optimize_hook=true) #and solve
    end    

    varArray = JuMP.all_variables(m)
    #Store results in vOptData
    push!(vd.Y_N, map(JuMP.value, objs))
    push!(vd.X_E, JuMP.value.(varArray))

    for var in cstr_rhs
        JuMP.is_fixed(var) && JuMP.unfix.(var) #errors if calling unfix on an unfixed variable
    end

end

function solve_eps(m::JuMP.Model, ϵ::Float64, round_results, verbose ; args...)
    @info "ϵ = $ϵ"
    #Retrieve objectives and their senses from vOptData
    vd = getvOptData(m)
    empty!(vd.Y_N) ; empty!(vd.X_E)# ; empty!(vd.logObjs)
    f1,f2 = vd.objs
    f1Sense,f2Sense = vd.objSenses
    varArray = JuMP.all_variables(m)
    
    #Set the first objective as an objective in the JuMP Model
    JuMP.set_objective(m, f1Sense, f1)
    
    R1 = f1Sense==MOI.MIN_SENSE ? (<=) : (>=)
    R2 = f2Sense==MOI.MIN_SENSE ? (<=) : (>=)
    weak_dom(a, b) = R1(a[1], b[1]) && R2(a[2], b[2])

    #Declare the epsilon-constraint (RHS will be set later)
    RHS = JuMP.@variable(m)
    eps = (f2Sense == MOI.MIN_SENSE) ? JuMP.@constraint(m, f2 <= RHS) : JuMP.@constraint(m, f2 >= RHS)

    #Solve with that objective
    JuMP.optimize!(m, ignore_optimize_hook=true)
    status = JuMP.termination_status(m)

    if status == MOI.OPTIMAL
        #While a solution exists
        while status == MOI.OPTIMAL
            #Get the score on the objectives
            f1Val = JuMP.value(f1)
            f2Val = JuMP.value(f2)
    
            #If last solution found is dominated by this one
            if length(vd.Y_N) > 0
                if weak_dom((f1Val, f2Val), vd.Y_N[end])
                    verbose && println(vd.Y_N[end], "dominated by ($f1Val, $f2Val)")
                    pop!(vd.Y_N) ; pop!(vd.X_E) #Remove last solution from Y_N and X_E
                end
            end

            #Store results in vOptData
            push!(vd.X_E, JuMP.value.(varArray))
            push!(vd.Y_N, [f1Val, f2Val])

            verbose && print("z1 = ", f1Val, ", z2 = ", f2Val)

            #Set the RHS of the epsilon-constraint
            if f2Sense == MOI.MIN_SENSE
                JuMP.fix(RHS, f2Val - ϵ)
                verbose && println(". Solving with f2 <= ", f2Val - ϵ)
            else
                JuMP.fix(RHS, f2Val + ϵ)
                verbose && println(". Solving with f2 >= ", f2Val + ϵ)
            end

            #And solve again
            JuMP.optimize!(m, ignore_optimize_hook=true)
            status = JuMP.termination_status(m)
        end

        #Sort X_E and Y_N
        s = sortperm(vd.Y_N, by = first)
        vd.X_E, vd.Y_N = vd.X_E[s], vd.Y_N[s]

        #Remove the epsilon constraint
        JuMP.delete(m, eps) ; JuMP.delete(m, RHS)
    else
        return status
    end
    return MOI.OPTIMAL
end

global varArray = Array{JuMP.VariableRef}
global x_star = []
global model 

function callback_noCuts(cb_data::CPLEX.CallbackContext, context_id) 
    global varArray
    global x_star
    global model 

    node_statut = callback_node_status(cb_data, model)
     
    if node_statut == MOI.CALLBACK_NODE_STATUS_FRACTIONAL || node_statut == MOI.CALLBACK_NODE_STATUS_INTEGER
        CPLEX.load_callback_variable_primal(cb_data, context_id)
        x_star = callback_value.(cb_data, varArray)
    end
end

function stock_all_primal_sols(m::JuMP.Model, f1, f2, varArray)
    Y_integer = Vector{Vector{Float64}}() ; X_integer = Vector{Vector{Float64}}()
    sols_nb = result_count(m)
    for i =1:sols_nb 
        if JuMP.has_values(m, result = Int(i))
            push!(Y_integer, [JuMP.value(f1, result = Int(i)), JuMP.value(f2, result = Int(i))])
            push!(X_integer, JuMP.value.(varArray, result = Int(i)))
        end
    end
    return Y_integer, X_integer
end

function solve_dicho_callback(m::JuMP.Model, round_results, verbose; args...)
    global varArray
    global x_star
    global model = m

    vd = getvOptData(m)
    empty!(vd.Y_N) ; empty!(vd.X_E) ; empty!(vd.lambda) # ; empty!(vd.logObjs)
    f1, f2 = vd.objs
    f1Sense, f2Sense = vd.objSenses
    varArray = JuMP.all_variables(m)

    Y_integer = Vector{Vector{Float64}}() ; X_integer = Vector{Vector{Float64}}()

    R1 = f1Sense==MOI.MIN_SENSE ? (<=) : (>=)
    R2 = f2Sense==MOI.MIN_SENSE ? (<=) : (>=)
    weak_dom(a, b) = R1(a[1], b[1]) && R2(a[2], b[2]) && a[1]!= b[1] && a[2]!= b[2]

    function sorting()
        # #Sort X_E and Y_N
        s = sortperm(vd.Y_N, by = x -> (-x[1], x[2]))
        vd.Y_N = vd.Y_N[s] ; vd.X_E = vd.X_E[s]
        vd.lambda = vd.lambda[s]

        i = 1
        while i < length(vd.Y_N)
            j = i+1
            while j<= length(vd.Y_N)
                if weak_dom(vd.Y_N[i], vd.Y_N[j])
                    # if vd.Y_N[j] == [-101062.0, -315795.0]
                    #     error( "dicho deleting $(vd.Y_N[j]) !")
                    # end

                    deleteat!(vd.Y_N, j) ; deleteat!(vd.X_E, j)
                    deleteat!(vd.lambda, j)
                elseif weak_dom(vd.Y_N[j], vd.Y_N[i])
                    # if vd.Y_N[i] == [-101062.0, -315795.0]
                    #     error( "dicho deleting $(vd.Y_N[i]) !")
                    # end

                    deleteat!(vd.Y_N, i) ; deleteat!(vd.X_E, i)
                    deleteat!(vd.lambda, i)
                    j -= 1 ; break
                else
                    j += 1
                end
            end
            if j > length(vd.Y_N) i += 1 end 
        end

        #Sort X_integer and Y_integer
        s = sortperm(Y_integer, by = x -> (-x[1], x[2]))
        Y_integer = Y_integer[s] ; X_integer = X_integer[s]

        i = 1
        while i < length(Y_integer)
            j = i+1
            while j<= length(Y_integer)
                if weak_dom(Y_integer[i], Y_integer[j])
                    # if Y_integer[j] == [-101062.0, -315795.0]
                    #     error("dicho Y_integer deleting $(Y_integer[j]) !")
                    # end

                    deleteat!(Y_integer, j) ; deleteat!(X_integer, j)
                elseif weak_dom(Y_integer[j], Y_integer[i])
                    # if Y_integer[i]== [-101062.0, -315795.0]
                    #     error( "dicho Y_integer deleting $(Y_integer[i]) !")
                    # end

                    deleteat!(Y_integer, i) ; deleteat!(X_integer, i)
                    j -= 1 ; break
                else
                    j += 1
                end
            end
            if j > length(Y_integer) i += 1 end 
        end
    end

    #Set the first objective as an objective in the JuMP JuMP.Model
    JuMP.set_objective(m, f1Sense, f1)
    verbose && println("solving for z1")
    
    #Solve with that objective
    x_star = []
    MOI.set(m, MOI.NumberOfThreads(), 1)
    MOI.set(model, CPLEX.CallbackFunction(), callback_noCuts)
    JuMP.optimize!(m, ignore_optimize_hook=true)
    status = JuMP.termination_status(m)

    #If a solution exists
    yr_1 = 0.0 ; yr_2 = 0.0 
    # @info "                 status = $status "
    if status == MOI.INFEASIBLE 
        sorting()
        return Y_integer, X_integer 
    end
    if status == MOI.OPTIMAL
        # stock heur sol 
        Y, X = stock_all_primal_sols(m, f1, f2, varArray)
        append!(Y_integer, Y) ; append!(X_integer, X)

        yr_1 = JuMP.value(f1)
        yr_2 = JuMP.value(f2)
        #Store results in vOptData
        push!(vd.Y_N, round_results ? round.([yr_1, yr_2]) : [yr_1, yr_2])
        push!(vd.X_E, JuMP.value.(varArray))
        push!(vd.lambda, [1.0, 0.0])

    elseif status == MOI.NODE_LIMIT && length(x_star) > 0
        if has_values(m)
            # stock heur sol 
            Y, X = stock_all_primal_sols(m, f1, f2, varArray)
            append!(Y_integer, Y) ; append!(X_integer, X)
        end
        yr_1 = x_star'* collect(values(f1.terms)) + f1.constant
        yr_2 = x_star'* collect(values(f2.terms)) + f2.constant
        #Store results in vOptData
        push!(vd.Y_N, round_results ? round.([yr_1, yr_2]) : [yr_1, yr_2])
        push!(vd.X_E, x_star)
        push!(vd.lambda, [1.0, 0.0])

    else
        println("has primal ? $(JuMP.has_values(m))")
        error("Condition  status $status ")
    end

    #Set the second objective as an objective in the JuMP JuMP.Model
    JuMP.set_objective(m, f2Sense, f2)
    verbose && println("solving for z2")

    #Solve with that objective
    x_star = []
    JuMP.optimize!(m, ignore_optimize_hook=true)
    status = JuMP.termination_status(m)

    ys_1 = 0.0 ; ys_2 = 0.0 
    # @info "                 status = $status "

    if status == MOI.INFEASIBLE
        sorting()
        return Y_integer, X_integer
    end
    if status == MOI.OPTIMAL
        # stock heur sol 
        Y, X = stock_all_primal_sols(m, f1, f2, varArray)
        append!(Y_integer, Y) ; append!(X_integer, X)

        ys_1 = JuMP.value(f1)
        ys_2 = JuMP.value(f2)

        if !isapprox(yr_1, ys_1, atol=1e-3) || !isapprox(yr_2, ys_2, atol=1e-3)
            #Store results in vOptData
            push!(vd.Y_N, round_results ? round.([ys_1, ys_2]) : [ys_1, ys_2])
            push!(vd.X_E, JuMP.value.(varArray))
            push!(vd.lambda, [0.0, 1.0])

            Y, X = dichoRecursion_callback(m, yr_1, yr_2, ys_1, ys_2, varArray, round_results, verbose ; args...)
            append!(Y_integer, Y) ; append!(X_integer, X)
        end

    elseif status == MOI.NODE_LIMIT && length(x_star) > 0
        if has_values(m)
            # stock heur sol 
            Y, X = stock_all_primal_sols(m, f1, f2, varArray)
            append!(Y_integer, Y) ; append!(X_integer, X)
        end
        ys_1 = x_star'* collect(values(f1.terms)) + f1.constant
        ys_2 = x_star'* collect(values(f2.terms)) + f2.constant

        if !isapprox(yr_1, ys_1, atol=1e-3) || !isapprox(yr_2, ys_2, atol=1e-3)
            #Store results in vOptData
            push!(vd.Y_N, round_results ? round.([ys_1, ys_2]) : [ys_1, ys_2])
            push!(vd.X_E, x_star)
            push!(vd.lambda, [0.0, 1.0])

            Y, X = dichoRecursion_callback(m, yr_1, yr_2, ys_1, ys_2, varArray, round_results, verbose ; args...)
            append!(Y_integer, Y) ; append!(X_integer, X)
        end

    else
        println("has primal ? $(JuMP.has_values(m))")
        error("Condition  status $status ")
    end

    # println("\n\n----------- \n before : ", vd.Y_N)
    # println("\n\n----------- \n before : ", Y_integer)
    sorting()
    # ----------------------
    #todo: verifying
    for i = 1:length(vd.Y_N)-1
        for j = i+1:length(vd.Y_N)
            if weak_dom(vd.Y_N[i], vd.Y_N[j]) || weak_dom(vd.Y_N[j], vd.Y_N[i])
                println("----------- \n after : ", vd.Y_N)
                error("vd.Y_N error in dicho ! ")
            end

            if vd.Y_N[i][1] < vd.Y_N[j][1] || vd.Y_N[i][2] > vd.Y_N[j][2]
                error("NATURAL ORDER vd.Y_N error in dicho")
            end
        end
    end

    for i = 1:length(Y_integer)-1
        for j = i+1:length(Y_integer)
            if weak_dom(Y_integer[i], Y_integer[j]) || weak_dom(Y_integer[j], Y_integer[i])
                println("----------- \n after : ", Y_integer)
                error("Y_integer error in dicho ! ")
            end

            if Y_integer[i][1] < Y_integer[j][1] || Y_integer[i][2] > Y_integer[j][2]
                error("NATURAL ORDER Y_integer error in dicho")
            end
        end
    end
    return Y_integer, X_integer
end

function dichoRecursion_callback(m::JuMP.Model, yr_1, yr_2, ys_1, ys_2, varArray, round_results, verbose ; args...)
    global varArray
    global x_star

    Y_integer = Vector{Vector{Float64}}() ; X_integer = Vector{Vector{Float64}}()

    vd = getvOptData(m)
    f1, f2 = vd.objs
    f1Sense, f2Sense = vd.objSenses

    λ1 = abs(yr_2 - ys_2)
    λ2 = abs(ys_1 - yr_1)

    f = AffExpr(0.0)

    if f1Sense==f2Sense
        lb = λ1*yr_1 + λ2*yr_2
        JuMP.set_objective(m, f1Sense, λ1*f1 + λ2*f2)
        verbose && println("solving for $λ1*f1 + $λ2*f2")    
        f = λ1*f1 + λ2*f2
    else
        lb = λ1*yr_1 - λ2*yr_2
        JuMP.set_objective(m, f1Sense, λ1*f1 - λ2*f2)
        verbose && println("solving for $λ1*f1 - $λ2*f2") 
        f = λ1*f1 - λ2*f2
    end

    x_star = []
    JuMP.optimize!(m, ignore_optimize_hook=true)
    status = JuMP.termination_status(m)

    #If a solution exists
    yt_1 = 0.0 ; yt_2 = 0.0 
    # @info "                 status = $status "

    if status == MOI.INFEASIBLE return Y_integer, X_integer end
    if status == MOI.OPTIMAL 
        # stock heur sol 
        Y, X = stock_all_primal_sols(m, f1, f2, varArray)
        append!(Y_integer, Y) ; append!(X_integer, X)

        yt_1 = JuMP.value(f1)
        yt_2 = JuMP.value(f2)
        val = f1Sense == f2Sense ? λ1*yt_1 + λ2*yt_2 : λ1*yt_1 - λ2*yt_2

        if (val < lb - 1e-4)
            push!(vd.Y_N, round_results ? round.([yt_1, yt_2]) : [yt_1, yt_2])
            push!(vd.X_E, JuMP.value.(varArray))
            push!(vd.lambda, [λ1, λ2])

            # push!(vd.logObjs, f)
            if yt_1 > ys_1+1e-4 && yt_2 > yr_2+1e-4
                Y, X = dichoRecursion_callback(m, yr_1, yr_2, yt_1, yt_2, varArray, round_results, verbose ; args...)
                append!(Y_integer, Y) ; append!(X_integer, X)
                Y, X = dichoRecursion_callback(m, yt_1, yt_2, ys_1, ys_2, varArray, round_results, verbose ; args...)
                append!(Y_integer, Y) ; append!(X_integer, X)
            end
        end

    elseif status == MOI.NODE_LIMIT && length(x_star) > 0
        if has_values(m)
            # stock heur sol 
            Y, X = stock_all_primal_sols(m, f1, f2, varArray)
            append!(Y_integer, Y) ; append!(X_integer, X)
        end
        yt_1 = x_star'* collect(values(f1.terms)) + f1.constant
        yt_2 = x_star'* collect(values(f2.terms)) + f2.constant
        val = f1Sense == f2Sense ? λ1*yt_1 + λ2*yt_2 : λ1*yt_1 - λ2*yt_2

        if ( val < lb - 1e-4)
            push!(vd.Y_N, round_results ? round.([yt_1, yt_2]) : [yt_1, yt_2])
            push!(vd.X_E, x_star)
            push!(vd.lambda, [λ1, λ2])

            # push!(vd.logObjs, f)
            if yt_1 > ys_1 +1e-4 && yt_2 > yr_2 +1e-4 
                Y, X = dichoRecursion_callback(m, yr_1, yr_2, yt_1, yt_2, varArray, round_results, verbose ; args...)
                append!(Y_integer, Y) ; append!(X_integer, X)
                Y, X = dichoRecursion_callback(m, yt_1, yt_2, ys_1, ys_2, varArray, round_results, verbose ; args...)
                append!(Y_integer, Y) ; append!(X_integer, X)
            end
        end

    else
        println("has primal ? $(JuMP.has_values(m))")

        error("Condition  status $status ")
    end
    return Y_integer, X_integer

end

function solve_Chalmet(m::JuMP.Model, step, verbose ; args...)
    vd = getvOptData(m)
    empty!(vd.Y_N) ; empty!(vd.X_E)# ; empty!(vd.logObjs)
    f1,f2 = vd.objs[1],vd.objs[2]
    f1Sense,f2Sense = vd.objSenses[1],vd.objSenses[2]
    varArray = JuMP.all_variables(m)

    #Set the first objective as an objective in the JuMP JuMP.Model
    JuMP.set_objective(m, f1Sense, f1)
    verbose && println("solving for z1")

    #Solve with that objective
    JuMP.optimize!(m, ignore_optimize_hook=true)
    status = JuMP.termination_status(m)

    #If a solution exists
    if status == MOI.OPTIMAL

        yr_1 = JuMP.value(f1)
        yr_2 = JuMP.value(f2)

        #Store results in vOptData
        push!(vd.Y_N, [yr_1, yr_2])
        push!(vd.X_E, JuMP.value.(varArray))

        #Set the second objective as an objective in the JuMP JuMP.Model
        JuMP.set_objective(m, f2Sense, f2)
        verbose && println("solving for z2")


        #Solve with that objective
        JuMP.optimize!(m, ignore_optimize_hook=true)
        status = JuMP.termination_status(m)

        if status == MOI.OPTIMAL

            ys_1 = JuMP.value(f1)
            ys_2 = JuMP.value(f2)

            if !isapprox(yr_1, ys_1, atol=1e-3) || !isapprox(yr_2, ys_2, atol=1e-3)
                push!(vd.Y_N, [ys_1, ys_2])
                push!(vd.X_E, JuMP.value.(varArray))
                #Declare the constraints on z1 and z2 (RHS will be set later)
                rhs_z1 = JuMP.@variable(m)
                rhs_z2 = JuMP.@variable(m)
                cstr_z1 = f1Sense == MOI.MIN_SENSE ? JuMP.@constraint(m, f1 <= rhs_z1) : JuMP.@constraint(m, f1 >= rhs_z1)
                cstr_z2 = f2Sense == MOI.MIN_SENSE ? JuMP.@constraint(m, f2 <= rhs_z2) : JuMP.@constraint(m, f2 >= rhs_z2)

                ChalmetRecursion(m, yr_1, yr_2, ys_1, ys_2, varArray, rhs_z1, rhs_z2, step, verbose; args...)

                #Sort X_E and Y_N
                s = sortperm(vd.Y_N, by = first)
                vd.Y_N = vd.Y_N[s]
                vd.X_E = vd.X_E[s]

                R1 = f1Sense==MOI.MIN_SENSE ? (<=) : (>=)
                R2 = f2Sense==MOI.MIN_SENSE ? (<=) : (>=)
                weak_dom(a, b) = R1(a[1], b[1]) && R2(a[2], b[2])

                #Filter X_E and Y_N :
                inds = Int[]
                for i = 1:length(vd.Y_N)-1
                    if weak_dom(vd.Y_N[i], vd.Y_N[i+1])
                        push!(inds, i+1)
                    elseif weak_dom(vd.Y_N[i+1], vd.Y_N[i])
                        push!(inds, i)
                    end
                end
                deleteat!(vd.Y_N, inds)
                deleteat!(vd.X_E, inds)

                JuMP.delete(m, cstr_z1)
                JuMP.delete(m, cstr_z2)
                JuMP.delete(m, rhs_z1)
                JuMP.delete(m, rhs_z2)
            end
        end
    end

    return MOI.OPTIMAL
end

function ChalmetRecursion(m::JuMP.Model, yr_1, yr_2, ys_1, ys_2, varArray, rhs_z1, rhs_z2, ϵ, verbose ; args...)

    vd = getvOptData(m)
    f1, f2 = vd.objs
    f1Sense, f2Sense = vd.objSenses
    JuMP.set_objective(m, f1Sense, f1Sense==f2Sense ? f1 + f2 : f1 - f2)
    lbz1 = f1Sense==MOI.MAX_SENSE ? min(yr_1, ys_1) : max(yr_1, ys_1)
    lbz2 = f2Sense==MOI.MAX_SENSE ? min(yr_2, ys_2) : max(yr_2, ys_2)

    if f1Sense == MOI.MIN_SENSE
        JuMP.fix(rhs_z1, lbz1 - ϵ)
        verbose && print("solving with z1 <= ", lbz1 - ϵ)
    else
        JuMP.fix(rhs_z1, lbz1 + ϵ)
        verbose && print("solving with z1 >= ", lbz1 + ϵ)
    end

    if f2Sense == MOI.MIN_SENSE
        JuMP.fix(rhs_z2, lbz2 - ϵ)
        verbose && println(" and z2 <= ", lbz2 - ϵ)
    else
        JuMP.fix(rhs_z2, lbz2 + ϵ)
        verbose && println(" and z2 >= ", lbz2 + ϵ)
    end

    JuMP.optimize!(m, ignore_optimize_hook=true)
    status = JuMP.termination_status(m)

    if status == MOI.OPTIMAL
        yt_1 = JuMP.value(f1)
        yt_2 = JuMP.value(f2)
        push!(vd.Y_N, [yt_1, yt_2])
        push!(vd.X_E, JuMP.value.(varArray))
        ChalmetRecursion(m, yr_1, yr_2, yt_1, yt_2, varArray, rhs_z1, rhs_z2, ϵ, verbose; args...)
        ChalmetRecursion(m, yt_1, yt_2, ys_1, ys_2, varArray, rhs_z1, rhs_z2, ϵ, verbose ; args...)
    end
end



function solve_dicho(m::JuMP.Model, round_results, verbose; args...)
    vd = getvOptData(m)
    empty!(vd.Y_N) ; empty!(vd.X_E); empty!(vd.lambda)# ; empty!(vd.logObjs)
    f1, f2 = vd.objs
    f1Sense, f2Sense = vd.objSenses
    varArray = JuMP.all_variables(m)

    #Set the first objective as an objective in the JuMP JuMP.Model
    # MOI.set(m, MOI.NumberOfThreads(), 1)
    JuMP.set_objective(m, f1Sense, f1)
    verbose && println("solving for z1")
    
    #Solve with that objective
    JuMP.optimize!(m, ignore_optimize_hook=true)
    status = JuMP.termination_status(m)

    #If a solution exists
    if status == MOI.OPTIMAL

        yr_1 = JuMP.value(f1)
        yr_2 = JuMP.value(f2)

        #Store results in vOptData
        push!(vd.Y_N, round_results ? round.([yr_1, yr_2]) : [yr_1, yr_2])
        push!(vd.X_E, JuMP.value.(varArray))
        # push!(vd.logObjs, f1)

        #Set the second objective as an objective in the JuMP JuMP.Model
        JuMP.set_objective(m, f2Sense, f2)
        verbose && println("solving for z2")

        #Solve with that objective
        JuMP.optimize!(m, ignore_optimize_hook=true)
        status = JuMP.termination_status(m)

        if status == MOI.OPTIMAL

            ys_1 = JuMP.value(f1)
            ys_2 = JuMP.value(f2)

            if !isapprox(yr_1, ys_1, atol=1e-3) || !isapprox(yr_2, ys_2, atol=1e-3)
                push!(vd.Y_N, round_results ? round.([ys_1, ys_2]) : [ys_1, ys_2])
                push!(vd.X_E, JuMP.value.(varArray))
                # push!(vd.logObjs, f2)
                dichoRecursion(m, yr_1, yr_2, ys_1, ys_2, varArray, round_results, verbose ; args...)
            end
        
            #Sort X_E and Y_N
            s = sortperm(vd.Y_N, by = first)
            vd.Y_N = vd.Y_N[s]
            vd.X_E = vd.X_E[s]
            # vd.logObjs = vd.logObjs[s]

            R1 = f1Sense==MOI.MIN_SENSE ? (<=) : (>=)
            R2 = f2Sense==MOI.MIN_SENSE ? (<=) : (>=)
            weak_dom(a, b) = R1(a[1], b[1]) && R2(a[2], b[2])

            #Filter X_E and Y_N :
            inds = Int[]
            for i = 1:length(vd.Y_N)-1
                if weak_dom(vd.Y_N[i], vd.Y_N[i+1])
                    push!(inds, i+1)
                elseif weak_dom(vd.Y_N[i+1], vd.Y_N[i])
                    push!(inds, i)
                end
            end
            deleteat!(vd.Y_N, inds)
            deleteat!(vd.X_E, inds)
            # deleteat!(vd.logObjs, inds)
        end
    end

    status
end

function dichoRecursion(m::JuMP.Model, yr_1, yr_2, ys_1, ys_2, varArray, round_results, verbose ; args...)

    vd = getvOptData(m)
    f1, f2 = vd.objs
    f1Sense, f2Sense = vd.objSenses

    λ1 = abs(yr_2 - ys_2)
    λ2 = abs(ys_1 - yr_1)

    f = AffExpr(0.0)

    if f1Sense==f2Sense
        lb = λ1*yr_1 + λ2*yr_2
        JuMP.set_objective(m, f1Sense, λ1*f1 + λ2*f2)
        verbose && println("solving for $λ1*f1 + $λ2*f2")    
        f = λ1*f1 + λ2*f2
    else
        lb = λ1*yr_1 - λ2*yr_2
        JuMP.set_objective(m, f1Sense, λ1*f1 - λ2*f2)
        verbose && println("solving for $λ1*f1 - $λ2*f2") 
        f = λ1*f1 - λ2*f2
    end

    JuMP.optimize!(m, ignore_optimize_hook=true)

    yt_1 = JuMP.value(f1)
    yt_2 = JuMP.value(f2)

    val = f1Sense == f2Sense ? λ1*yt_1 + λ2*yt_2 : λ1*yt_1 - λ2*yt_2

    if (f1Sense == MOI.MIN_SENSE && val < lb - 1e-4) || val > lb + 1e-4
        push!(vd.Y_N, round_results ? round.([yt_1, yt_2]) : [yt_1, yt_2])
        push!(vd.X_E, JuMP.value.(varArray))
        # push!(vd.logObjs, f)
        dichoRecursion(m, yr_1, yr_2, yt_1, yt_2, varArray, round_results, verbose ; args...)
        dichoRecursion(m, yt_1, yt_2, ys_1, ys_2, varArray, round_results, verbose ; args...)
    end

end