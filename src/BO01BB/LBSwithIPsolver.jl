using JuMP
include("struct.jl")


global varArray = Array{JuMP.VariableRef}
global x_star = []
global model 
global bst_val = -Inf 
global curr_λ = []
global C 

function callback_noCuts(cb_data)
    global varArray
    global x_star
    global model 
    global C
    global bst_val
    global curr_λ

    node_statut = callback_node_status(cb_data, model)
     
    # if node_statut == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
        tmp = callback_value.(cb_data, varArray)
        val = tmp'* C[1, 2:end] + C[1, 1] + tmp'* C[2, 2:end] + C[2, 1]
        if val > bst_val
            bst_val = val ; x_star = callback_value.(cb_data, varArray)
        end
    # end
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

"""
New correction iterative algorithm for LBS 
"""
function LBSinvokingIPsolveer(m::JuMP.Model, lp_copied::JuMP.Model, c, round_results, verbose; args...)
    global varArray
    global x_star
    global model = m
    global C = c
    global curr_λ

    L = RelaxedBoundSet()
    vd = getvOptData(m)
    empty!(vd.Y_N) ; empty!(vd.X_E) ; empty!(vd.lambda) 
    f1, f2 = vd.objs
    f1Sense, f2Sense = vd.objSenses
    varArray = JuMP.all_variables(m)
    Gap = 0.0

    varArray_copied = JuMP.all_variables(lp_copied)
    f1_copied = varArray_copied'* c[1, 2:end] + c[1, 1]
    f2_copied = varArray_copied'* c[2, 2:end] + c[2, 1]

    Y_integer = Vector{Vector{Float64}}() ; X_integer = Vector{Vector{Float64}}()

    R1 = f1Sense==MOI.MIN_SENSE ? (<=) : (>=)
    R2 = f2Sense==MOI.MIN_SENSE ? (<=) : (>=)
    weak_dom(a, b) = R1(a[1], b[1]) && R2(a[2], b[2]) && a[1]!= b[1] && a[2]!= b[2]

    # set up callback 
    MOI.set(m, MOI.NumberOfThreads(), 1) ; MOI.set(m, MOI.UserCutCallback(), callback_noCuts)

    # -------------------------------------------
    # step 1 : calculate the left extreme point
    # ------------------------------------------- 
    JuMP.set_objective(m, f1Sense, f1) ; JuMP.set_objective(lp_copied, f1Sense, f1_copied)
    verbose && println("solving for z1")
    
    
    x_star = [] ; bst_val = -Inf 
    curr_λ = [1.0, 0.0]
    JuMP.optimize!(m, ignore_optimize_hook=true) ; status = JuMP.termination_status(m)

    # in case of infeasibility 
    yr_1 = 0.0 ; yr_2 = 0.0 
    if status == MOI.INFEASIBLE 
        return Y_integer, X_integer, Gap #todo : vd.Y_N empty 
    end

    # in case of optimality 
    if status == MOI.OPTIMAL
        # stock heur sol 
        Y, X = stock_all_primal_sols(m, f1, f2, varArray)
        append!(Y_integer, Y) ; append!(X_integer, X)

        yr_1 = JuMP.value(f1) ; yr_2 = JuMP.value(f2)
        # store results in vOptData # todo : use natural order list 
        # push!(vd.Y_N, round_results ? round.([yr_1, yr_2]) : [yr_1, yr_2])
        # push!(vd.X_E, JuMP.value.(varArray))
        # push!(vd.lambda, curr_λ)
        push!(L.natural_order_vect, Solution(JuMP.value.(varArray), [yr_1, yr_2], curr_λ ), filtered=true)

    # otherwise, take the best primal sol so far 
    elseif status == MOI.NODE_LIMIT || status == TIME_LIMIT
        # stock heuristic sol 
        if has_values(m)
            Y, X = stock_all_primal_sols(m, f1, f2, varArray) ; append!(Y_integer, Y) ; append!(X_integer, X)
        end

        if length(x_star) > 0
            nothing
        else
            best_bound = objective_bound(m)
            ctr_bound = JuMP.@constraint(lp_copied, f1_copied >= best_bound)
            JuMP.optimize!(lp_copied, ignore_optimize_hook=true)
            x_star = JuMP.value.(varArray_copied)

            if JuMP.is_valid(lp_copied, ctr_bound)
                JuMP.delete(lp_copied, ctr_bound) ; JuMP.unregister(lp_copied, :ctr_bound)
            end
        end

        yr_1 = x_star'* c[1, 2:end] + c[1, 1]
        yr_2 = x_star'* c[2, 2:end] + c[2, 1]
        # store results in vOptData # todo : use natural order list 
        # push!(vd.Y_N, round_results ? round.([yr_1, yr_2]) : [yr_1, yr_2])
        # push!(vd.X_E, x_star) ; push!(vd.lambda, curr_λ)

        push!(L.natural_order_vect, Solution(x_star, [yr_1, yr_2], curr_λ ), filerted=true)
    else
        println("has primal ? $(JuMP.has_values(m))")
        error("Condition  status $status ")
    end

    # -------------------------------------------
    # step 2 : calculate the right extreme point
    # ------------------------------------------- 
    JuMP.set_objective(m, f2Sense, f2) ; JuMP.set_objective(lp_copied, f2Sense, f2_copied)
    verbose && println("solving for z2")

    x_star = [] ; bst_val = -Inf 
    curr_λ = [0.0, 1.0]
    JuMP.optimize!(m, ignore_optimize_hook=true) ; status = JuMP.termination_status(m)

    ys_1 = 0.0 ; ys_2 = 0.0 
    if status == MOI.INFEASIBLE
        return Y_integer, X_integer, Gap #todo : vd.Y_N empty
    end

    if status == MOI.OPTIMAL
        # stock heur sol 
        Y, X = stock_all_primal_sols(m, f1, f2, varArray)
        append!(Y_integer, Y) ; append!(X_integer, X)

        #Store results in vOptData # todo : use natural order list 
        push!(vd.Y_N, round_results ? round.([ys_1, ys_2]) : [ys_1, ys_2])
        push!(vd.X_E, JuMP.value.(varArray))
        push!(vd.lambda, curr_λ)

        ys_1 = JuMP.value(f1) ; ys_2 = JuMP.value(f2)
        if isapprox(yr_1, ys_1, atol=1e-3) || isapprox(yr_2, ys_2, atol=1e-3)


            # # todo : stop 
            return Y_integer, X_integer, Gap #todo : vd.Y_N empty
        end

    elseif status == MOI.NODE_LIMIT || status == TIME_LIMIT
        if has_values(m)
            Y, X = stock_all_primal_sols(m, f1, f2, varArray)
            append!(Y_integer, Y) ; append!(X_integer, X)
        end

        if length(x_star) > 0
            nothing
        else
            best_bound = objective_bound(m)
            ctr_bound = JuMP.@constraint(lp_copied, f2_copied >= best_bound)
            JuMP.optimize!(lp_copied, ignore_optimize_hook=true)

            x_star = JuMP.value.(varArray_copied)

            if JuMP.is_valid(lp_copied, ctr_bound)
                JuMP.delete(lp_copied, ctr_bound) ; JuMP.unregister(lp_copied, :ctr_bound)
            end
        end

        ys_1 = x_star'* c[1, 2:end] + c[1, 1]
        ys_2 = x_star'* c[2, 2:end] + c[2, 1]

        # store results in vOptData # todo : use natural order list 
        # push!(vd.Y_N, round_results ? round.([ys_1, ys_2]) : [ys_1, ys_2])
        # push!(vd.X_E, x_star)
        # push!(vd.lambda, curr_λ)
        push!(L.natural_order_vect, Solution(x_star, [ys_1, ys_2], curr_λ ), filtered=true)

        if isapprox(yr_1, ys_1, atol=1e-3) || isapprox(yr_2, ys_2, atol=1e-3)
            # # todo : stop 
            return Y_integer, X_integer, Gap #todo : vd.Y_N empty
        end
    else
        println("has primal ? $(JuMP.has_values(m))")
        error("Condition  status $status ")
    end

    # -----------------------------------------
    # step 3 : fix the next search direction 
    # -----------------------------------------
    #todo : using natural order list 

    return Y_integer, X_integer, Gap #todo : vd.Y_N empty
end



function dichoRecursion_callback(m::JuMP.Model, lp_copied::JuMP.Model, c, yr_1, yr_2, ys_1, ys_2, varArray, round_results, verbose ; args...)
    global varArray
    global x_star
    global curr_λ


    Y_integer = Vector{Vector{Float64}}() ; X_integer = Vector{Vector{Float64}}()

    vd = getvOptData(m)
    f1, f2 = vd.objs
    f1Sense, f2Sense = vd.objSenses

    varArray_copied = JuMP.all_variables(lp_copied)

    f1_copied = varArray_copied'* c[1, 2:end] + c[1, 1]
    f2_copied = varArray_copied'* c[2, 2:end] + c[2, 1]

    λ1 = abs(yr_2 - ys_2)
    λ2 = abs(ys_1 - yr_1)

    f = AffExpr(0.0)
    Gap = 0.0

    lb = λ1*yr_1 + λ2*yr_2
    JuMP.set_objective(m, f1Sense, λ1*f1 + λ2*f2) ; JuMP.set_objective(lp_copied, f1Sense, λ1*f1_copied + λ2*f2_copied)
    verbose && println("solving for $λ1*f1 + $λ2*f2")    
    f = λ1*f1 + λ2*f2


    x_star = [] ; bst_val = -Inf 
    curr_λ = [λ1, λ2]
    JuMP.optimize!(m, ignore_optimize_hook=true) ; status = JuMP.termination_status(m)
    #If a solution exists
    yt_1 = 0.0 ; yt_2 = 0.0 

    if status == MOI.INFEASIBLE return Y_integer, X_integer, Gap end
    # Gap += MOI.get(m, MOI.RelativeGap())*100

    if status == MOI.OPTIMAL 
        # stock heur sol 
        Y, X = stock_all_primal_sols(m, f1, f2, varArray)
        append!(Y_integer, Y) ; append!(X_integer, X)

        yt_1 = JuMP.value(f1) ; yt_2 = JuMP.value(f2)
        val = f1Sense == f2Sense ? λ1*yt_1 + λ2*yt_2 : λ1*yt_1 - λ2*yt_2

        if (val < lb - 1e-4)
            push!(vd.Y_N, round_results ? round.([yt_1, yt_2]) : [yt_1, yt_2])
            push!(vd.X_E, JuMP.value.(varArray))
            push!(vd.lambda, curr_λv)

            if yt_1 > ys_1+1e-4 && yt_2 > yr_2+1e-4
                Y, X, G = dichoRecursion_callback(m, lp_copied, c, yr_1, yr_2, yt_1, yt_2, varArray, round_results, verbose ; args...)
                append!(Y_integer, Y) ; append!(X_integer, X)
                Gap += G
                Y, X, G = dichoRecursion_callback(m, lp_copied, c, yt_1, yt_2, ys_1, ys_2, varArray, round_results, verbose ; args...)
                append!(Y_integer, Y) ; append!(X_integer, X)
                Gap += G
            end
        end

    elseif status == MOI.NODE_LIMIT || status == TIME_LIMIT
        
        if has_values(m)
            Y, X = stock_all_primal_sols(m, f1, f2, varArray)
            append!(Y_integer, Y) ; append!(X_integer, X)
        end

        if length(x_star) > 0
            nothing
        else
            best_bound = objective_bound(m)
            ctr_bound = JuMP.@constraint(lp_copied, λ1*f1_copied + λ2*f2_copied >= best_bound)
            JuMP.optimize!(lp_copied, ignore_optimize_hook=true)

            x_star = JuMP.value.(varArray_copied)

            if JuMP.is_valid(lp_copied, ctr_bound)
                JuMP.delete(lp_copied, ctr_bound) ; JuMP.unregister(lp_copied, :ctr_bound)
            end
        end

        yt_1 = x_star'* c[1, 2:end] + c[1, 1]
        yt_2 = x_star'* c[2, 2:end] + c[2, 1]
        val = f1Sense == f2Sense ? λ1*yt_1 + λ2*yt_2 : λ1*yt_1 - λ2*yt_2

        if ( val < lb - 1e-4)
            push!(vd.Y_N, round_results ? round.([yt_1, yt_2]) : [yt_1, yt_2])
            push!(vd.X_E, x_star) ; push!(vd.lambda, curr_λ)

            if yt_1 > ys_1 +1e-4 && yt_2 > yr_2 +1e-4 
                Y, X, G = dichoRecursion_callback(m,lp_copied, c, yr_1, yr_2, yt_1, yt_2, varArray, round_results, verbose ; args...)
                append!(Y_integer, Y) ; append!(X_integer, X)
                Gap += G
                Y, X, G = dichoRecursion_callback(m, lp_copied, c, yt_1, yt_2, ys_1, ys_2, varArray, round_results, verbose ; args...)
                append!(Y_integer, Y) ; append!(X_integer, X)
                Gap += G
            end
        end

    else
        println("has primal ? $(JuMP.has_values(m))")
        error("Condition  status $status ")
    end
    return Y_integer, X_integer, Gap

end