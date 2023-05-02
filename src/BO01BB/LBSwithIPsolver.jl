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

    # node_statut = callback_node_status(cb_data, model)
     
    # if node_statut == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
    #     tmp = callback_value.(cb_data, varArray)
    #     val = (tmp'* C[1, 2:end] + C[1, 1]) * curr_λ[1] + (tmp'* C[2, 2:end] + C[2, 1] ) *curr_λ[2]
    #     if val > bst_val
    #         bst_val = val ; x_star = callback_value.(cb_data, varArray)
    #     end
    x_star = callback_value.(cb_data, varArray)
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
straight line λ1 z1 = λ2 z2 + ct.
"""
function updateCT(s::Solution)
    if s.ct ≠ 0.0 || length(s.λ) < 2 || s.λ == [0.0, 0.0] return end 

    s.ct = s.λ[1] * s.y[1] + s.λ[2] * s.y[2]
end

"""
Supprime all lower bounds under current straight line.

    complexity : O(|L|) # todo : improve
"""
function filtering(lb::Float64, L::RelaxedBoundSet, λ)
    # remove all points under current line 
    to_delete = Int64[] ; i = 1
    for s in L.natural_order_vect.sols
        if s.y[1]*λ[1] + s.y[2]*λ[2] < lb - 1e-4 push!(to_delete, i) end
        i += 1
    end

    deleteat!(L.natural_order_vect.sols, to_delete)
end


"""
Find the next search direction from current point idx.

    complexity : O(|L|)
"""
function next_direc(idx::Int64, L::RelaxedBoundSet, todo)
    i = idx ; j = idx ; l = 0 ; r = length(L.natural_order_vect.sols) + 1
    while i > 1
        i -= 1
        if length(L.natural_order_vect.sols[i].xEquiv) > 0 l = i ; break end 
    end

    while j < length(L.natural_order_vect.sols)
        j += 1
        if length(L.natural_order_vect.sols[j].xEquiv) > 0 r = j ; break end 
    end

    # prepare next directions 
    l ≥ 1 ? push!(todo, [L.natural_order_vect.sols[l].y, L.natural_order_vect.sols[idx].y]) : nothing 
    r ≤ length(L.natural_order_vect.sols) ? push!(todo, [L.natural_order_vect.sols[idx].y, L.natural_order_vect.sols[r].y]) : nothing 
end


"""
Computing intersected points between current point idx ans LBS.

    complexity : O(|L|^2)  #todo : improve complexity 
"""
function intersectionPts(L::RelaxedBoundSet, idx::Int)::Set{Solution}
    # @info "searching intersection pts for idx $idx "
    s2 = L.natural_order_vect.sols[idx] ; res = Set{Solution}()

    for i = 1:length(L.natural_order_vect.sols)
        if i== idx continue end 

        s1 = L.natural_order_vect.sols[i]

        det = s1.λ[1] * s2.λ[2] - s2.λ[1] * s1.λ[2]
        if abs(det) < 1e-3 continue end 
        
        y = [ (s1.ct * s2.λ[2] - s2.ct * s1.λ[2])/det ,
            (s1.λ[1] * s2.ct - s2.λ[1] * s1.ct )/det
        ]

        # println("$y  found !! ")
        # check 
        valid = true
        for j = 1:length(L.natural_order_vect.sols)
            t = L.natural_order_vect.sols[j]
            if y[1] * t.λ[1] + y[2] * t.λ[2] < t.y[1] * t.λ[1] + t.y[2] * t.λ[2] - 1e-4
                valid = false ; break
            end
        end

        if valid 
            new = Solution(Vector{Vector{Float64}}(), y, false, s1.λ, 0.0) ; updateCT(new)
            push!(res, new )
        end
    end
    
    return res
end

"""
New correction iterative algorithm for LBS 
    taking account into intersection with odd LBS
"""
function LBSinvokingIPsolveer(L::RelaxedBoundSet , m::JuMP.Model, lp_copied::JuMP.Model, c, verbose; args...)
    global varArray
    global x_star
    global model = m
    global C = c
    global curr_λ
    global bst_val

    # L = RelaxedBoundSet()
    vd = getvOptData(m)
    # empty!(vd.Y_N) ; empty!(vd.X_E) ; empty!(vd.lambda) 
    f1, f2 = vd.objs
    f1Sense, f2Sense = vd.objSenses
    varArray = JuMP.all_variables(m)

    varArray_copied = JuMP.all_variables(lp_copied)
    f1_copied = varArray_copied'* c[1, 2:end] + c[1, 1]
    f2_copied = varArray_copied'* c[2, 2:end] + c[2, 1]

    Y_integer = Vector{Vector{Float64}}() ; X_integer = Vector{Vector{Float64}}()
    
    # set up callback 
    MOI.set(m, MOI.NumberOfThreads(), 1) ; MOI.set(m, MOI.UserCutCallback(), callback_noCuts)

    # start from the extreme point 
    if length(L.natural_order_vect.sols) < 2
        println("\n\t calculate extreme pts ")
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
        val = -Inf ; idx = -1

        if status == MOI.INFEASIBLE 
            return Y_integer, X_integer
        end

        # in case of optimality 
        if status == MOI.OPTIMAL
            # stock heur sol 
            Y, X = stock_all_primal_sols(m, f1, f2, varArray)
            append!(Y_integer, Y) ; append!(X_integer, X)

            yr_1 = JuMP.value(f1) ; yr_2 = JuMP.value(f2)

            idx = push!(L.natural_order_vect, Solution(JuMP.value.(varArray), [yr_1, yr_2], curr_λ ), filtered=true )
            idx > 0 ? updateCT(L.natural_order_vect.sols[idx]) : nothing

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

            idx = push!(L.natural_order_vect, Solution(x_star, [yr_1, yr_2], curr_λ ), filtered = true)
            idx > 0 ? updateCT(L.natural_order_vect.sols[idx]) : nothing
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
        val = -Inf ; idx = -1

        if status == MOI.INFEASIBLE
            return Y_integer, X_integer
        end

        if status == MOI.OPTIMAL
            # stock heur sol 
            Y, X = stock_all_primal_sols(m, f1, f2, varArray)
            append!(Y_integer, Y) ; append!(X_integer, X)
            ys_1 = JuMP.value(f1) ; ys_2 = JuMP.value(f2)

            idx = push!(L.natural_order_vect, Solution(JuMP.value.(varArray), [ys_1, ys_2], curr_λ ), filtered=true)
            idx > 0 ? updateCT(L.natural_order_vect.sols[idx]) : nothing

            if isapprox(yr_1, ys_1, atol=1e-3) || isapprox(yr_2, ys_2, atol=1e-3)
                return Y_integer, X_integer
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

            idx = push!(L.natural_order_vect, Solution(x_star, [ys_1, ys_2], curr_λ ), filtered=true)
            idx > 0 ? updateCT(L.natural_order_vect.sols[idx]) : nothing

            if isapprox(yr_1, ys_1, atol=1e-3) || isapprox(yr_2, ys_2, atol=1e-3)
                return Y_integer, X_integer
            end
        else
            println("has primal ? $(JuMP.has_values(m))")
            error("Condition  status $status ")
        end
    end

    # -----------------------------------------
    # step 3 : fix the next search direction 
    # -----------------------------------------
    if length(L.natural_order_vect.sols) < 2 
        return Y_integer, X_integer
    end

    todo = []; 
    
    push!(todo, [L.natural_order_vect.sols[1].y, L.natural_order_vect.sols[end].y] )

    # ----------------------------------------------------------------
    # repeat the same procedure until no more direction in todo list 
    # ----------------------------------------------------------------
    # println("start with LBS : ", L.natural_order_vect)
    iter = 1
    while length(todo) > 0 # && iter ≤ 5
        # println("----------------------")
        println( "\t iter = $iter ..." )
        iter += 1
        # iter ≥ 15 ? @info "iter $ter " : nothing
        # println("----------------------")

        p = popfirst!(todo) ; yl = p[1] ;  yr = p[2]
        λ = [ abs(yr[2] - yl[2]) , abs(yl[1] - yr[1]) ] 

        # println("yl = $yl \t yr = $yr \t λ = $λ ")

        # solve the mono scalarization problem 
        f = AffExpr(0.0)    
        lb = λ[1] * yl[1] + λ[2] * yl[2]  
        JuMP.set_objective(m, f1Sense, λ[1]*f1 + λ[2]*f2) ; JuMP.set_objective(lp_copied, f1Sense, λ[1]*f1_copied + λ[2]*f2_copied)
        verbose && println("solving for $(λ[1])*f1 + $(λ[2])*f2")    
        f = λ[1]*f1 + λ[2]*f2

        x_star = [] ; bst_val = -Inf 
        curr_λ = λ
        JuMP.optimize!(m, ignore_optimize_hook=true) ; status = JuMP.termination_status(m)
    
        yt_1 = 0.0 ; yt_2 = 0.0 
        val = -Inf ; idx = -1

        if status == MOI.INFEASIBLE 
            # @info "scalarization infeasible ! "
            continue 
        end

        if status == MOI.OPTIMAL 
            # @info "scalarization optimal ! "

            # stock heur sol 
            Y, X = stock_all_primal_sols(m, f1, f2, varArray)
            append!(Y_integer, Y) ; append!(X_integer, X)
    
            yt_1 = JuMP.value(f1) ; yt_2 = JuMP.value(f2)
            val = λ[1]*yt_1 + λ[2]*yt_2 

            # same point 
            if (isapprox(yt_1, yl[1], atol=1e-3) && isapprox(yt_2, yl[2], atol=1e-3) ) || 
                (isapprox(yt_1, yr[1], atol=1e-3) && isapprox(yt_2, yr[2], atol=1e-3) )
                continue
            end

            # add new sol in LBS without filtering 
            idx = push!(L.natural_order_vect, Solution(JuMP.value.(varArray), [yt_1, yt_2], curr_λ ) )
            idx > 0 ? updateCT(L.natural_order_vect.sols[idx]) : nothing

            # println("idx = $idx , LBS : ", L.natural_order_vect)
    
        elseif status == MOI.NODE_LIMIT || status == TIME_LIMIT
            # @info "scalarization limit ! "

            if has_values(m)
                Y, X = stock_all_primal_sols(m, f1, f2, varArray)
                append!(Y_integer, Y) ; append!(X_integer, X)
            end
    
            if length(x_star) > 0
                nothing
            else
                best_bound = objective_bound(m)
                ctr_bound = JuMP.@constraint(lp_copied, λ[1]*f1_copied + λ[2]*f2_copied >= best_bound)
                JuMP.optimize!(lp_copied, ignore_optimize_hook=true)
    
                x_star = JuMP.value.(varArray_copied)
    
                if JuMP.is_valid(lp_copied, ctr_bound)
                    JuMP.delete(lp_copied, ctr_bound) ; JuMP.unregister(lp_copied, :ctr_bound)
                end
            end
    
            yt_1 = x_star'* c[1, 2:end] + c[1, 1]
            yt_2 = x_star'* c[2, 2:end] + c[2, 1]
            val = λ[1]*yt_1 + λ[2]*yt_2

            # same point
            if (isapprox(yt_1, yl[1], atol=1e-3) && isapprox(yt_2, yl[2], atol=1e-3) ) || 
                (isapprox(yt_1, yr[1], atol=1e-3) && isapprox(yt_2, yr[2], atol=1e-3) )
                continue
            end

            # add new sol in LBS without filtering 
            idx = push!(L.natural_order_vect, Solution(x_star, [yt_1, yt_2], curr_λ ) )
            idx > 0 ? updateCT(L.natural_order_vect.sols[idx]) : nothing

            # println("idx = $idx , LBS : ", L.natural_order_vect)

        else
            println("has primal ? $(JuMP.has_values(m))")
            error("Condition  status $status ")
        end


        # -----------------------------
        # case : equality    
        # -----------------------------
        if (abs(val - lb) ≤ 1e-4) || idx < 1
            continue
        end

        # find point intersection 
        intersection = intersectionPts(L, idx)
        # println("intersection pts: ", intersection )

        valid = true
        for s in L.natural_order_vect.sols
            if s.λ[1] * yt_1 + s.λ[2] * yt_2 < s.y[1] * s.λ[1] + s.y[2] * s.λ[2]-1e-4
                valid = false ; break
            end
        end

        # under the current LBS 
        if !valid 
            # @info "out of box !"
            deleteat!(L.natural_order_vect.sols, idx) ; 
            # add new intersection points 
            for s in intersection
                push!(L.natural_order_vect, s)
            end

            # println("update LBS : ", L.natural_order_vect)
            continue 
        end 

        bckp = L.natural_order_vect.sols[idx].y 
        # filter lower bounds under current line 
        filtering(val, L, λ)

        # add new intersection points 
        for s in intersection
            push!(L.natural_order_vect, s)
        end

        # println("update LBS : ", L.natural_order_vect)

        # define new search direction 
        idx = 0 ; located = false
        for s in L.natural_order_vect.sols
            idx +=1
            if s.y[1] == bckp[1] && s.y[2] == bckp[2] located = true ; break end
        end
        located ? next_direc(idx, L, todo) : nothing    
        
    end

    return Y_integer, X_integer
end


#todo : improve for re-optimization 
function opt_scalar_callbackalt(L::RelaxedBoundSet , m::JuMP.Model, lp_copied::JuMP.Model, c, λ, verbose; args...)
    global varArray
    global x_star
    global model = m
    global C = c
    global curr_λ = λ
    global bst_val

    vd = getvOptData(m)
    # empty!(vd.Y_N) ; empty!(vd.X_E) ; empty!(vd.lambda) 
    f1, f2 = vd.objs
    f1Sense, f2Sense = vd.objSenses
    varArray = JuMP.all_variables(m)

    varArray_copied = JuMP.all_variables(lp_copied)
    f1_copied = varArray_copied'* c[1, 2:end] + c[1, 1]
    f2_copied = varArray_copied'* c[2, 2:end] + c[2, 1]

    Y_integer = Vector{Vector{Float64}}() ; X_integer = Vector{Vector{Float64}}()
    
    # set up callback 
    MOI.set(m, MOI.NumberOfThreads(), 1) ; MOI.set(m, MOI.UserCutCallback(), callback_noCuts)

    JuMP.set_objective(m, f1Sense, λ[1]*f1 + λ[2]*f2) ; JuMP.set_objective(lp_copied, f1Sense, λ[1]*f1_copied + λ[2]*f2_copied)
    
    x_star = [] ; bst_val = -Inf 
    JuMP.optimize!(m, ignore_optimize_hook=true) ; status = JuMP.termination_status(m)

    yt_1 = 0.0 ; yt_2 = 0.0 
    val = -Inf ; idx = -1

    # in case of infeasibility 
    if status == MOI.INFEASIBLE 
        return Y_integer, X_integer
    end

    if status == MOI.OPTIMAL 
        # @info "scalarization optimal ! "

        # stock heur sol 
        Y, X = stock_all_primal_sols(m, f1, f2, varArray)
        append!(Y_integer, Y) ; append!(X_integer, X)

        yt_1 = JuMP.value(f1) ; yt_2 = JuMP.value(f2)
        val = λ[1]*yt_1 + λ[2]*yt_2 

        # add new sol in LBS without filtering 
        idx = push!(L.natural_order_vect, Solution(JuMP.value.(varArray), [yt_1, yt_2], curr_λ ) )
        updateCT(L.natural_order_vect.sols[idx])

        # println("idx = $idx , LBS : ", L.natural_order_vect)

    elseif status == MOI.NODE_LIMIT || status == TIME_LIMIT
        # @info "scalarization limit ! "

        if has_values(m)
            Y, X = stock_all_primal_sols(m, f1, f2, varArray)
            append!(Y_integer, Y) ; append!(X_integer, X)
        end

        if length(x_star) > 0
            nothing
        else
            best_bound = objective_bound(m)
            ctr_bound = JuMP.@constraint(lp_copied, λ[1]*f1_copied + λ[2]*f2_copied >= best_bound)
            JuMP.optimize!(lp_copied, ignore_optimize_hook=true)

            x_star = JuMP.value.(varArray_copied)

            if JuMP.is_valid(lp_copied, ctr_bound)
                JuMP.delete(lp_copied, ctr_bound) ; JuMP.unregister(lp_copied, :ctr_bound)
            end
        end

        yt_1 = x_star'* c[1, 2:end] + c[1, 1]
        yt_2 = x_star'* c[2, 2:end] + c[2, 1]
        val = λ[1]*yt_1 + λ[2]*yt_2

        # add new sol in LBS without filtering 
        idx = push!(L.natural_order_vect, Solution(x_star, [yt_1, yt_2], curr_λ ) )
        updateCT(L.natural_order_vect.sols[idx])

        # println("idx = $idx , LBS : ", L.natural_order_vect)

    else
        println("has primal ? $(JuMP.has_values(m))")
        error("Condition  status $status ")
    end


    # find point intersection 
    intersection = intersectionPts(L, idx)
    # println("intersection pts: ", intersection )

    valid = true
    for s in L.natural_order_vect.sols
        if s.λ[1] * yt_1 + s.λ[2] * yt_2 < s.y[1] * s.λ[1] + s.y[2] * s.λ[2]-1e-4
            valid = false ; break
        end
    end

    # under the current LBS 
    if !valid 
        # @info "out of box !"
        deleteat!(L.natural_order_vect.sols, idx) ; 
        # add new intersection points 
        for s in intersection
            push!(L.natural_order_vect, s)
        end

        # println("update LBS : ", L.natural_order_vect)
        return Y_integer, X_integer
    end 

    # filter lower bounds under current line 
    filtering(val, L, λ)

    # add new intersection points 
    for s in intersection
        push!(L.natural_order_vect, s)
    end

    # println("update LBS : ", L.natural_order_vect)

    return Y_integer, X_integer
end


