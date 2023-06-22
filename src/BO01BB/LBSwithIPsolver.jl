using JuMP
include("struct.jl")
TOL = 1e-3


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
     
    if node_statut == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
    #     tmp = callback_value.(cb_data, varArray)
    #     val = (tmp'* C[1, 2:end] + C[1, 1]) * curr_λ[1] + (tmp'* C[2, 2:end] + C[2, 1] ) *curr_λ[2]
    #     if val > bst_val
    #         bst_val = val ; x_star = callback_value.(cb_data, varArray)
    #     end
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


"""
Supprime all lower bounds under current straight line.

    complexity : O(|L|) # todo : improve
"""
function filtering(lb::Float64, L::RelaxedBoundSet, λ)
    # remove all points under current line 
    to_delete = Int64[] ; i = 1
    for s in L.natural_order_vect.sols
        if s.y[1]*λ[1] + s.y[2]*λ[2] < lb - TOL push!(to_delete, i) end
        i += 1
    end

    deleteat!(L.natural_order_vect.sols, to_delete)
end


"""
Find the next search direction from current point idx.

    complexity : O(|L|)             #todo :  searching from extreme pts 
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
Computing intersected points between current point `idx` and LBS.
straight line       a x + b y = c 
            <=>   -Δz2 z1 + Δz1 z2 = Δz1 c 

a = -Δz2 = -λ1   b = Δz1 = -λ2    ct = Δz1 c 

line1 ->     a1 z11 + b1 z12 = ct1
             a1 = -λ11 b1 = -λ12

line2 ->     a2 z21 + b2 z22 = ct2
             a2 = -λ21 b2 = -λ22

 ----------- inter -------------
det = a1 b2 - a2 b1 (!!! ϵ près )
z1 = (ct1 b2 - ct2 b1)/det 
z2 = (a1 ct2 - a2 ct1)/det 

    complexity : O(|L|^2)  #todo : improve complexity 
"""
function intersectionPts(L::RelaxedBoundSet, idx::Int64)::Set{Solution}
    s2 = L.natural_order_vect.sols[idx] ; res = Set{Solution}()

    for i = 1:length(L.natural_order_vect.sols)
        if i== idx continue end 

        s1 = L.natural_order_vect.sols[i]

        a1 = -s1.λ[1] ; b1 = -s1.λ[2]
        a2 = -s2.λ[1] ; b2 = -s2.λ[2]

        det = a1 * b2 - a2 * b1
        if abs(det) ≤ TOL continue end 
        
        y = [ (s1.ct*b2 - s2.ct*b1)/det ,
              (a1* s2.ct - a2 * s1.ct)/det
            ]
        valid = true
        for j = 1:length(L.natural_order_vect.sols)
            t = L.natural_order_vect.sols[j]                # todo : check compared with segements
            if y[1] * t.λ[1] + y[2] * t.λ[2] < t.y[1] * t.λ[1] + t.y[2] * t.λ[2] - TOL
                valid = false ; break
            end
        end

        if valid 
            new = Solution(Vector{Vector{Float64}}(), y, false, s1.λ, Inf) ; updateCT(new)
            push!(res, new )
        end
    end
    
    return res
end

function updateLBS(L::RelaxedBoundSet, idx::Int, val::Float64, curr_λ, yt)
    intersection = intersectionPts(L, idx)

    valid = true
    for s in L.natural_order_vect.sols         # todo : compared with segments if the current pt is under LBS ?
        if s.λ'* yt < s.y'*s.λ -TOL
            valid = false ; break
        end
    end

    # under the current LBS 
    if !valid 
        deleteat!(L.natural_order_vect.sols, idx)
    end

    # filter lower bounds under current line 
    filtering(val, L, curr_λ)

    # add new intersection points 
    for s in intersection
        push!(L.natural_order_vect, s)
    end
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

    vd = getvOptData(m)
    f1, f2 = vd.objs
    f1Sense, f2Sense = vd.objSenses
    varArray = JuMP.all_variables(m)

    varArray_copied = JuMP.all_variables(lp_copied)
    f1_copied = varArray_copied'* c[1, 2:end] + c[1, 1]
    f2_copied = varArray_copied'* c[2, 2:end] + c[2, 1]

    Y_integer = Vector{Vector{Float64}}() ; X_integer = Vector{Vector{Float64}}()
    
    # set up callback 
    MOI.set(m, MOI.NumberOfThreads(), 1) ; MOI.set(m, MOI.UserCutCallback(), callback_noCuts)
    pureL = RelaxedBoundSet()

    # todo : 
    if verbose
        print("parent LBS = [ ")
        for s in L.natural_order_vect.sols
            print("$(s.y) , ")
        end
        println("] ")
        println()
        print("pureL = [ ")
        for s in pureL.natural_order_vect.sols
            print("$(s.y) , ")
        end
        println("] ")
    end

    # -------------------------------------------
    # step 1 : calculate the left extreme point
    # ------------------------------------------- 
    JuMP.set_objective(m, f1Sense, f1) ; JuMP.set_objective(lp_copied, f1Sense, f1_copied)
    
    x_star = [] ; bst_val = -Inf 
    idx = -1 ; val = -Inf
    curr_λ = [1.0, 0.0] ; newPt = false
    JuMP.optimize!(m, ignore_optimize_hook=true) ; status = JuMP.termination_status(m)

    # in case of infeasibility => !! no single extreme point 
    yr_1 = 0.0 ; yr_2 = 0.0 
    ext_l = Solution()
    if status == MOI.INFEASIBLE 
        empty!(L.natural_order_vect.sols)
        return Y_integer, X_integer
    end

    # in case of optimality 
    if status == MOI.OPTIMAL
        # stock heur sol 
        Y, X = stock_all_primal_sols(m, f1, f2, varArray)
        append!(Y_integer, Y) ; append!(X_integer, X)

        yr_1 = JuMP.value(f1) ; yr_2 = JuMP.value(f2)
        val = curr_λ[1]*yr_1 + curr_λ[2]*yr_2

        ext_l = Solution(JuMP.value.(varArray), [yr_1, yr_2], curr_λ ) ; updateCT(ext_l)
        idx, newPt = push!(L.natural_order_vect, ext_l) ; updateCT(L.natural_order_vect.sols[idx])
        updateLBS(L, idx, val, curr_λ, [yr_1, yr_2])

        idx, newPt = push!(pureL.natural_order_vect, ext_l) ; updateCT(pureL.natural_order_vect.sols[idx])
        updateLBS(pureL, idx, val, curr_λ, [yr_1, yr_2])

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
        val = curr_λ[1]*yr_1 + curr_λ[2]*yr_2

        ext_l = Solution(x_star, [yr_1, yr_2], curr_λ ) ; updateCT(ext_l)
        idx, newPt = push!(L.natural_order_vect, ext_l) ;  updateCT(L.natural_order_vect.sols[idx])
        updateLBS(L, idx, val, curr_λ, [yr_1, yr_2])

        idx, newPt = push!(pureL.natural_order_vect, ext_l) ; updateCT(pureL.natural_order_vect.sols[idx]) 
        updateLBS(pureL, idx, val, curr_λ, [yr_1, yr_2])

    else
        println("has primal ? $(JuMP.has_values(m))")
        error("Condition  status $status ")
    end

    # todo : 
    if verbose
        println("extl ", ext_l)
        print("parent LBS = [ ")
        for s in L.natural_order_vect.sols
            print("$(s.y) , ")
        end
        println("] ")
        println()
        print("pureL = [ ")
        for s in pureL.natural_order_vect.sols
            print("$(s.y) , ")
        end
        println("] ")
    end


    # -------------------------------------------
    # step 2 : calculate the right extreme point
    # ------------------------------------------- 
    JuMP.set_objective(m, f2Sense, f2) ; JuMP.set_objective(lp_copied, f2Sense, f2_copied)

    x_star = [] ; bst_val = -Inf 
    idx = -1 ; val = -Inf
    curr_λ = [0.0, 1.0] ; newPt = false
    JuMP.optimize!(m, ignore_optimize_hook=true) ; status = JuMP.termination_status(m)

    ys_1 = 0.0 ; ys_2 = 0.0 
    ext_r = Solution()
    if status == MOI.INFEASIBLE
        # todo : return single point OR the last updated LBS ??
        empty!(L.natural_order_vect.sols) ; push!(L.natural_order_vect, ext_l)
        # L.natural_order_vect.sols = deepcopy(pureL.natural_order_vect.sols)
        return Y_integer, X_integer
    end

    if status == MOI.OPTIMAL
        # stock heur sol 
        Y, X = stock_all_primal_sols(m, f1, f2, varArray)
        append!(Y_integer, Y) ; append!(X_integer, X)

        ys_1 = JuMP.value(f1) ; ys_2 = JuMP.value(f2)
        val = curr_λ[1]*ys_1 + curr_λ[2]*ys_2

        ext_r = Solution(JuMP.value.(varArray), [ys_1, ys_2], curr_λ ) ; updateCT(ext_r)
        idx, newPt = push!(L.natural_order_vect, ext_r) ; updateCT(L.natural_order_vect.sols[idx])
        updateLBS(L, idx, val, curr_λ, [ys_1, ys_2])

        idx, newPt = push!(pureL.natural_order_vect, ext_r) ; updateCT(pureL.natural_order_vect.sols[idx])
        updateLBS(pureL, idx, val, curr_λ, [ys_1, ys_2])

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
        val = curr_λ[1]*ys_1 + curr_λ[2]*ys_2

        ext_r = Solution(x_star, [ys_1, ys_2], curr_λ ) ; updateCT(ext_r)
        idx, newPt = push!(L.natural_order_vect, ext_r) ; updateCT(L.natural_order_vect.sols[idx])
        updateLBS(L, idx, val, curr_λ, [ys_1, ys_2])

        idx, newPt = push!(pureL.natural_order_vect, ext_r) ; updateCT(pureL.natural_order_vect.sols[idx])
        updateLBS(pureL, idx, val, curr_λ, [ys_1, ys_2])

    else
        println("has primal ? $(JuMP.has_values(m))")
        error("Condition  status $status ")
    end

    # todo : 
    if verbose
        println("extr ", ext_r)
        print("parent LBS = [ ")
        for s in L.natural_order_vect.sols
            print("$(s.y) , ")
        end
        println("] ")
        println()
        print("pureL = [ ")
        for s in pureL.natural_order_vect.sols
            print("$(s.y) , ")
        end
        println("] ")
    end

    # -----------------------------------------
    # step 3 : fix the next search direction 
    # -----------------------------------------
    if length(pureL.natural_order_vect.sols) < 2 
        return Y_integer, X_integer
    end

    l = 0 ; r = length(pureL.natural_order_vect.sols) +1
    while l < length(pureL.natural_order_vect.sols)
        l += 1
        if length(pureL.natural_order_vect.sols[l].xEquiv) > 0 break end 
    end

    while r > 1
        r -= 1
        if length(pureL.natural_order_vect.sols[r].xEquiv) > 0 break end 
    end

    if l ≥ r return Y_integer, X_integer end 

    todo = []; 
    
    push!(todo, [pureL.natural_order_vect.sols[l].y, pureL.natural_order_vect.sols[r].y] )

    # ----------------------------------------------------------------
    # repeat the same procedure until no more direction in todo list 
    # ----------------------------------------------------------------
    while length(todo) > 0

        p = popfirst!(todo) ; yl = p[1] ;  yr = p[2]
        # todo 
        # Δ2 = abs(yr[2] - yl[2])  ; Δ1 = abs(yl[1] - yr[1]) 
        # w = round(Δ2/(Δ2+Δ1), digits = 4)
        # λ = [w, round(1-w, digits=4)]      # normal to the segment
        λ = [ abs(yr[2] - yl[2]) , abs(yl[1] - yr[1]) ] 

        if verbose
            println("yl = ", yl , " yr = ", yr )
        end
        # solve the mono scalarization problem 
        f = AffExpr(0.0)    
        lb = λ[1] * yl[1] + λ[2] * yl[2]  
        JuMP.set_objective(m, f1Sense, λ[1]*f1 + λ[2]*f2) ; JuMP.set_objective(lp_copied, f1Sense, λ[1]*f1_copied + λ[2]*f2_copied)
        f = λ[1]*f1 + λ[2]*f2

        x_star = [] ; bst_val = -Inf 
        curr_λ = λ ; newPt = false 
        JuMP.optimize!(m, ignore_optimize_hook=true) ; status = JuMP.termination_status(m)
    
        yt_1 = 0.0 ; yt_2 = 0.0 
        val = -Inf ; idx = -1
        idxL = -1 ; newPtL = false 
        pt = Solution()
        if status == MOI.INFEASIBLE continue end

        if status == MOI.OPTIMAL 
            # stock heur sol 
            Y, X = stock_all_primal_sols(m, f1, f2, varArray)
            append!(Y_integer, Y) ; append!(X_integer, X)
    
            yt_1 = JuMP.value(f1) ; yt_2 = JuMP.value(f2)
            val = λ[1]*yt_1 + λ[2]*yt_2 

            # add new sol in LBS without filtering 
            pt = Solution(JuMP.value.(varArray), [yt_1, yt_2], curr_λ ) ; updateCT(pt)
            idxL, newPtL = push!(L.natural_order_vect, pt ) ; updateCT(L.natural_order_vect.sols[idxL])

            idx, newPt = push!(pureL.natural_order_vect, pt ) ; updateCT(pureL.natural_order_vect.sols[idx])

        elseif status == MOI.NODE_LIMIT || status == TIME_LIMIT
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
            pt = Solution(x_star, [yt_1, yt_2], curr_λ ); updateCT(pt)
            idxL, newPtL = push!(L.natural_order_vect, pt ) ; updateCT(L.natural_order_vect.sols[idxL])

            idx, newPt = push!(pureL.natural_order_vect, pt ) ; updateCT(pureL.natural_order_vect.sols[idx])

        else
            println("has primal ? $(JuMP.has_values(m))")
            error("Condition  status $status ")
        end

        # todo : 
        if verbose
            println("pt ", pt)
        end
        
        updateLBS(L, idxL, val, curr_λ, [yt_1, yt_2])

        # -----------------------------
        # case : equality    # todo : in case equality, filterage skipped 
        # -----------------------------
        if (abs(val - lb) ≤ TOL) ||!newPt continue end # 

        # todo : 
        if verbose
            print("parent LBS = [ ")
            for s in L.natural_order_vect.sols
                print("$(s.y) , ")
            end
            println("] ")
        end

        # find point intersection 
        intersection = intersectionPts(pureL, idx)

        valid = true
        for s in pureL.natural_order_vect.sols # todo : compared with segments if yt is under LBS 
            if s.λ[1] * yt_1 + s.λ[2] * yt_2 < s.y[1] * s.λ[1] + s.y[2] * s.λ[2]-TOL
                valid = false ; break
            end
        end

        # under the current LBS 
        if !valid 
            deleteat!(pureL.natural_order_vect.sols, idx) ; 
            filtering(val, pureL, λ)

            # add new intersection points 
            for s in intersection
                push!(pureL.natural_order_vect, s)
            end
            if verbose
                println()
                print("pureL = [ ")
                for s in pureL.natural_order_vect.sols
                    print("$(s.y) , ")
                end
                println("] ")
            end
            continue 
        end 

        bckp = pureL.natural_order_vect.sols[idx].y 
        # filter lower bounds under current line 
        filtering(val, pureL, λ)

        # add new intersection points 
        for s in intersection
            push!(pureL.natural_order_vect, s)
        end

        # define new search direction 
        idx = 0 ; located = false
        for s in pureL.natural_order_vect.sols
            idx +=1
            if s.y[1] == bckp[1] && s.y[2] == bckp[2] located = true ; break end
        end
        located ? next_direc(idx, pureL, todo) : nothing 
        if verbose
            println()
            print("pureL = [ ")
            for s in pureL.natural_order_vect.sols
                print("$(s.y) , ")
            end
            println("] ")
        end
    end

    # todo : 
    if verbose
        print("pureL = [ ")
        for s in pureL.natural_order_vect.sols
            print("$(s.y) , ")
        end
        println("] ")
        println()
        print("parent LBS = [ ")
        for s in L.natural_order_vect.sols
            print("$(s.y) , ")
        end
        println("] ")
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
    newPt = false ;

    # in case of infeasibility 
    if status == MOI.INFEASIBLE 
        return Y_integer, X_integer
    end

    if status == MOI.OPTIMAL 
        # stock heur sol 
        Y, X = stock_all_primal_sols(m, f1, f2, varArray)
        append!(Y_integer, Y) ; append!(X_integer, X)

        yt_1 = JuMP.value(f1) ; yt_2 = JuMP.value(f2)
        val = λ[1]*yt_1 + λ[2]*yt_2 

        # add new sol in LBS without filtering 
        idx, newPt = push!(L.natural_order_vect, Solution(JuMP.value.(varArray), [yt_1, yt_2], curr_λ ) )
        updateCT(L.natural_order_vect.sols[idx])

    elseif status == MOI.NODE_LIMIT || status == TIME_LIMIT
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
        idx, newPt = push!(L.natural_order_vect, Solution(x_star, [yt_1, yt_2], curr_λ ) )
        updateCT(L.natural_order_vect.sols[idx])

    else
        println("has primal ? $(JuMP.has_values(m))")
        error("Condition  status $status ")
    end
    
    updateLBS(L, idx, val, λ, [yt_1, yt_2])

    return Y_integer, X_integer
end




# --------------------------------------------
# explose here
# --------------------------------------------
# function LBSinvokingIPsolveer(L::RelaxedBoundSet , m::JuMP.Model, lp_copied::JuMP.Model, c, verbose; args...)
#     global varArray
#     global x_star
#     global model = m
#     global C = c
#     global curr_λ
#     global bst_val

#     vd = getvOptData(m)
#     f1, f2 = vd.objs
#     f1Sense, f2Sense = vd.objSenses
#     varArray = JuMP.all_variables(m)

#     varArray_copied = JuMP.all_variables(lp_copied)
#     f1_copied = varArray_copied'* c[1, 2:end] + c[1, 1]
#     f2_copied = varArray_copied'* c[2, 2:end] + c[2, 1]

#     Y_integer = Vector{Vector{Float64}}() ; X_integer = Vector{Vector{Float64}}()
    
#     # set up callback 
#     MOI.set(m, MOI.NumberOfThreads(), 1) ; MOI.set(m, MOI.UserCutCallback(), callback_noCuts)

#     # -------------------------------------------
#     # step 1 : calculate the left extreme point
#     # ------------------------------------------- 
#     JuMP.set_objective(m, f1Sense, f1) ; JuMP.set_objective(lp_copied, f1Sense, f1_copied)
#     verbose && println("solving for z1")
    
#     x_star = [] ; bst_val = -Inf 
#     idx = -1 ; val = -Inf
#     curr_λ = [1.0, 0.0] ; intersection = Set{Solution}()
#     JuMP.optimize!(m, ignore_optimize_hook=true) ; status = JuMP.termination_status(m)

#     # in case of infeasibility => !! no single extreme point 
#     yr_1 = 0.0 ; yr_2 = 0.0 
#     if status == MOI.INFEASIBLE 
#         empty!(L.natural_order_vect.sols)
#         return Y_integer, X_integer
#     end

#     # in case of optimality 
#     if status == MOI.OPTIMAL
#         # stock heur sol 
#         Y, X = stock_all_primal_sols(m, f1, f2, varArray)
#         append!(Y_integer, Y) ; append!(X_integer, X)

#         yr_1 = JuMP.value(f1) ; yr_2 = JuMP.value(f2)
#         val = curr_λ[1]*yr_1 + curr_λ[2]*yr_2

#         s = Solution(JuMP.value.(varArray), [yr_1, yr_2], curr_λ ) ; updateCT(s)
#         intersection = intersectionPts(L, s) ; idx = push!(L.natural_order_vect, s, filtered = true)
#         # idx > 0 ? updateCT(L.natural_order_vect.sols[idx]) : nothing

#     # otherwise, take the best primal sol so far 
#     elseif status == MOI.NODE_LIMIT || status == TIME_LIMIT
#         # stock heuristic sol 
#         if has_values(m)
#             Y, X = stock_all_primal_sols(m, f1, f2, varArray) ; append!(Y_integer, Y) ; append!(X_integer, X)
#         end

#         if length(x_star) > 0
#             nothing
#         else
#             best_bound = objective_bound(m)
#             ctr_bound = JuMP.@constraint(lp_copied, f1_copied >= best_bound)
#             JuMP.optimize!(lp_copied, ignore_optimize_hook=true)
#             x_star = JuMP.value.(varArray_copied)

#             if JuMP.is_valid(lp_copied, ctr_bound)
#                 JuMP.delete(lp_copied, ctr_bound) ; JuMP.unregister(lp_copied, :ctr_bound)
#             end
#         end

#         yr_1 = x_star'* c[1, 2:end] + c[1, 1]
#         yr_2 = x_star'* c[2, 2:end] + c[2, 1]
#         val = curr_λ[1]*yr_1 + curr_λ[2]*yr_2

#         s = Solution(x_star, [yr_1, yr_2], curr_λ ) ; updateCT(s)
#         intersection = intersectionPts(L, s) ; idx = push!(L.natural_order_vect, s, filtered = true)
#         # idx > 0 ? updateCT(L.natural_order_vect.sols[idx]) : nothing
#     else
#         println("has primal ? $(JuMP.has_values(m))")
#         error("Condition  status $status ")
#     end

#     # find intersections and filtering 
#     if idx > 0
#         # intersection = intersectionPts(L, idx)

#         # filter lower bounds under current line 
#         filtering(val, L, curr_λ)

#         # add new intersection points 
#         for s in intersection
#             push!(L.natural_order_vect, s, filtered = true)
#         end

#         valid = true
#         for s in L.natural_order_vect.sols         # todo : compared with segments if the current pt is under LBS ?
#             if s.λ[1] * yr_1 + s.λ[2] * yr_2 < s.y[1] * s.λ[1] + s.y[2] * s.λ[2]-1e-4
#                 valid = false ; break
#             end
#         end

#         # under the current LBS 
#         if !valid 
#             deleteat!(L.natural_order_vect.sols, idx) ; 
#         end
#     end

#     # -------------------------------------------
#     # step 2 : calculate the right extreme point
#     # ------------------------------------------- 
#     JuMP.set_objective(m, f2Sense, f2) ; JuMP.set_objective(lp_copied, f2Sense, f2_copied)
#     verbose && println("solving for z2")

#     x_star = [] ; bst_val = -Inf 
#     idx = -1 ; val = -Inf
#     curr_λ = [0.0, 1.0] ; intersection = Set{Solution}()
#     JuMP.optimize!(m, ignore_optimize_hook=true) ; status = JuMP.termination_status(m)

#     ys_1 = 0.0 ; ys_2 = 0.0 
#     if status == MOI.INFEASIBLE
#         # todo : return single point OR the last updated LBS ??
#         return Y_integer, X_integer
#     end

#     if status == MOI.OPTIMAL
#         # stock heur sol 
#         Y, X = stock_all_primal_sols(m, f1, f2, varArray)
#         append!(Y_integer, Y) ; append!(X_integer, X)

#         ys_1 = JuMP.value(f1) ; ys_2 = JuMP.value(f2)
#         val = curr_λ[1]*ys_1 + curr_λ[2]*ys_2

#         s = Solution(JuMP.value.(varArray), [ys_1, ys_2], curr_λ ) ; updateCT(s)
#         intersection = intersectionPts(L, s) ; idx = push!(L.natural_order_vect, s, filtered = true)
#         # idx > 0 ? updateCT(L.natural_order_vect.sols[idx]) : nothing

#     elseif status == MOI.NODE_LIMIT || status == TIME_LIMIT
#         if has_values(m)
#             Y, X = stock_all_primal_sols(m, f1, f2, varArray)
#             append!(Y_integer, Y) ; append!(X_integer, X)
#         end

#         if length(x_star) > 0
#             nothing
#         else
#             best_bound = objective_bound(m)
#             ctr_bound = JuMP.@constraint(lp_copied, f2_copied >= best_bound)
#             JuMP.optimize!(lp_copied, ignore_optimize_hook=true)

#             x_star = JuMP.value.(varArray_copied)

#             if JuMP.is_valid(lp_copied, ctr_bound)
#                 JuMP.delete(lp_copied, ctr_bound) ; JuMP.unregister(lp_copied, :ctr_bound)
#             end
#         end

#         ys_1 = x_star'* c[1, 2:end] + c[1, 1]
#         ys_2 = x_star'* c[2, 2:end] + c[2, 1]
#         val = curr_λ[1]*ys_1 + curr_λ[2]*ys_2

#         s = Solution(x_star, [ys_1, ys_2], curr_λ ) ; updateCT(s)
#         intersection = intersectionPts(L, s) ; idx = push!(L.natural_order_vect, s, filtered = true)
#         # idx > 0 ? updateCT(L.natural_order_vect.sols[idx]) : nothing
#     else
#         println("has primal ? $(JuMP.has_values(m))")
#         error("Condition  status $status ")
#     end

#     # find intersections and filtering 
#     if idx > 0
#         # intersection = intersectionPts(L, idx)

#         # filter lower bounds under current line 
#         filtering(val, L, curr_λ)

#         # add new intersection points 
#         for s in intersection
#             push!(L.natural_order_vect, s, filtered = true)
#         end

#         valid = true
#         for s in L.natural_order_vect.sols         # todo : compared with segments if the current pt is under LBS ?
#             if s.λ[1] * ys_1 + s.λ[2] * ys_2 < s.y[1] * s.λ[1] + s.y[2] * s.λ[2]-1e-4
#                 valid = false ; break
#             end
#         end

#         # under the current LBS 
#         if !valid 
#             deleteat!(L.natural_order_vect.sols, idx) ; 
#         end

#     end

#     # -----------------------------------------
#     # step 3 : fix the next search direction 
#     # -----------------------------------------
#     if length(L.natural_order_vect.sols) < 2 
#         return Y_integer, X_integer
#     end

#     l = 0 ; r = length(L.natural_order_vect.sols) +1
#     while l < length(L.natural_order_vect.sols)
#         l += 1
#         if length(L.natural_order_vect.sols[l].xEquiv) > 0 break end 
#     end

#     while r > 1
#         r -= 1
#         if length(L.natural_order_vect.sols[r].xEquiv) > 0 break end 
#     end

#     if l ≥ r return Y_integer, X_integer end 

#     todo = []; 
    
#     push!(todo, [L.natural_order_vect.sols[l].y, L.natural_order_vect.sols[r].y] )

#     # ----------------------------------------------------------------
#     # repeat the same procedure until no more direction in todo list 
#     # ----------------------------------------------------------------
#     iter = 1
#     while length(todo) > 0
#         # println("----------------------")
#         # @info "iter = $iter "
#         # iter += 1
#         # println("----------------------")

#         p = popfirst!(todo) ; yl = p[1] ;  yr = p[2]
#         λ = [ abs(yr[2] - yl[2]) , abs(yl[1] - yr[1]) ] 

#         # println("yl = $yl \t yr = $yr \t λ = $λ ")

#         # solve the mono scalarization problem 
#         f = AffExpr(0.0)    
#         lb = λ[1] * yl[1] + λ[2] * yl[2]  
#         JuMP.set_objective(m, f1Sense, λ[1]*f1 + λ[2]*f2) ; JuMP.set_objective(lp_copied, f1Sense, λ[1]*f1_copied + λ[2]*f2_copied)
#         verbose && println("solving for $(λ[1])*f1 + $(λ[2])*f2")    
#         f = λ[1]*f1 + λ[2]*f2

#         x_star = [] ; bst_val = -Inf 
#         curr_λ = λ ; intersection = Set{Solution}()
#         JuMP.optimize!(m, ignore_optimize_hook=true) ; status = JuMP.termination_status(m)
    
#         yt_1 = 0.0 ; yt_2 = 0.0 
#         val = -Inf ; idx = -1

#         if status == MOI.INFEASIBLE continue end

#         if status == MOI.OPTIMAL 
#             # stock heur sol 
#             Y, X = stock_all_primal_sols(m, f1, f2, varArray)
#             append!(Y_integer, Y) ; append!(X_integer, X)
    
#             yt_1 = JuMP.value(f1) ; yt_2 = JuMP.value(f2)
#             val = λ[1]*yt_1 + λ[2]*yt_2 

#             # add new sol in LBS without filtering 
#             s = Solution(JuMP.value.(varArray), [yt_1, yt_2], curr_λ ) ; updateCT(s)
#             intersection = intersectionPts(L, s) ; idx = push!(L.natural_order_vect, s, filtered = true)
#             # idx > 0 ? updateCT(L.natural_order_vect.sols[idx]) : nothing

#         elseif status == MOI.NODE_LIMIT || status == TIME_LIMIT
#             if has_values(m)
#                 Y, X = stock_all_primal_sols(m, f1, f2, varArray)
#                 append!(Y_integer, Y) ; append!(X_integer, X)
#             end
    
#             if length(x_star) > 0
#                 nothing
#             else
#                 best_bound = objective_bound(m)
#                 ctr_bound = JuMP.@constraint(lp_copied, λ[1]*f1_copied + λ[2]*f2_copied >= best_bound)
#                 JuMP.optimize!(lp_copied, ignore_optimize_hook=true)
    
#                 x_star = JuMP.value.(varArray_copied)
    
#                 if JuMP.is_valid(lp_copied, ctr_bound)
#                     JuMP.delete(lp_copied, ctr_bound) ; JuMP.unregister(lp_copied, :ctr_bound)
#                 end
#             end
    
#             yt_1 = x_star'* c[1, 2:end] + c[1, 1]
#             yt_2 = x_star'* c[2, 2:end] + c[2, 1]
#             val = λ[1]*yt_1 + λ[2]*yt_2

#             # add new sol in LBS without filtering 
#             s = Solution(x_star, [yt_1, yt_2], curr_λ ) ; updateCT(s)
#             intersection = intersectionPts(L, s) ; idx = push!(L.natural_order_vect, s, filtered = true )
#             # idx > 0 ? updateCT(L.natural_order_vect.sols[idx]) : nothing
#         else
#             println("has primal ? $(JuMP.has_values(m))")
#             error("Condition  status $status ")
#         end

#         # -----------------------------
#         # case : equality    # todo : in case equality, filterage skipped 
#         # -----------------------------
#         if (abs(val - lb) ≤ 1e-4) || idx < 0 continue end

#         # # find point intersection 
#         # intersection = intersectionPts(L, idx)

#         bckp = L.natural_order_vect.sols[idx].y 
#         # filter lower bounds under current line 
#         filtering(val, L, λ)

#         # add new intersection points 
#         for s in intersection
#             push!(L.natural_order_vect, s, filtered = true)
#         end

#         valid = true
#         for s in L.natural_order_vect.sols # todo : compared with segments if yt is under LBS 
#             if s.λ[1] * yt_1 + s.λ[2] * yt_2 < s.y[1] * s.λ[1] + s.y[2] * s.λ[2]-1e-4
#                 valid = false ; break
#             end
#         end

#         # under the current LBS 
#         if !valid 
#             deleteat!(L.natural_order_vect.sols, idx) ; 
#             continue 
#         end 


#         # define new search direction 
#         idx = 0 ; located = false
#         for s in L.natural_order_vect.sols
#             idx +=1
#             if s.y[1] == bckp[1] && s.y[2] == bckp[2] located = true ; break end
#         end
#         located ? next_direc(idx, L, todo) : nothing 
#     end

#     return Y_integer, X_integer
# end


# # todo : improve for re-optimization 
# function opt_scalar_callbackalt(L::RelaxedBoundSet , m::JuMP.Model, lp_copied::JuMP.Model, c, λ, verbose; args...)
#     global varArray
#     global x_star
#     global model = m
#     global C = c
#     global curr_λ = λ
#     global bst_val

#     vd = getvOptData(m)
#     f1, f2 = vd.objs
#     f1Sense, f2Sense = vd.objSenses
#     varArray = JuMP.all_variables(m)

#     varArray_copied = JuMP.all_variables(lp_copied)
#     f1_copied = varArray_copied'* c[1, 2:end] + c[1, 1]
#     f2_copied = varArray_copied'* c[2, 2:end] + c[2, 1]

#     Y_integer = Vector{Vector{Float64}}() ; X_integer = Vector{Vector{Float64}}()
    
#     # set up callback 
#     MOI.set(m, MOI.NumberOfThreads(), 1) ; MOI.set(m, MOI.UserCutCallback(), callback_noCuts)

#     JuMP.set_objective(m, f1Sense, λ[1]*f1 + λ[2]*f2) ; JuMP.set_objective(lp_copied, f1Sense, λ[1]*f1_copied + λ[2]*f2_copied)
    
#     x_star = [] ; bst_val = -Inf 
#     JuMP.optimize!(m, ignore_optimize_hook=true) ; status = JuMP.termination_status(m)

#     yt_1 = 0.0 ; yt_2 = 0.0 
#     val = -Inf ; idx = -1
#     intersection = Set{Solution}()

#     # in case of infeasibility 
#     if status == MOI.INFEASIBLE 
#         return Y_integer, X_integer
#     end

#     if status == MOI.OPTIMAL 
#         # stock heur sol 
#         Y, X = stock_all_primal_sols(m, f1, f2, varArray)
#         append!(Y_integer, Y) ; append!(X_integer, X)

#         yt_1 = JuMP.value(f1) ; yt_2 = JuMP.value(f2)
#         val = λ[1]*yt_1 + λ[2]*yt_2 

#         # add new sol in LBS without filtering 
#         s = Solution(JuMP.value.(varArray), [yt_1, yt_2], curr_λ ) ; updateCT(s)
#         intersection = intersectionPts(L, s) ; idx = push!(L.natural_order_vect, s, filtered = true )
#         # idx > 0 ? updateCT(L.natural_order_vect.sols[idx]) : nothing

#     elseif status == MOI.NODE_LIMIT || status == TIME_LIMIT
#         if has_values(m)
#             Y, X = stock_all_primal_sols(m, f1, f2, varArray)
#             append!(Y_integer, Y) ; append!(X_integer, X)
#         end

#         if length(x_star) > 0
#             nothing
#         else
#             best_bound = objective_bound(m)
#             ctr_bound = JuMP.@constraint(lp_copied, λ[1]*f1_copied + λ[2]*f2_copied >= best_bound)
#             JuMP.optimize!(lp_copied, ignore_optimize_hook=true)

#             x_star = JuMP.value.(varArray_copied)

#             if JuMP.is_valid(lp_copied, ctr_bound)
#                 JuMP.delete(lp_copied, ctr_bound) ; JuMP.unregister(lp_copied, :ctr_bound)
#             end
#         end

#         yt_1 = x_star'* c[1, 2:end] + c[1, 1]
#         yt_2 = x_star'* c[2, 2:end] + c[2, 1]
#         val = λ[1]*yt_1 + λ[2]*yt_2

#         # add new sol in LBS without filtering 
#         s = Solution(x_star, [yt_1, yt_2], curr_λ ) ; updateCT(s)
#         intersection = intersectionPts(L, s) ; idx = push!(L.natural_order_vect, s , filtered = true)
#         # idx > 0 ? updateCT(L.natural_order_vect.sols[idx]) : nothing
#     else
#         println("has primal ? $(JuMP.has_values(m))")
#         error("Condition  status $status ")
#     end

#     if idx > 0 
#         # find point intersection 
#         # intersection = intersectionPts(L, idx)

#         # filter lower bounds under current line 
#         filtering(val, L, λ)

#         # add new intersection points 
#         for s in intersection
#             push!(L.natural_order_vect, s, filtered = true)
#         end

#         valid = true
#         for s in L.natural_order_vect.sols # todo :  compared with segments if yt is under LBS 
#             if s.λ[1] * yt_1 + s.λ[2] * yt_2 < s.y[1] * s.λ[1] + s.y[2] * s.λ[2]-1e-4
#                 valid = false ; break
#             end
#         end

#         # under the current LBS 
#         if !valid 
#             deleteat!(L.natural_order_vect.sols, idx) ; 
#         end 
        
#     end

#     return Y_integer, X_integer
# end