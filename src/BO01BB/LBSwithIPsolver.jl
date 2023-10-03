using JuMP
include("struct.jl")
TOL = 1e-4


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
            if isBinary(JuMP.value.(varArray, result = Int(i)))
                # push!(Y_integer, round.([JuMP.value(f1, result = Int(i)), JuMP.value(f2, result = Int(i))], digits = 4) )
                # push!(X_integer, round.(JuMP.value.(varArray, result = Int(i)), digits = 4) )
                push!(Y_integer, [JuMP.value(f1, result = Int(i)), JuMP.value(f2, result = Int(i))] )
                push!(X_integer, JuMP.value.(varArray, result = Int(i)) )
            else
                println("heuristic ", [JuMP.value(f1, result = Int(i)), JuMP.value(f2, result = Int(i))] )
                println(JuMP.value.(varArray, result = Int(i)))
            end

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
        if s.y'*λ ≤ lb-TOL push!(to_delete, i) end
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
        if length(L.natural_order_vect.sols[i].xEquiv[1]) > 0 l = i ; break end 
    end

    while j < length(L.natural_order_vect.sols)
        j += 1
        if length(L.natural_order_vect.sols[j].xEquiv[1]) > 0 r = j ; break end 
    end

    # prepare next directions # todo : !! verify deepcopy !! 
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
            t = L.natural_order_vect.sols[j]               
            if y'* t.λ ≤ t.y'* t.λ -TOL
                valid = false ; break
            end
        end

        if valid 
            new = Solution(Vector{Float64}(), y, [s1.λ[1], s1.λ[2]]) ; updateCT(new)
            push!(res, new )
        end
    end
    
    return res
end


function updateLBS(L::RelaxedBoundSet, idx::Int, val::Float64, curr_λ, yt)
    intersection = intersectionPts(L, idx)

    valid = true
    for s in L.natural_order_vect.sols        
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
straight line       a x + b y = c 
            <=>   -Δz2 z1 + Δz1 z2 = Δz1 c 

a = -Δz2 = -λ1   b = Δz1 = -λ2    ct = Δz1 c 
"""
function updateLBSwithEPB(node::Node)
    @assert length(node.nadirPt) > 1

    # 1er bounding 
    pt = Solution(Vector{Float64}(), node.nadirPt, [0.0, 1.0]) ; updateCT(pt)
    ptl = Solution() ; ptr = Solution()

    for i = 1:length(node.RBS.natural_order_vect.sols)
        s = node.RBS.natural_order_vect.sols[i]

        a1 = -s.λ[1] ; b1 = -s.λ[2]
        a2 = -pt.λ[1] ; b2 = -pt.λ[2]

        det = a1 * b2 - a2 * b1
        if abs(det) ≤ TOL continue end 
        
        y = [ (s.ct*b2 - pt.ct*b1)/det ,
              (a1* pt.ct - a2 * s.ct)/det
            ]
        valid = true
        for j = 1:length(node.RBS.natural_order_vect.sols)
            t = node.RBS.natural_order_vect.sols[j]               
            if y'* t.λ ≤ t.y'* t.λ -TOL
                valid = false ; break
            end
        end

        if valid 
            ptl = Solution(Vector{Float64}(), y, [1.0, 0.0]) ; updateCT(ptl)
            break
        end
    end

    # todo verify if bugs exist 
    if length(ptl.y) != 2 || ptl.y[1] == Inf || ptl.y[2] == Inf
        # println("\n --------------------- ")
        # println("ptl ", ptl)
        # println("nadri ", node.nadirPt, " \t boundz2 ", node.duplicationBound)
        # println("LBS ", node.RBS.natural_order_vect.sols)
        # error("EPB bounding error ! ")
        empty!(node.RBS.natural_order_vect.sols) ; return
    end
    
    # 2-th bounding 
    pt.λ = [1.0, 0.0] ; updateCT(pt)

    for i = 1:length(node.RBS.natural_order_vect.sols)
        s = node.RBS.natural_order_vect.sols[i]

        a1 = -s.λ[1] ; b1 = -s.λ[2]
        a2 = -pt.λ[1] ; b2 = -pt.λ[2]

        det = a1 * b2 - a2 * b1
        if abs(det) ≤ TOL continue end 
        
        y = [ (s.ct*b2 - pt.ct*b1)/det ,
              (a1* pt.ct - a2 * s.ct)/det
            ]
        valid = true
        for j = 1:length(node.RBS.natural_order_vect.sols)
            t = node.RBS.natural_order_vect.sols[j]               
            if y'* t.λ ≤ t.y'* t.λ -TOL
                valid = false ; break
            end
        end

        if valid 
            ptr = Solution(Vector{Float64}(), y, [0.0, 1.0]) ; updateCT(ptr)
            break
        end
    end

    # todo verify if bugs exist 
    if length(ptr.y) != 2 || ptr.y[1] == Inf || ptr.y[2] == Inf
        # println("\n --------------------- ")
        # println("ptr ", ptr )
        # println("nadri ", node.nadirPt, " \t boundz2 ", node.duplicationBound)
        # println("LBS ", node.RBS.natural_order_vect.sols)
        # error("EPB bounding error ! ")
        empty!(node.RBS.natural_order_vect.sols) ; return

    end

    # 3-th bounding 
    if node.duplicationBound != Inf && node.duplicationBound > ptr.y[2]
        pt = Solution(Vector{Float64}(), [node.nadirPt[1], node.duplicationBound], [0.0, 1.0]) ; updateCT(pt)
        for i = 1:length(node.RBS.natural_order_vect.sols)
            s = node.RBS.natural_order_vect.sols[i]
    
            a1 = -s.λ[1] ; b1 = -s.λ[2]
            a2 = -pt.λ[1] ; b2 = -pt.λ[2]
    
            det = a1 * b2 - a2 * b1
            if abs(det) ≤ TOL continue end 
            
            y = [ (s.ct*b2 - pt.ct*b1)/det ,
                  (a1* pt.ct - a2 * s.ct)/det
                ]
            valid = true
            for j = 1:length(node.RBS.natural_order_vect.sols)
                t = node.RBS.natural_order_vect.sols[j]               
                if y'* t.λ ≤ t.y'* t.λ -TOL
                    valid = false ; break
                end
            end
    
            if valid 
                ptr = Solution(Vector{Float64}(), y, [0.0, 1.0]) ; updateCT(ptr)
                break
            end
        end
    end

    # todo verify if bugs exist 
    if length(ptr.y) != 2 || ptr.y[1] == Inf || ptr.y[2] == Inf
        # println("\n --------------------- ")
        # println("ptr ", ptr )
        # println("nadri ", node.nadirPt, " \t boundz2 ", node.duplicationBound)
        # println("LBS ", node.RBS.natural_order_vect.sols)
        # error("EPB bounding error ! ")
        empty!(node.RBS.natural_order_vect.sols) ; return

    end

    # remove all points under current line 
    to_delete = Int64[] ; i = 1
    for s in node.RBS.natural_order_vect.sols
        if s.y'*ptl.λ ≤ ptl.y'*ptl.λ -TOL || s.y'*ptr.λ ≤ ptr.y'*ptr.λ -TOL ||
            s.y[2] ≥ ptl.y[2] + TOL || s.y[1] ≥ ptr.y[1] + TOL
            push!(to_delete, i) 
        end
        i += 1
    end

    deleteat!(node.RBS.natural_order_vect.sols, to_delete)
    push!(node.RBS.natural_order_vect, ptl); push!(node.RBS.natural_order_vect, ptr)
end


"""
New correction iterative algorithm for LBS 
    taking account into intersection with odd LBS
"""
function LBSinvokingIPsolver(pb::BO01Problem , L::RelaxedBoundSet , K::Int64; args...)
    global varArray
    global x_star
    global model = pb.m
    global C = pb.c
    global curr_λ
    global bst_val
    iter_count = 0

    vd = getvOptData(pb.m)
    f1, f2 = vd.objs
    f1Sense, f2Sense = vd.objSenses
    varArray = JuMP.all_variables(pb.m)

    varArray_copied = JuMP.all_variables(pb.lp_copied)
    f1_copied = varArray_copied'* pb.c[1, 2:end] + pb.c[1, 1]
    f2_copied = varArray_copied'* pb.c[2, 2:end] + pb.c[2, 1]

    Y_integer = Vector{Vector{Float64}}() ; X_integer = Vector{Vector{Float64}}()
    
    # set up callback 
    MOI.set(pb.m, MOI.NumberOfThreads(), 1) ; MOI.set(pb.m, MOI.UserCutCallback(), callback_noCuts)
    pureL = RelaxedBoundSet()

    # -------------------------------------------
    # step 1 : calculate the left extreme point
    # ------------------------------------------- 
    JuMP.set_objective(pb.m, f1Sense, f1) ; JuMP.set_objective(pb.lp_copied, f1Sense, f1_copied)
    
    x_star = [] ; bst_val = -Inf 
    idx = -1 ; val = -Inf
    curr_λ = [1.0, 0.0] ; newPt = false
    JuMP.optimize!(pb.m, ignore_optimize_hook=true) ; status = JuMP.termination_status(pb.m)
    iter_count += 1

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
        Y, X = stock_all_primal_sols(pb.m, f1, f2, varArray)
        append!(Y_integer, Y) ; append!(X_integer, X)

        x_round = JuMP.value.(varArray)
        yr_1 = x_round'*pb.c[1, 2:end] + pb.c[1, 1] ; yr_2 = x_round'*pb.c[2, 2:end] + pb.c[2, 1]
        val = curr_λ[1]*yr_1 + curr_λ[2]*yr_2

        ext_l = Solution(x_round, [yr_1, yr_2], [curr_λ[1], curr_λ[2]] ) ; updateCT(ext_l)
        roundSol(pb, ext_l)

        idx, newPt = push!(L.natural_order_vect, ext_l) 
        updateLBS(L, idx, val, [curr_λ[1], curr_λ[2]], [yr_1, yr_2]) 

        idx, newPt = push!(pureL.natural_order_vect, ext_l) 
        updateLBS(pureL, idx, val, [curr_λ[1], curr_λ[2]], [yr_1, yr_2])

    # otherwise, take the best primal sol so far 
    elseif status == MOI.NODE_LIMIT || status == TIME_LIMIT
        # stock heuristic sol 
        if has_values(pb.m)
            Y, X = stock_all_primal_sols(pb.m, f1, f2, varArray) ; append!(Y_integer, Y) ; append!(X_integer, X)
        end

        if length(x_star) > 0
            nothing
        else
            best_bound = objective_bound(pb.m)
            ctr_bound = JuMP.@constraint(pb.lp_copied, f1_copied >= best_bound)
            JuMP.optimize!(pb.lp_copied, ignore_optimize_hook=true)
            x_star = JuMP.value.(varArray_copied)

            if JuMP.is_valid(pb.lp_copied, ctr_bound)
                JuMP.delete(pb.lp_copied, ctr_bound) ; JuMP.unregister(pb.lp_copied, :ctr_bound)
            end
        end

        yr_1 = x_star'* pb.c[1, 2:end] + pb.c[1, 1]
        yr_2 = x_star'* pb.c[2, 2:end] + pb.c[2, 1]
        val = curr_λ[1]*yr_1 + curr_λ[2]*yr_2

        ext_l = Solution(x_star, [yr_1, yr_2], [curr_λ[1], curr_λ[2]] ) ; updateCT(ext_l)
        roundSol(pb, ext_l)

        idx, newPt = push!(L.natural_order_vect, ext_l)
        updateLBS(L, idx, val, [curr_λ[1], curr_λ[2]], [yr_1, yr_2])

        idx, newPt = push!(pureL.natural_order_vect, ext_l)
        updateLBS(pureL, idx, val, [curr_λ[1], curr_λ[2]], [yr_1, yr_2])

    else
        println("has primal ? $(JuMP.has_values(m))")
        error("Condition  status $status ")
    end

    if iter_count ≥ K  return Y_integer, X_integer end
    # -------------------------------------------
    # step 2 : calculate the right extreme point
    # ------------------------------------------- 
    JuMP.set_objective(pb.m, f2Sense, f2) ; JuMP.set_objective(pb.lp_copied, f2Sense, f2_copied)

    x_star = [] ; bst_val = -Inf 
    idx = -1 ; val = -Inf
    curr_λ = [0.0, 1.0] ; newPt = false
    JuMP.optimize!(pb.m, ignore_optimize_hook=true) ; status = JuMP.termination_status(pb.m)
    iter_count += 1

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
        Y, X = stock_all_primal_sols(pb.m, f1, f2, varArray)
        append!(Y_integer, Y) ; append!(X_integer, X)

        x_round = JuMP.value.(varArray)
        ys_1 = x_round'*pb.c[1, 2:end] + pb.c[1, 1]; ys_2 = x_round'*pb.c[2, 2:end] + pb.c[2, 1]
        val = curr_λ[1]*ys_1 + curr_λ[2]*ys_2

        ext_r = Solution(x_round, [ys_1, ys_2], [curr_λ[1], curr_λ[2]] ) ; updateCT(ext_r)
        roundSol(pb, ext_r)

        idx, newPt = push!(L.natural_order_vect, ext_r)
        updateLBS(L, idx, val, [curr_λ[1], curr_λ[2]], [ys_1, ys_2])

        idx, newPt = push!(pureL.natural_order_vect, ext_r) 
        updateLBS(pureL, idx, val, [curr_λ[1], curr_λ[2]], [ys_1, ys_2])

    elseif status == MOI.NODE_LIMIT || status == TIME_LIMIT
        if has_values(pb.m)
            Y, X = stock_all_primal_sols(pb.m, f1, f2, varArray)
            append!(Y_integer, Y) ; append!(X_integer, X)
        end

        if length(x_star) > 0
            nothing
        else
            best_bound = objective_bound(pb.m)
            ctr_bound = JuMP.@constraint(pb.lp_copied, f2_copied >= best_bound)
            JuMP.optimize!(pb.lp_copied, ignore_optimize_hook=true)

            x_star = JuMP.value.(varArray_copied)

            if JuMP.is_valid(pb.lp_copied, ctr_bound)
                JuMP.delete(pb.lp_copied, ctr_bound) ; JuMP.unregister(pb.lp_copied, :ctr_bound)
            end
        end

        ys_1 = x_star'* pb.c[1, 2:end] + pb.c[1, 1]
        ys_2 = x_star'* pb.c[2, 2:end] + pb.c[2, 1]
        val = curr_λ[1]*ys_1 + curr_λ[2]*ys_2

        ext_r = Solution(x_star, [ys_1, ys_2], [curr_λ[1], curr_λ[2]] ) ; updateCT(ext_r)
        roundSol(pb, ext_r)

        idx, newPt = push!(L.natural_order_vect, ext_r) 
        updateLBS(L, idx, val, [curr_λ[1], curr_λ[2]], [ys_1, ys_2]) 

        idx, newPt = push!(pureL.natural_order_vect, ext_r)
        updateLBS(pureL, idx, val, [curr_λ[1], curr_λ[2]], [ys_1, ys_2])

    else
        println("has primal ? $(JuMP.has_values(m))")
        error("Condition  status $status ")
    end

    if iter_count ≥ K return Y_integer, X_integer end

    # -----------------------------------------
    # step 3 : fix the next search direction 
    # -----------------------------------------
    if length(pureL.natural_order_vect.sols) < 2 
        return Y_integer, X_integer
    end

    l = 0 ; r = length(pureL.natural_order_vect.sols) +1
    while l < length(pureL.natural_order_vect.sols)
        l += 1
        if length(pureL.natural_order_vect.sols[l].xEquiv[1]) > 0 break end 
    end

    while r > 1
        r -= 1
        if length(pureL.natural_order_vect.sols[r].xEquiv[1]) > 0 break end 
    end

    if l ≥ r return Y_integer, X_integer end 

    todo = []; 
    
    push!(todo, [pureL.natural_order_vect.sols[l].y, pureL.natural_order_vect.sols[r].y] )

    # ----------------------------------------------------------------
    # repeat the same procedure until no more direction in todo list 
    # ----------------------------------------------------------------
    while length(todo) > 0
        if iter_count ≥ K return Y_integer, X_integer end

        p = popfirst!(todo) ; yl = p[1] ;  yr = p[2]

        λ = [ abs(yr[2] - yl[2]) , abs(yl[1] - yr[1]) ] 

        # solve the mono scalarization problem 
        f = AffExpr(0.0)    
        lb = λ'* yl
        JuMP.set_objective(pb.m, f1Sense, λ[1]*f1 + λ[2]*f2) ; JuMP.set_objective(pb.lp_copied, f1Sense, λ[1]*f1_copied + λ[2]*f2_copied)
        f = λ[1]*f1 + λ[2]*f2

        x_star = [] ; bst_val = -Inf 
        curr_λ = λ ; newPt = false 
        JuMP.optimize!(pb.m, ignore_optimize_hook=true) ; status = JuMP.termination_status(pb.m)
        iter_count += 1

        yt_1 = 0.0 ; yt_2 = 0.0 
        val = -Inf ; idx = -1
        idxL = -1 ; newPtL = false 
        pt = Solution()
        if status == MOI.INFEASIBLE continue end

        if status == MOI.OPTIMAL 
            # stock heur sol 
            Y, X = stock_all_primal_sols(pb.m, f1, f2, varArray)
            append!(Y_integer, Y) ; append!(X_integer, X)
    
            x_round = JuMP.value.(varArray)
            yt_1 = x_round'*pb.c[1, 2:end] + pb.c[1, 1] ; yt_2 = x_round'*pb.c[2, 2:end] + pb.c[2, 1]
            val = λ[1]*yt_1 + λ[2]*yt_2 
            if (isapprox(yt_1, yl[1], atol=TOL) && isapprox(yt_2, yl[2], atol=TOL) ) || 
                (isapprox(yr[1], yt_1, atol=TOL) && isapprox(yr[2], yt_2, atol=TOL) )
                continue
            end

            # add new sol in LBS without filtering 
            pt = Solution(x_round, [yt_1, yt_2], [curr_λ[1], curr_λ[2]] ) ; updateCT(pt)
            roundSol(pb, pt)

            idxL, newPtL = push!(L.natural_order_vect, pt)

            idx, newPt = push!(pureL.natural_order_vect, pt)

        elseif status == MOI.NODE_LIMIT || status == TIME_LIMIT
            if has_values(pb.m)
                Y, X = stock_all_primal_sols(pb.m, f1, f2, varArray)
                append!(Y_integer, Y) ; append!(X_integer, X)
            end
    
            if length(x_star) > 0
                nothing
            else
                best_bound = objective_bound(pb.m)
                ctr_bound = JuMP.@constraint(pb.lp_copied, λ[1]*f1_copied + λ[2]*f2_copied >= best_bound)
                JuMP.optimize!(pb.lp_copied, ignore_optimize_hook=true)
    
                x_star = JuMP.value.(varArray_copied)
    
                if JuMP.is_valid(pb.lp_copied, ctr_bound)
                    JuMP.delete(pb.lp_copied, ctr_bound) ; JuMP.unregister(pb.lp_copied, :ctr_bound)
                end
            end
    
            yt_1 = x_star'* pb.c[1, 2:end] + pb.c[1, 1]
            yt_2 = x_star'* pb.c[2, 2:end] + pb.c[2, 1]
            val = λ[1]*yt_1 + λ[2]*yt_2
            if (isapprox(yt_1, yl[1], atol=TOL) && isapprox(yt_2, yl[2], atol=TOL) ) || 
                (isapprox(yr[1], yt_1, atol=TOL) && isapprox(yr[2], yt_2, atol=TOL) )
                continue
            end

            # add new sol in LBS without filtering 
            pt = Solution(x_star, [yt_1, yt_2], [curr_λ[1], curr_λ[2]] ); updateCT(pt)
            roundSol(pb, pt)

            idxL, newPtL = push!(L.natural_order_vect, pt) 

            idx, newPt = push!(pureL.natural_order_vect, pt ) 

        else
            println("has primal ? $(JuMP.has_values(m))")
            error("Condition  status $status ")
        end
        
        updateLBS(L, idxL, val, [curr_λ[1], curr_λ[2]], [yt_1, yt_2])

        # -----------------------------
        # case : equality    
        # -----------------------------
        if !newPt continue end 

        # find point intersection 
        intersection = intersectionPts(pureL, idx,)

        valid = true
        for s in pureL.natural_order_vect.sols
            if s.λ[1] * yt_1 + s.λ[2] * yt_2 < s.y'* s.λ - TOL
                valid = false ; break
            end
        end

        # under the current LBS 
        if !valid 
            deleteat!(pureL.natural_order_vect.sols, idx)
            filtering(val, pureL, [λ[1], λ[2] ])

            # add new intersection points 
            for s in intersection
                push!(pureL.natural_order_vect, s)
            end
            continue 
        end 

        bckp = pureL.natural_order_vect.sols[idx].y 
        # filter lower bounds under current line 
        filtering(val, pureL, [λ[1], λ[2] ])

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
    end

    return Y_integer, X_integer

end


#todo : improve for re-optimization 
function opt_scalar_callbackalt(pb::BO01Problem, L::RelaxedBoundSet, λ; args...)
    global varArray
    global x_star
    global model = pb.m
    global C = pb.c
    global curr_λ = λ
    global bst_val

    vd = getvOptData(pb.m)
    f1, f2 = vd.objs
    f1Sense, f2Sense = vd.objSenses
    varArray = JuMP.all_variables(pb.m)

    varArray_copied = JuMP.all_variables(lp_copied)
    f1_copied = varArray_copied'* pb.c[1, 2:end] + pb.c[1, 1]
    f2_copied = varArray_copied'* pb.c[2, 2:end] + pb.c[2, 1]

    Y_integer = Vector{Vector{Float64}}() ; X_integer = Vector{Vector{Float64}}()
    
    # set up callback 
    MOI.set(pb.m, MOI.NumberOfThreads(), 1) ; MOI.set(pb.m, MOI.UserCutCallback(), callback_noCuts)

    JuMP.set_objective(pb.m, f1Sense, λ[1]*f1 + λ[2]*f2) ; JuMP.set_objective(pb.lp_copied, f1Sense, λ[1]*f1_copied + λ[2]*f2_copied)
    
    x_star = [] ; bst_val = -Inf 
    JuMP.optimize!(pb.m, ignore_optimize_hook=true) ; status = JuMP.termination_status(pb.m)

    yt_1 = 0.0 ; yt_2 = 0.0 
    val = -Inf ; idx = -1
    newPt = false ; pt = Solution()

    # in case of infeasibility 
    if status == MOI.INFEASIBLE 
        return Y_integer, X_integer
    end

    if status == MOI.OPTIMAL 
        # stock heur sol 
        Y, X = stock_all_primal_sols(pb.m, f1, f2, varArray)
        append!(Y_integer, Y) ; append!(X_integer, X)

        x_round = JuMP.value.(varArray)
        yt_1 = x_round'*pb.c[1, 2:end] + pb.c[1, 1] ; yt_2 = x_round'*pb.c[2, 2:end] + pb.c[2, 1]
        val = λ[1]*yt_1 + λ[2]*yt_2 

        # add new sol in LBS without filtering 
        pt = Solution(x_round, [yt_1, yt_2], [curr_λ[1], curr_λ[2]] ) ; updateCT(pt)
        roundSol(pb, pt)

        idx, newPt = push!(L.natural_order_vect, pt )

    elseif status == MOI.NODE_LIMIT || status == TIME_LIMIT
        if has_values(m)
            Y, X = stock_all_primal_sols(m, f1, f2, varArray)
            append!(Y_integer, Y) ; append!(X_integer, X)
        end

        if length(x_star) > 0
            nothing
        else
            best_bound = objective_bound(pb.m)
            ctr_bound = JuMP.@constraint(pb.lp_copied, λ[1]*f1_copied + λ[2]*f2_copied >= best_bound)
            JuMP.optimize!(pb.lp_copied, ignore_optimize_hook=true)

            x_star = JuMP.value.(varArray_copied)

            if JuMP.is_valid(pb.lp_copied, ctr_bound)
                JuMP.delete(pb.lp_copied, ctr_bound) ; JuMP.unregister(pb.lp_copied, :ctr_bound)
            end
        end

        yt_1 = x_star'* pb.c[1, 2:end] + pb.c[1, 1]
        yt_2 = x_star'* pb.c[2, 2:end] + pb.c[2, 1]
        val = λ[1]*yt_1 + λ[2]*yt_2

        # add new sol in LBS without filtering 
        pt = Solution(x_star, [yt_1, yt_2], [curr_λ[1], curr_λ[2]] ) ; updateCT(pt)
        roundSol(pb, pt)

        idx, newPt = push!(L.natural_order_vect,  pt)

    else
        println("has primal ? $(JuMP.has_values(m))")
        error("Condition  status $status ")
    end
    
    updateLBS(L, idx, val, [λ[1], λ[2] ], [yt_1, yt_2])

    return Y_integer, X_integer
end



function chordalImptovLBS(L::RelaxedBoundSet , m::JuMP.Model, lp_copied::JuMP.Model, c, K::Int64=1024; args...)
    step = max(0, length(L.natural_order_vect.sols) - K) 
    X = Vector{Vector{Float64}}() ; Y = Vector{Vector{Float64}}() 

    if K == 1024 || step == 0 return LBSinvokingIPsolveer(L, m, lp_copied, c, K=K; args...) end 

    idx_l = 1; idx_r = idx_l + step
    lambdas = []

    while idx_r ≤ length(L.natural_order_vect.sols) 

        λ = [abs(L.natural_order_vect.sols[idx_l].y[2] - L.natural_order_vect.sols[idx_r].y[2]), 
                abs(L.natural_order_vect.sols[idx_r].y[1] - L.natural_order_vect.sols[idx_l].y[1])
            ]
        push!(lambdas, λ)
        idx_l += 1 ; idx_r = idx_l + step
    end

    for λ in lambdas
        Y_integer, X_integer = opt_scalar_callbackalt(L, m, lp_copied, c, [λ[1], λ[2] ] ; args...)      
        append!(Y, Y_integer) ; append!(X, X_integer)
    end

    return X, Y
end



function dynamicImptovLBS(L::RelaxedBoundSet , m::JuMP.Model, lp_copied::JuMP.Model, c, K::Int64=1024; args...)
    X = Vector{Vector{Float64}}() ; Y = Vector{Vector{Float64}}() 

    if K == 1024 return LBSinvokingIPsolveer(L, m, lp_copied, c, K=K; args...) end 

    lambdas = []

    # {1/K, 2/K, ... K-1/K}
    for i=1:K
        w = i/(K+1)
        push!(lambdas, [w, 1-w])
    end

    for λ in lambdas
        Y_integer, X_integer = opt_scalar_callbackalt(L, m, lp_copied, c, [λ[1], λ[2] ] ; args...)      
        append!(Y, Y_integer) ; append!(X, X_integer)
    end

    return X, Y
end