# This file contains functions of cutting planes algorithm.
include("BBtree.jl")
include("../algorithms.jl")
include("struct.jl")
include("separators.jl")
include("cutPool.jl")

using JuMP, CPLEX 

const max_step = 2
const loop_limit = 5


"""
Compute the lower bound set of the LP polyhedron by dichotomy method.

Return `true` if this node is fathomed by infeasibility.
"""
function compute_LBS(node::Node, pb::BO01Problem, round_results, verbose ; args...)
    #------------------------------------------------------------------------------
    # solve the LP relaxation by dichotomy method including the partial assignment
    #------------------------------------------------------------------------------
    solve_dicho(pb.m, round_results, false ; args...)
    vd_LP = getvOptData(pb.m)

    #-------------------------------------------------------------------------------
    # in case of the LP relaxed (sub) problem is infeasible, prune the actual node
    #-------------------------------------------------------------------------------
    if size(vd_LP.Y_N, 1) == 0
        prune!(node, INFEASIBILITY)
        if verbose
            @info "node $(node.num) is unfeasible !"
        end
        pb.info.nb_nodes_pruned += 1
        # pb.info.status = MOI.INFEASIBLE
        return true
    end

    # construct/complete the relaxed bound set
    node.RBS = RelaxedBoundSet()
    for i = 1:length(vd_LP.Y_N)
        push!(node.RBS.natural_order_vect, Solution(vd_LP.X_E[i], vd_LP.Y_N[i])) #TODO: filtered=true
    end
    # for i=1:length(node.RBS.natural_order_vect)-1
    #     node.RBS.segments[node.RBS.natural_order_vect.sols[i]] = true
    # end

    return false
end


"""
A single point cut off algorithm. 

Returns :
    - the fractional/integer point
    - true, if new cut(s) has benn found 
"""
function SP_cut_off(sol::Solution, node::Node, pb::BO01Problem, round_results, verbose ; args...)
    x_star = sol.xEquiv[1]
    if sol.is_binary return (x_star, false) end 

    # call generator
    # @info "Calling SP_CG_separator ..."
    # start_sep = time()
    # (isValidCut, α, viol) = SP_CG_separator(x_star, pb.A, pb.b)
    # pb.info.cuts_infos.times_calling_separators += (time() - start_sep)

    # if isValidCut
    #     @info " ------------------------- cut found"
    #     start_pool = time()
    #     push!(pb.cpool, α) && push_cutScore(node.cuts_ref, CutScore(size(pb.cpool.tab, 1), viol, 0)) # push cut reference
    #     pb.info.cuts_infos.times_oper_cutPool += (time() - start_pool)

    #     pb.info.cuts_infos.cuts_applied += 1 ; pb.info.cuts_infos.sp_cuts += 1
    #     con = JuMP.@constraint(pb.m, α[2:end]'*pb.varArray ≤ α[1]) ; push!(node.con_cuts, con)
    #     return (x_star, true)
    # end

    # @info "Calling SP_KP_heurSeparator ..."
    start_sep = time()
    cuts = SP_KP_heurSeparator2(x_star, pb.A, pb.b)
    pb.info.cuts_infos.times_calling_separators += (time() - start_sep)

    if length(cuts) > 0
        for cut in cuts
            viol = sum(cut[j+1]*(1-x_star[j]) for j = 1:length(cut)-1 )
            start_pool = time()
            ineq = Cut(cut)
            if push!(node.cutpool, ineq)
                pb.info.cuts_infos.cuts_applied += 1 ; pb.info.cuts_infos.sp_cuts += 1
                con = JuMP.@constraint(pb.m, cut[2:end]'*pb.varArray ≤ cut[1]) ; push!(node.con_cuts, con)
            end
            pb.info.cuts_infos.times_oper_cutPool += (time() - start_pool)

        end
        return (x_star, true)
    end

    return (x_star, false)
end


"""
Given two points `yₗ` and `yᵣ`, compute the extreme points between `yₗ` and `yᵣ`.

Return a list of extreme solutions in a natural order.
"""
function local_dichotomy(pb::BO01Problem, LBS::Vector{Solution}, idx_l::Int64, idx_r::Int64, round_results, verbose; args...)
    vd = getvOptData(pb.m)
    empty!(vd.Y_N) ; empty!(vd.X_E)
    f1, f2 = vd.objs
    f1Sense, f2Sense = vd.objSenses
    localSols = NaturalOrderVector()

    # store the beginning extreme pts 
    if idx_l == 0
        JuMP.set_objective(pb.m, f1Sense, f1) ; JuMP.optimize!(pb.m, ignore_optimize_hook=true)
        status = JuMP.termination_status(pb.m)
        if status == MOI.OPTIMAL
            yr_1 = JuMP.value(f1) ; yr_2 = JuMP.value(f2)
            push!(vd.Y_N, round_results ? round.([yr_1, yr_2]) : [yr_1, yr_2]) ; push!(vd.X_E, JuMP.value.(pb.varArray))
        else
            return localSols
        end
    else
        yr_1 = LBS[idx_l].y[1] ; yr_2 = LBS[idx_l].y[2]
    end

    if idx_r > length(LBS)
        JuMP.set_objective(pb.m, f2Sense, f2) ; JuMP.optimize!(pb.m, ignore_optimize_hook=true)
        status = JuMP.termination_status(pb.m)
        if status == MOI.OPTIMAL
            ys_1 = JuMP.value(f1) ; ys_2 = JuMP.value(f2)
            push!(vd.Y_N, round_results ? round.([ys_1, ys_2]) : [ys_1, ys_2]) ; push!(vd.X_E, JuMP.value.(pb.varArray))
        else
            return localSols
        end
    else
        ys_1 = LBS[idx_r].y[1] ; ys_2 = LBS[idx_r].y[2]
    end

    dichoRecursion(pb.m, yr_1, yr_2, ys_1, ys_2, pb.varArray, round_results, verbose ; args...)

    #Sort X_E and Y_N
    s = sortperm(vd.Y_N, by = first)
    vd.Y_N = vd.Y_N[s] ; vd.X_E = vd.X_E[s]
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
    deleteat!(vd.Y_N, inds) ; deleteat!(vd.X_E, inds)

    for i = 1:length(vd.Y_N)
        push!(localSols, Solution(vd.X_E[i], vd.Y_N[i])) # TODO : , filtered=true
    end
    empty!(vd.Y_N) ; empty!(vd.X_E)
    return localSols
end


function MP_cutting_planes(node::Node, pb::BO01Problem, round_results, verbose ; args...)  
    # println("----------------")
    # @info "node $(node.num)"
    # println("----------------")

    componentLBS = Vector{Vector{Solution}}() ; push!(componentLBS, node.RBS.natural_order_vect.sols)
    componentColored = Vector{Vector{Bool}}() ; node.RBS = RelaxedBoundSet()

    ite = 0
    while ite < loop_limit 
        ite += 1 ; pb.info.cuts_infos.ite_total += 1 
        ptsCounter = 0 ; ptsTotal = 0
        # ------------------------------------------------------------------------------
        # 1. generate multi-point cuts for each subset of LBS
        # ------------------------------------------------------------------------------
        for LBS in componentLBS
            colored = BO_cutting_planes(node, pb, LBS, round_results, verbose ; args...) ; ptsTotal += length(LBS)
            push!(componentColored, colored) ; ptsCounter += sum(colored)

            for i=1:length(LBS)
                if !colored[i]      # pt not cut off
                    push!(node.RBS.natural_order_vect, LBS[i]) # , filtered=true
                end
            end
        end

        # --------------------------------------------------
        # 2. stop if no more valid cut can be found
        # --------------------------------------------------
        if ptsCounter/ptsTotal < 0.4
            return false 
        end

        # ---------------------------------------------------
        # 3. otherwise, re-optimize by solving dicho -> LBS
        # ---------------------------------------------------
        start_dicho = time() ; len = size(componentLBS, 1)
        for _ in 1:len
            LBS = popfirst!(componentLBS) ; colored = popfirst!(componentColored)
            idx_l = colored[1] ? 0 : -1 ; idx_r = -1
            # println("début idx_l = $idx_l  ->  idx_r = $idx_r")
            # println("colored : ", colored)
            for i = 1:length(colored)
                if colored[i] continue end 
        
                if idx_l==-1 && (i+1)≤ length(colored) && colored[i+1]
                    idx_l = i ; continue
                end
        
                if idx_l != -1 && idx_r==-1
                    idx_r = i ; localSols = local_dichotomy(pb, LBS, idx_l, idx_r, round_results, verbose; args...)
                    if length(localSols) > 0 push!(componentLBS, localSols.sols) end 
                    # println("find interval idx_l = $idx_l  ->  idx_r = $idx_r")
                    # println("new sols : ", localSols)
                    if (i+1)≤ length(colored) && colored[i+1]
                        idx_l = i ; idx_r = -1
                    else
                        idx_l = -1; idx_r = -1
                    end
                end
            end
            if idx_l != -1 && idx_r == -1
                idx_r = length(colored)+1
                # println("find interval idx_l = $idx_l  ->  idx_r = $idx_r")
                localSols = local_dichotomy(pb, LBS, idx_l, idx_r, round_results, verbose; args...)
                if length(localSols) > 0 push!(componentLBS, localSols.sols) end 
                # println("new sols : ", localSols)
            end
        end
        
        pb.info.cuts_infos.times_calling_dicho += (time() - start_dicho)

        # in case of infeasibility #TODO : infeasible 
        if length(node.RBS.natural_order_vect)==0 && size(componentLBS, 1)==0 return true end

        # in case of integrity # TODO : integrity
        all_integers = true
        for LBS in componentLBS
            for sol in LBS
                if !sol.is_binary 
                    all_integers = false ; break 
                end
            end
        end
        for sol in node.RBS.natural_order_vect.sols
            if !sol.is_binary 
                all_integers = false ; break 
            end
        end
        if all_integers
            for LBS in componentLBS
                for sol in LBS  
                    push!(node.RBS.natural_order_vect, sol) # , filtered=true
                end
            end
            return false 
        end

    end # end loop while limit 
    for LBS in componentLBS
        for sol in LBS  
            push!(node.RBS.natural_order_vect, sol) # , filtered=true
        end
    end
    return false 
end




"""
Cutting planes scheme for the multi-point cuts. For now we assume that every point `y` in criteria space has only
one corresponding vector `x` in decision space. 

Return ture if the node is infeasible after adding cuts.
"""
function BO_cutting_planes(node::Node, pb::BO01Problem, LBS::Vector{Solution}, round_results, verbose ; args...)
    l = 1 #; cut_counter = 0
    colored = [false for _ in LBS]

    while l ≤ length(LBS)
        if LBS[l].is_binary 
            l += 1 ; continue    
        end

        for ∇ = max_step:-1:0 
            if ∇ == 0
                # @info "MP_cutting_planes calling `SP_cut_off`"
                (_, new_cut) = SP_cut_off(LBS[l], node, pb, round_results, verbose ; args...) 
                if new_cut colored[l] = true end # cut_counter += 1 
                l += 1
            else 
                r = l+∇
                if r > length(LBS) || LBS[r].is_binary continue end

                # @info "Calling MP_CG_separator ..."
                # start_sep = time()
                # (isValidCut, α, viol) = MP_CG_separator(LBS[l].xEquiv[1], LBS[r].xEquiv[1], pb.A, pb.b)
                # pb.info.cuts_infos.times_calling_separators += (time() - start_sep)

                # if isValidCut
                #     cut_counter += (∇+1)
                #     @info " ---------------------------- cut found "
                #     start_pool = time()
                #     push!(pb.cpool, α) && push_cutScore(node.cuts_ref, CutScore(size(pb.cpool.tab, 1), viol, 0)) # push cut reference
                #     pb.info.cuts_infos.times_oper_cutPool += (time() - start_pool)

                #     pb.info.cuts_infos.cuts_applied += 1 ; pb.info.cuts_infos.mp_cuts += 1
                #     con = JuMP.@constraint(pb.m, α[2:end]'*pb.varArray ≤ α[1]) ; push!(node.con_cuts, con)
                #     l = r + 1 ; break
                # end

                # @info "Calling MP_KP_heurSeparator ..."
                start_sep = time()
                cuts = MP_KP_heurSeparator2(LBS[l].xEquiv[1], LBS[r].xEquiv[1], pb.A, pb.b)
                pb.info.cuts_infos.times_calling_separators += (time() - start_sep)

                if length(cuts) > 0
                    for i =l:r colored[i] = true end 
                    # cut_counter += (∇+1)
                    for cut in cuts
                        viol_l = sum(cut[j+1]*(1-LBS[l].xEquiv[1][j]) for j = 1:length(cut)-1 )
                        viol_r = sum(cut[j+1]*(1-LBS[r].xEquiv[1][j]) for j = 1:length(cut)-1 )
                        viol = maximum([viol_l, viol_r])
                        
                        start_pool = time()
                        ineq = Cut(cut)
                        if push!(node.cutpool, ineq)
                            pb.info.cuts_infos.cuts_applied += 1 ; pb.info.cuts_infos.mp_cuts += 1
                            con = JuMP.@constraint(pb.m, cut[2:end]'*pb.varArray ≤ cut[1]) ; push!(node.con_cuts, con)
                        end
                        pb.info.cuts_infos.times_oper_cutPool += (time() - start_pool)

                    end
                    l = r + 1 ; break
                end

            end

        end
    end
    return colored
    # return cut_counter
end

