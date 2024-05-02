# This file contains functions of cutting planes algorithm.
include("BBtree.jl")
include("../algorithms.jl")
include("separators.jl")
include("cutPool.jl")

using JuMP 
include("LBSwithIPsolver.jl")

const max_step = 2


"""
Compute the lower bound set of the LP polyhedron by dichotomy method.

Return `true` if this node is fathomed by infeasibility.
"""
function compute_LBS(node::Node, pb::BO01Problem, incumbent::IncumbentSet, round_results, verbose ; args...)
    # # todo : analyse to be deleted 
    # pb.info.λ_iter = 0.0 

    #------------------------------------------------------------------------------
    # solve the LP relaxation by dichotomy method including the partial assignment
    #------------------------------------------------------------------------------
    if pb.param.root_relax
        limits = 2^20
         
        # # todo (option) : EPB no LBS computation (bounding directly) seems not helpful 
        # if node.EPB && length(node.RBS.natural_order_vect.sols) > 0
        #     return false
        # end

        # todo (option) : λ limit tuning decreasing in depth (adaptive) not sure 
        if !pb.info.LBSexhaustive && !isRoot(node) && length(pb.varArray)- length(node.assignment) >10
            limits = pb.info.λ_limit
            # limits = round(Int, max(3, ceil(Int64 , pb.info.rootLBS / 5 ) ) )
        end

        # (option 1) : dichtomic-like concave-convex algorithm (default unlimited λ)
        if pb.info.λ_strategy == 0
            start = time()
            Y_integer, X_integer = LBSinvokingIPsolver(pb, node.RBS, limits; args...)      
            pb.info.relaxation_time += (time() - start)
            
        # (option 2) : chordal improvement 
        elseif pb.info.λ_strategy == 1
            start = time()
            Y_integer, X_integer = chordalImptovLBS(pb, node.RBS, limits; args...)
            pb.info.relaxation_time += (time() - start)

        # (option 3) : dynamic/equitable directions {1/K, 2/K, ... K-1/K}
        elseif pb.info.λ_strategy == 2
            start = time()
            Y_integer, X_integer = dynamicImptovLBS(pb, node.RBS, limits; args...)
            pb.info.relaxation_time += (time() - start)
        else
            @warn("λ searching strategy unknown with $(pb.info.λ_strategy)")
        end

        # # todo : analyse to be deleted 
        # pb.info.λ_acc += pb.info.λ_iter
        # if isRoot(node)
        #     pb.info.λ_max = pb.info.λ_iter ; pb.info.λ_min = pb.info.λ_iter 
        # else
        #     pb.info.λ_min = pb.info.λ_min < pb.info.λ_iter ? pb.info.λ_min : pb.info.λ_iter
        #     pb.info.λ_max = pb.info.λ_max > pb.info.λ_iter ? pb.info.λ_max : pb.info.λ_iter 
        # end

        start = time()
        for i = 1:length(Y_integer) 
            s = Solution(X_integer[i], Y_integer[i]) ; roundSol(pb, s)
            if s.is_binary push!(incumbent.natural_order_vect, s, filtered=true) end
        end
        pb.info.update_incumb_time += (time() - start) 

        if length(node.RBS.natural_order_vect.sols) == 0
            prune!(node, INFEASIBILITY)
            if verbose
                @info "node $(node.num) is unfeasible !"
            end
            return true
        end
    else
        # # todo (option) : EPB no LBS computation (bounding directly)
        # if node.EPB && length(node.RBS.natural_order_vect.sols) > 0
        #     return false
        # end
        start = time()
        solve_dicho(pb.m, round_results, false ; args...)
        pb.info.relaxation_time += (time() - start)

        vd_LP = getvOptData(pb.m)
    
        #-------------------------------------------------------------------------------
        # in case of the LP relaxed (sub) problem is infeasible, prune the actual node
        #-------------------------------------------------------------------------------
        if size(vd_LP.Y_N, 1) == 0
            prune!(node, INFEASIBILITY)
            if verbose
                @info "node $(node.num) is unfeasible !"
            end
            return true
        end
    
        # construct/complete the relaxed bound set
        node.RBS = RelaxedBoundSet()
        for i = 1:length(vd_LP.Y_N)
            s = Solution(vd_LP.X_E[i], vd_LP.Y_N[i], vd_LP.lambda[i]) ; roundSol(pb, s)
            push!(node.RBS.natural_order_vect, s, filtered=true )
        end
    end

    return false
end

# todo : partial reoptimize for root relaxation LBS only 
function reoptimize_LBS(node::Node, pb::BO01Problem, incumbent::IncumbentSet, cut_off, round_results, verbose ; args...)
    n = length(node.RBS.natural_order_vect.sols) ; lambdas = []

    for indx in cut_off
        push!(lambdas, [node.RBS.natural_order_vect.sols[indx].λ[1], node.RBS.natural_order_vect.sols[indx].λ[2] ] )
    end

    # delete all points cut off 
    # deleteat!(node.RBS.natural_order_vect.sols, cut_off)

    # in each direction, re-optimize 
    for λ in lambdas
        if pb.param.root_relax

            start = time() 
            Y_integer, X_integer = opt_scalar_callbackalt(pb, node.RBS, [λ[1], λ[2] ] ; args...)      
            pb.info.relaxation_time += (time() - start)
    
            start = time()
            for i = 1:length(Y_integer) 
                s = Solution(X_integer[i], Y_integer[i]) ; roundSol(pb, s)
                if s.is_binary push!(incumbent.natural_order_vect, s, filtered=true) end
            end
            pb.info.update_incumb_time += (time() - start)
        end
    end

end




"""
A single point cut off algorithm. 

Returns :
    - the fractional/integer point
    - true, if new cut(s) has benn found 
"""
function SP_cut_off(i::Int64, node::Node, pb::BO01Problem, round_results, verbose ; args...)
    x_star = node.RBS.natural_order_vect.sols[i].xEquiv[1]
    if node.RBS.natural_order_vect.sols[i].is_binary return (x_star, false) end 

    start_sep = time()
    cuts = SP_KP_heurSeparator2(x_star, pb.A, pb.b, node.assignment)
    pb.info.cuts_infos.times_calling_separators += (time() - start_sep)

    if length(cuts) > 0
        for cut in cuts
            start_pool = time()
            ineq = Cut(cut)
            if push!(node.cutpool, ineq)
                pb.info.cuts_infos.cuts_applied += 1 ; pb.info.cuts_infos.sp_cuts += 1
                con = JuMP.@constraint(pb.m, cut[2:end]'*pb.varArray ≤ cut[1]) ; push!(node.con_cuts, con)
                con = JuMP.@constraint(pb.lp_copied, cut[2:end]'*pb.varArray_copied ≤ cut[1]) ; push!(node.con_cuts_copied, con)
            end
            pb.info.cuts_infos.times_oper_cutPool += (time() - start_pool)

        end
        return (x_star, true)
    end

    # # call generator
    # start_sep = time()
    # (isValidCut, α, _) = SP_CG_separator(x_star, pb)
    # pb.info.cuts_infos.times_calling_separators += (time() - start_sep)

    # if isValidCut
    #     # @info "CG cut found ! "
    #     start_pool = time()
    #     ineq = Cut(α)
    #     if push!(node.cutpool, ineq)
    #         pb.info.cuts_infos.cuts_applied += 1 ; pb.info.cuts_infos.sp_cuts += 1
    #         con = JuMP.@constraint(pb.m, α[2:end]'*pb.varArray ≤ α[1]) ; push!(node.con_cuts, con)
    #     end
    #     pb.info.cuts_infos.times_oper_cutPool += (time() - start_pool)
    #     return (x_star, true)
    # end

    return (x_star, false)
end


function isCutable(node::Node, i::Int64, j::Int64, incumbent::IncumbentSet)::Bool
    for t = i+1:j-1 
        l = node.RBS.natural_order_vect.sols[t] 
        if l.is_binary
            return false 
        end
    end
    return true
end

"""
Cutting planes scheme for the multi-point cuts. For now we assume that every point `y` in criteria space has only
one corresponding vector `x` in decision space. 

Return ture if the node is infeasible after adding cuts.
"""
function MP_cutting_planes(node::Node, pb::BO01Problem, incumbent::IncumbentSet, loop_limit::Int64, round_results, verbose ; args...)
    numVars = length(pb.varArray) ; numRows = size(pb.A, 1)
    LBS = node.RBS.natural_order_vect.sols #; loop_limit = 5

    # # todo : 
    # println("-----------------------")
    # println("---- node $(node.num) ----")
    # println("-----------------------")


    # ------------------------------------------------------------------------------
    # 1. generate multi-point cuts if has any, or single-point cut off
    # ------------------------------------------------------------------------------
    ite = 0 ; if !isRoot(node) loop_limit = 1 end 
    while ite < loop_limit 
        ite += 1 ; pb.info.cuts_infos.ite_total += 1 

        l = 1 ; cut_counter = 0
        cut_off = []

        # # todo : 
        # TEST_MC = false
        # println(" iter $ite \n")

        while l ≤ length(LBS)
            if LBS[l].is_binary || length(LBS[l].xEquiv[1]) == 0 
                l += 1 ; continue    
            end
            
            for ∇ = max_step:-1:0 
                if ∇ == 0
                    (_, new_cut) = SP_cut_off(l, node, pb, round_results, verbose ; args...) 
                    if new_cut 
                        cut_counter += 1 ; push!(cut_off, l) 
                        # # todo : example 
                        # TEST_MC = true
                        # if TEST_MC @info "sp-cut (pts $l ) " end
                    end 
                    l += 1
                else 
                    r = l+∇
                    continue
                    
                    if r > length(LBS) || LBS[r].is_binary || length(LBS[r].xEquiv[1]) == 0 continue end

                    if ∇ > 1 && !isCutable(node, l, r, incumbent) continue end 

                    start_sep = time()
                    cuts = MP_KP_heurSeparator2(LBS[l].xEquiv[1], LBS[r].xEquiv[1], pb.A, pb.b, node.assignment)
                    pb.info.cuts_infos.times_calling_separators += (time() - start_sep)

                    if length(cuts) > 0
                        cut_counter += (∇+1) 
                        # # todo : example 
                        # TEST_MC = true
                        # if TEST_MC @info "multi-cut ($cut_counter pts $l to $r) " end

                        for i=l:r push!(cut_off, i) end 
                        for cut in cuts
                            start_pool = time()
                            ineq = Cut(cut)
                            if push!(node.cutpool, ineq)
                                pb.info.cuts_infos.cuts_applied += 1 ; pb.info.cuts_infos.mp_cuts += 1
                                con = JuMP.@constraint(pb.m, cut[2:end]'*pb.varArray ≤ cut[1]) ; push!(node.con_cuts, con)
                                con = JuMP.@constraint(pb.lp_copied, cut[2:end]'*pb.varArray_copied ≤ cut[1]) ; push!(node.con_cuts_copied, con)
                            end
                            pb.info.cuts_infos.times_oper_cutPool += (time() - start_pool)

                        end
                        l = r + 1 ; break
                    end

                    # start_sep = time()
                    # (isValidCut, α, _) = MP_CG_separator(LBS[l].xEquiv[1], LBS[r].xEquiv[1], pb)
                    # pb.info.cuts_infos.times_calling_separators += (time() - start_sep)

                    # if isValidCut
                    #     # @info "CG cut found ! "
                    #     cut_counter += (∇+1)
                    #     start_pool = time()
                    #     ineq = Cut(α)
                    #     if push!(node.cutpool, ineq)
                    #         pb.info.cuts_infos.cuts_applied += 1 ; pb.info.cuts_infos.mp_cuts += 1
                    #         con = JuMP.@constraint(pb.m, α[2:end]'*pb.varArray ≤ α[1]) ; push!(node.con_cuts, con)
                    #     end
                    #     pb.info.cuts_infos.times_oper_cutPool += (time() - start_pool)
                    #     l = r + 1 ; break
                    # end

                end

            end
        end

        # ---------------------------------------------------
        # 3. otherwise, re-optimize by solving dicho -> LBS
        # ---------------------------------------------------
        if cut_counter > 0
            # # todo : 
            # if TEST_MC
            #     println("--------------- before cut")
            #     println(LBS)
            # end

            if pb.param.root_relax
                reoptimize_LBS(node, pb, incumbent, cut_off, round_results, verbose; args)
            else
                pruned = compute_LBS(node, pb, incumbent, round_results, verbose; args)
                # # todo : 
                # if TEST_MC
                #     println("--------------- after cut")
                #     LBS = node.RBS.natural_order_vect.sols
                #     println(LBS)
                # end
                if pruned return true end
            end

            LBS = node.RBS.natural_order_vect.sols

            # in case of infeasibility
            if length(LBS)==0 return true end
        end

        # --------------------------------------------------
        # 2. stop if no more valid cut can be found
        # --------------------------------------------------
        if cut_counter/length(LBS) < 0.3
            return false 
        end
    end
    return false 
end
