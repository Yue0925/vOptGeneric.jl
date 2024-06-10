## This file contains functions related to node fathoming.

include("cuttingPlanes.jl")
include("GM.jl")

TOL = 1e-4
"""
Given a point `x_star`, iterate all valid cuts of parent node and stock the references 
in the current node.
"""
function loadingCutInPool(node::Node, pb::BO01Problem)
    if isRoot(node) return end 

    l = 1 ; LBS = node.RBS.natural_order_vect.sols

    while l ≤ length(LBS)
        if LBS[l].is_binary || length(LBS[l].xEquiv[1]) == 0 
            l += 1 ; continue
        end

        xₗ_star = LBS[l].xEquiv[1]

        for ∇ = max_step:-1:0 
            if ∇ == 0
                # single-point cut 
                for (k, cuts) in node.pred.cutpool.hashMap
                    for cut in cuts
                        α = cut.row
                        violationₗ = maximum([ (xₗ_star'*α[2:end] - α[1]), 0.0 ])
                        if violationₗ > 0.0
                            # ineq = Cut(α)
                            pb.info.cuts_infos.cuts_applied += 1 ; pb.info.cuts_infos.sp_cuts += 1
                            if push!(node.cutpool, cut)# && push_cutScore(node.cuts_ref, CutScore(length(node.cutpool.hashMap[k]), violationₗ, k))
                                con = JuMP.@constraint(pb.m, α[2:end]'*pb.varArray ≤ α[1]) ; push!(node.con_cuts, con)
                                con = JuMP.@constraint(pb.lp_copied, α[2:end]'*pb.varArray_copied ≤ α[1]) ; push!(node.con_cuts_copied, con)
                            end
                        end
                    end
                end
                l += 1
            else 
                applied = false
                r = l+∇
                if r > length(LBS) || LBS[r].is_binary || length(LBS[r].xEquiv[1]) == 0 continue end

                xᵣ_star = LBS[r].xEquiv[1]
                # multi-point cut 
                for (k, cuts) in node.pred.cutpool.hashMap
                    for cut in cuts 
                        α = cut.row
                        violationₗ = maximum([ (xₗ_star'*α[2:end] - α[1]), 0.0 ])
                        violationᵣ = maximum([ (xᵣ_star'*α[2:end] - α[1]), 0.0 ])
                        viol = maximum([violationₗ, violationᵣ])
                        if viol > 0.0
                            applied = true
                            # ineq = Cut(α)
                            pb.info.cuts_infos.cuts_applied += 1 ; pb.info.cuts_infos.mp_cuts += 1
                            if push!(node.cutpool, cut)# && push_cutScore(node.cuts_ref, CutScore(length(node.cutpool.hashMap[k]), viol, k))
                                con = JuMP.@constraint(pb.m, α[2:end]'*pb.varArray ≤ α[1]) ; push!(node.con_cuts, con)
                                con = JuMP.@constraint(pb.lp_copied, α[2:end]'*pb.varArray_copied ≤ α[1]) ; push!(node.con_cuts_copied, con)
                            end
                        end
                    end
                end
                if applied 
                    l = r+1 ; break
                end 
            end
        end
    end

end

function GM_heuristic(problem::BO01Problem, incumbent::IncumbentSet)
    nb_feas = 0 ; null_pt = 0
    GMtime = time()
    vg, nbgen = GM(problem.lp_copied, problem.varArray_copied, problem.c, 20, 30, 20)
    GMtime = time() - GMtime

    println(" GMtime = $GMtime \n total try = $nbgen ")
    problem.info.heur_time += GMtime
    # ----------------------------------------------------------

    for k = 1:nbgen 
        if vg[k].sFea
            nb_feas += 1
            (_, flag) = push!(incumbent.natural_order_vect, Solution(vg[k].sInt.x .*1.0, vg[k].sInt.y .* 1.0), filtered=true)
            flag ? nothing : null_pt += 1
        end
    end

    println("new feas = $nb_feas \t null_pt = $null_pt")

end

"""
Compute and stock the relaxed bound set (i.e. the LP relaxation) of the (sub)problem defined by the given node.
Return `true` if the node is pruned by infeasibility.
"""
function LPRelaxByDicho(node::Node, pb::BO01Problem, incumbent::IncumbentSet, round_results, verbose ; args...)::Bool
    objcons, objcons_copied = setVarObjBounds(node, pb)
    num_var = length(pb.varArray)
    start = time()

    # step 1 : calculate LBS of the actual sub-problem
    pruned = compute_LBS(node, pb, incumbent, round_results, verbose; args)
    
    # ------------------------
    # apply valid cuts 
    # ------------------------
    if pb.param.cp_activated && !pruned #&& node.depth < num_var/3
        @assert length(node.RBS.natural_order_vect) > 0 "valid LBS is empty !"

        loop_limit = 5
        # step 2 : add valid cuts constraints then re-optimize 
        start_processing = time()
        loadingCutInPool( node, pb)      # complexity O(pt ⋅ cuts)
        pb.info.cuts_infos.times_add_retrieve_cuts += (time() - start_processing)

        # todo : test re-opt
        if !pb.param.EPB && length(node.cutpool.hashMap) > 0
            pruned = compute_LBS(node, pb, incumbent, round_results, verbose; args)
            if pruned return true end
            loop_limit -= 1
        end
        # ---------------------------

        pruned = MP_cutting_planes(node, pb, incumbent, loop_limit, round_results, verbose ; args...)

        # # ----------------------------------------------------------
        # # todo : heuristics Gravity machine
        # if node.depth %10 == 0 && length(pb.varArray)- length(node.assignment) >10
        #     println("node $(node.depth)")
        #     GM_heuristic(pb, incumbent)
        # end

        # step 3 : retrieve applied valid cuts 
        start_processing = time()
        for con in node.con_cuts
            if JuMP.is_valid( pb.m, con)
                JuMP.delete( pb.m, con) ; JuMP.unregister( pb.m, :con) # remove the symbolic reference
            end
        end
        for con in node.con_cuts_copied
            if JuMP.is_valid( pb.lp_copied, con)
                JuMP.delete( pb.lp_copied, con) ; JuMP.unregister( pb.lp_copied, :con) # remove the symbolic reference
            end
        end
        pb.info.cuts_infos.times_add_retrieve_cuts += (time() - start_processing)

        pb.info.cuts_infos.times_total_for_cuts += (time() - start)        
    end

        # # ----------------------------------------------------------
        # # todo : heuristics Gravity machine
        # if node.depth %10 == 0 && length(pb.varArray)- length(node.assignment) >10
        #     println("node $(node.depth)")
        #     GM_heuristic(pb, incumbent)
        # end

    removeVarObjBounds(node, pb, objcons, objcons_copied) ; return pruned
end



"""
At the given node, update (filtered by dominance) the global incumbent set.
Return `true` if the node is pruned by optimality.
"""
function updateIncumbent(node::Node, pb::BO01Problem, incumbent::IncumbentSet, verbose)::Bool
    start = time()
    #-----------------------------------------------------------
    # check optimality && update the incumbent set
    #-----------------------------------------------------------

    for i = 1:length(node.RBS.natural_order_vect)
        if node.RBS.natural_order_vect.sols[i].is_binary
            s = node.RBS.natural_order_vect.sols[i]
            push!(incumbent.natural_order_vect, s, filtered=true) 
        end
    end

    if length(node.RBS.natural_order_vect)==1 && node.RBS.natural_order_vect.sols[1].is_binary
        prune!(node, OPTIMALITY)
        if verbose
            @info "node $(node.num) is fathomed by optimality ! and length = $(length(node.RBS.natural_order_vect))"
        end
        pb.info.update_incumb_time += (time() - start) ; return true
    end
    pb.info.update_incumb_time += (time() - start) ; return false
end

"""
Return local nadir points (so-called corner points) of the given incumbent set, or the single point it contains.
"""
function getNadirPoints(incumbent::IncumbentSet) # , ptl, ptr
    nadir_pts = NaturalOrderVector()

    if length(incumbent.natural_order_vect) == 1 return incumbent.natural_order_vect end 

    for i = 1:length(incumbent.natural_order_vect)-1

        push!(nadir_pts, Solution(
            Vector{Float64}(),
            [incumbent.natural_order_vect.sols[i+1].y[1],
            incumbent.natural_order_vect.sols[i].y[2]
            ])#, filtered=true
        )
    end

    return nadir_pts
end


"""
A fully explicit dominance test, and prune the given node if it's fathomed by dominance.
(i.e. ∀ l∈L: ∃ u∈U s.t. λu ≤ λl )
Return `true` if the given node is fathomed by dominance.
"""
function fullyExplicitDominanceTest(node::Node, incumbent::IncumbentSet, worst_nadir_pt::Vector{Float64}, EPB::Bool)
    @assert length(node.RBS.natural_order_vect) > 0 "relaxed bound set is empty for node $(node.num)"

    # we can't compare the LBS and UBS if the incumbent set is empty
    if length(incumbent.natural_order_vect) == 0 return false end

    nadir_pts = getNadirPoints(incumbent)

    # if there exists an upper bound u s.t. u≦l
    function weak_dom(l)
        for u ∈ nadir_pts.sols
            if u ≤ l && u != l
                return true
            end
        end
        return false
    end

    # ------------------------------------------
    # if the LBS consists of a single point
    # ------------------------------------------
    if length(node.RBS.natural_order_vect) == 1
        return weak_dom(node.RBS.natural_order_vect.sols[1])
    end

    # ----------------------------------------------
    # if the LBS consists of segments
    # ----------------------------------------------
    # two extreme points of LBS
    ptl = node.RBS.natural_order_vect.sols[1] ; ptr = node.RBS.natural_order_vect.sols[end]

    # Case 1 :  if only one feasible point in UBS 
    if length(incumbent.natural_order_vect) == 1 
        u = incumbent.natural_order_vect.sols[1]
        return u.y[1] < ptl.y[1] && u.y[2] < ptr.y[2]
    end

    # test range condition necessary 1 : LBS ⊆ UBS 
    # i.e. UBS includes the LP lexico-optimum
    u_l = incumbent.natural_order_vect.sols[1] ; u_r = incumbent.natural_order_vect.sols[end]

    sufficient = (u_l.y[1] < ptl.y[1] && u_r.y[2] < ptr.y[2])

    if !sufficient return false end

    # test condition necessary 2 : LBS ≤/dominates UBS 
    fathomed = true# ; dist_naditPt = Vector{Float64}()

    # iterate of all local nadir points
    for u ∈ nadir_pts.sols
        existence = false ; compared = false

        # case 1 : if u is dominates the ideal point of LBS 
        if u.y[2] < ptr.y[2] && u.y[1] < ptl.y[1]
            return true
        end

        # case 3 : complete pairwise comparison
        for i=1:length(node.RBS.natural_order_vect)-1              # ∀ segment l ∈ LBS 

            sol_l = node.RBS.natural_order_vect.sols[i] ; sol_r = node.RBS.natural_order_vect.sols[i+1]

            # -------------------------------------- spacial case LBS 
            if sol_l.y[2] == sol_r.y[2] # horizon line 
                if u.y[1] < sol_l.y[1] || u.y[1] > sol_r.y[1]
                    continue
                end

                compared = true 
                if u.y[2] < sol_l.y[2]
                    existence = true ; break
                end
            end

            if sol_l.y[1] == sol_r.y[1] # vertical line 
                if u.y[2] < sol_r.y[2] || u.y[2] > sol_l.y[2]
                    continue
                end

                compared = true 
                if u.y[1] < sol_l.y[1]
                    existence = true ; break
                end
            end
            # -------------------------------------

            if (u.y[2] > sol_l.y[2] || u.y[2] < sol_r.y[2]) && (u.y[1] > sol_r.y[1] || u.y[1] < sol_l.y[1])
                continue
            end

            λ = [sol_l.y[2] - sol_r.y[2], sol_r.y[1] - sol_l.y[1]]      # normal to the segment

            compared = true

            if λ'*u.y < λ'*sol_r.y -TOL 
                existence = true ; break
            end
        end
        
        # case 4 : condition dominance violated, then stock the non-dominated local nadir pts to prepare EPB
        if compared && !existence 
            fathomed = false
            if EPB
                # if !isRoot(node) && isNadirPointDuplicated(node, u.y)   # the current local nadir pt is already branched 
                if !isRoot(node) && (u.y in node.pred.localNadirPts || u.y == node.pred.nadirPt || u.y == node.nadirPt)    # the current local nadir pt is already branched 
                    node.localNadirPts = Vector{Vector{Float64}}() ; return fathomed 

                # don't do EPB branching if the local nadir point is worse than the LBS's nadir point 
                elseif (u.y[2] ≥ ptl.y[2] && u.y[1] ≥ ptr.y[1]) 
                    node.localNadirPts = Vector{Vector{Float64}}() ; return fathomed   

                else 
                    push!(node.localNadirPts, u.y) #; push!(dist_naditPt, dist_ratio(worst_nadir_pt, u.y, ideal_pt))
                end 
            else
                return fathomed
            end
            
        end

        if !compared && (u.y[1] ≥ ptl.y[1] && u.y[2] ≥ ptr.y[2] )
            if EPB node.localNadirPts = Vector{Vector{Float64}}() end               # no need to (extended) pareto branching
            return false
        end
    end

    #todo : indicator not good ...
    # if sum(dist_naditPt)/length(dist_naditPt) > 1/2
    #     node.localNadirPts = Vector{Vector{Float64}}() ; return false
    # end

    return fathomed
end

#-----------------------------
#------ no good indicator ----
#-----------------------------

"""
Return the orthogonal distance between a point `p` and a segment defined by two points `extl` and `extr`.
"""
function orthogonal_dist(p::Vector{Float64}, extl::Vector{Float64}, extr::Vector{Float64})
    return abs((extr[1] - extl[1]) * (extl[2] - p[2]) - (extl[1] - p[1]) * (extr[2] - extl[2])) / sqrt((extr[1] - extl[1])^2 + (extr[2] - extl[2])^2) 
end


"""
Given a non-dominated nadir point, return `true` if the decider EP - branch on it, considering the average distance of nadir point from LBS.
"""
function dist_ratio( worst_nadir_pt::Vector{Float64}, nadir_pt::Vector{Float64}, ideal_pt::Vector{Float64})
    dist_LBS = sqrt((worst_nadir_pt[1] - ideal_pt[1])^2 + (worst_nadir_pt[2] - ideal_pt[2])^2) 
    dist_nadir = sqrt((nadir_pt[1] - ideal_pt[1])^2 + (nadir_pt[2] - ideal_pt[2])^2) 
    return dist_nadir/dist_LBS
end