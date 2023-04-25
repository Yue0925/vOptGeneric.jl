# This file contains functions related to node fathoming.

include("cuttingPlanes.jl")

# global TOL 

"""
Given a point `x_star`, iterate all valid cuts of parent node and stock the references 
in the current node.
"""
function loadingCutInPool(node::Node, pb::BO01Problem)
    if isRoot(node) return end 

    l = 1 ; LBS = node.RBS.natural_order_vect.sols

    while l ≤ length(LBS)
        if LBS[l].is_binary 
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
                if r > length(LBS) || LBS[r].is_binary continue end

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

        # step 2 : add valid cuts constraints then re-optimize 
        start_processing = time()
        loadingCutInPool( node, pb)         # complexity O(pt ⋅ cuts)
        pb.info.cuts_infos.times_add_retrieve_cuts += (time() - start_processing)

        pruned = MP_cutting_planes(node, pb, incumbent, round_results, verbose ; args...)

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

    # println("-------------------")
    # println("node $(node.num)")
    # println("-------------------")
    # # todo : local heuristic search
    # if true #and() > 0.5# true #node.depth < 10 
    #     start = time() 
    #     U_newfea = feasPumingJumping(node, pb, incumbent ; verbose=false)
    #     println("|U_newfea| = $(length(U_newfea.sols))")
    #     println("|incumbent| before = $(length(incumbent.natural_order_vect.sols))")
    #     for s in U_newfea.sols
    #         push!(incumbent.natural_order_vect, s, filtered=true)
    #     end
    #     pb.info.update_incumb_time += (time() - start) 
    #     println("|incumbent| after = $(length(incumbent.natural_order_vect.sols))")
    # end
    # todo 
    # @info "node $(node.num)  depth $(node.depth) \t |LBS| = $(length(node.RBS.natural_order_vect)) "
    
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

    if pb.param.root_relax 
        pb.info.update_incumb_time += (time() - start) ; return false 
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
    # @assert length(incumbent.natural_order_vect) > 1 "`getNadirPoints` requires at least two upper bounds in incumbent list."

    if length(incumbent.natural_order_vect) == 1 return incumbent.natural_order_vect end 

    for i = 1:length(incumbent.natural_order_vect)-1

        push!(nadir_pts, Solution(
            Vector{Vector{Float64}}(),
            [incumbent.natural_order_vect.sols[i].y[1],
            incumbent.natural_order_vect.sols[i+1].y[2]
            ],
            true, Vector{Float64}(), 0.0), filtered=true
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

    # if there exists an upper bound u s.t. u≦l
    function weak_dom(l)
        for u ∈ incumbent.natural_order_vect.sols
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
        l = node.RBS.natural_order_vect.sols[1]
        return weak_dom(l)
    end

    # ----------------------------------------------
    # if the LBS consists of segments
    # ----------------------------------------------
    # two extreme points of LBS
    ptl = node.RBS.natural_order_vect.sols[1] ; ptr = node.RBS.natural_order_vect.sols[end]

    # Case 1 :  if only one feasible point in UBS 
    if length(incumbent.natural_order_vect) == 1 
        u = incumbent.natural_order_vect.sols[1]
        if u.y[1] < ptr.y[1] && u.y[2] < ptl.y[2]
            return true
        else
            return false
        end
    end

    # Case 2 : otherwise, do the pairwise comparison of the local nadir points with LBS  
    nadir_pts = getNadirPoints(incumbent) # , ptl, ptr

    # test range condition necessary 1 : LBS ⊆ UBS 
    u_l = incumbent.natural_order_vect.sols[1] ; u_r = incumbent.natural_order_vect.sols[end]

    sufficient = (u_l.y[2] < ptl.y[2] && u_r.y[1] < ptr.y[1])

    if !sufficient return false end

    # test condition necessary 2 : LBS ≤/dominates UBS 
    fathomed = true# ; dist_naditPt = Vector{Float64}()

    # iterate of all local nadir points
    for u ∈ nadir_pts.sols
        existence = false ; compared = false

        # case 1 : if u is dominates the ideal point of LBS 
        if u.y[1] < ptr.y[1] && u.y[2] < ptl.y[2]
            return true
        end

        # case 3 : complete pairwise comparison
        for i=1:length(node.RBS.natural_order_vect)-1              # ∀ segment l ∈ LBS 

            sol_l = node.RBS.natural_order_vect.sols[i] ; sol_r = node.RBS.natural_order_vect.sols[i+1]

            if (u.y[1] > sol_l.y[1] || u.y[1] < sol_r.y[1]) && (u.y[2] > sol_r.y[2] || u.y[2] < sol_l.y[2])
                continue
            end
            
            λ = [sol_r.y[2] - sol_l.y[2], sol_l.y[1] - sol_r.y[1]]      # normal to the segment

            compared = true

            if λ'*u.y < λ'*sol_r.y#&& λ'*u.y < λ'*sol_l.y
                existence = true ; break
            end
        end
        
        # case 4 : condition dominance violated, then stock the non-dominated local nadir pts to prepare EPB
        if compared && !existence 
            fathomed = false
            if EPB
                if !isRoot(node) && (u.y in node.pred.localNadirPts || u.y == node.pred.nadirPt || u.y == node.nadirPt)    # the current local nadir pt is already branched 
                    node.localNadirPts = Vector{Vector{Float64}}() ; return fathomed 
                    # nothing 
                else 
                    push!(node.localNadirPts, u.y) #; push!(dist_naditPt, dist_ratio(worst_nadir_pt, u.y, ideal_pt))
                end 
            else
                return fathomed
            end
            
        end

        if !compared && (u.y[1] ≥ ptr.y[1] && u.y[2] ≥ ptl.y[2] )
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


"""
Dominance test designed for the LBS that is the convex intersection of the set of lines passing lower bound and perpendicular to it's normal.

    Return `True` if this node is pruned by dominance.
"""
function fullyExplicitDominanceTestByNormal(node::Node, incumbent::IncumbentSet, worst_nadir_pt::Vector{Float64}, EPB::Bool)
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
        l = node.RBS.natural_order_vect.sols[1]
        return weak_dom(l)
    end

    # ----------------------------------------------
    # if the LBS consists of segments
    # ----------------------------------------------
    # two extreme points of LBS
    ptl = node.RBS.natural_order_vect.sols[1] ; ptr = node.RBS.natural_order_vect.sols[end]

    # Case 1 :  if only one feasible point in UBS 
    if length(incumbent.natural_order_vect) == 1 
        u = incumbent.natural_order_vect.sols[1]
        if u.y[1] < ptr.y[1] && u.y[2] < ptl.y[2]
            return true
        else
            return false
        end
    end

    # test range condition necessary 1 : LBS ⊆ UBS 
    u_l = incumbent.natural_order_vect.sols[1] ; u_r = incumbent.natural_order_vect.sols[end]
    # u_l = nadir_pts.sols[1] ; u_r = nadir_pts.sols[end]

    sufficient = (u_l.y[2] < ptl.y[2] && u_r.y[1] < ptr.y[1])
    if !sufficient return false end

    fathomed = true 
    # iterate of all local nadir points
    for u ∈ nadir_pts.sols
        in_polygone = true

        # case 1 : if u is dominates the ideal point of LBS 
        if u.y[1] < ptr.y[1] && u.y[2] < ptl.y[2]
            return true
        end

        if u.y[2] < ptl.y[2] || u.y[1] < ptr.y[1]
            continue
        end

        # case 3 : complete pairwise comparison
        for sol in node.RBS.natural_order_vect.sols # i=1:length(node.RBS.natural_order_vect)              # ∀ segment l ∈ LBS 
            # #todo : ignore intersection pt 
            # if length(sol.xEquiv) == 0 continue end

            λ = sol.λ

            if λ[1] != 0.0 &&  λ[2] != 0.0 && λ'*u.y < λ'*sol.y # strictly inferior : case limit 
                in_polygone = false ; break
            end
        end

        # the end of comparison 
        if in_polygone
            fathomed = false
            if EPB
                if !isRoot(node) && (u.y in node.pred.localNadirPts || u.y == node.pred.nadirPt || u.y == node.nadirPt)    # the current local nadir pt is already branched 
                    node.localNadirPts = Vector{Vector{Float64}}() ; return fathomed 

                elseif (u.y[1] ≥ ptl.y[1] && u.y[2] ≥ ptr.y[2])
                    node.localNadirPts = Vector{Vector{Float64}}() ; return fathomed   
                else 
                    push!(node.localNadirPts, u.y)
                end 
            else
                return fathomed
            end
        end
    end

    return fathomed
end




# -------------------------------------------
# -------------------------------------------
# -------------------------------------------
"""
Dominance test designed for the LBS connected by consecutive local ideal points.

    Return `True` if this node is pruned by dominance.
"""
function fullyExplicitDominanceTestNonConvex(node::Node, incumbent::IncumbentSet, worst_nadir_pt::Vector{Float64}, EPB::Bool)::Bool
    @assert length(node.RBS.natural_order_vect) > 0 "relaxed bound set is empty for node $(node.num)"

    # we can't compare the LBS and UBS if the incumbent set is empty
    if length(incumbent.natural_order_vect) == 0 return false end

    # if there exists an upper bound u s.t. u≦l
    function weak_dom(l)
        for u ∈ incumbent.natural_order_vect.sols
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
        l = node.RBS.natural_order_vect.sols[1]
        return weak_dom(l)
    end

    # ----------------------------------------------
    # if the LBS consists of segments
    # ----------------------------------------------
    # two extreme points of LBS
    ptl = node.RBS.natural_order_vect.sols[1] ; ptr = node.RBS.natural_order_vect.sols[end]

    # Case 1 :  if only one feasible point in UBS 
    if length(incumbent.natural_order_vect) == 1 
        u = incumbent.natural_order_vect.sols[1]
        if u.y[1] < ptr.y[1]  && u.y[2] < ptl.y[2] 
            return true
        else
            return false
        end
    end

    # Case 2 : otherwise, do the pairwise comparison of the local nadir points with LBS  
    nadir_pts = getNadirPoints(incumbent) # , ptl, ptr

    # test range condition necessary 1 : LBS ⊆ UBS 
    u_l = incumbent.natural_order_vect.sols[1] ; u_r = incumbent.natural_order_vect.sols[end]

    sufficient = (u_l.y[2] < ptl.y[2] && u_r.y[1] < ptr.y[1])

    if !sufficient return false end

    # (ideal_pts, fathomed) = getIdealPoints(node.RBS, u_l, u_r) ; if fathomed return true end

    # test condition necessary 2 : LBS ≤/dominates UBS 
    fathomed = true 

    # iterate of all local nadir points
    for u ∈ nadir_pts.sols
        existence = false ; compared = false

        # case 1 : if u is dominates the ideal point of LBS 
        if u.y[1] < ptr.y[1] && u.y[2] < ptl.y[2]
            return true
        end

        # case 3 : complete pairwise comparison
        for i=1:length(node.RBS.natural_order_vect)-1              # ∀ segment l ∈ LBS 

            sol_l = node.RBS.natural_order_vect.sols[i] ; sol_r = node.RBS.natural_order_vect.sols[i+1]

            if (u.y[2] > sol_r.y[2] || u.y[2] < sol_l.y[2] )
                continue
            end
            
            compared = true

            if u.y[1] < sol_r.y[1] 
                existence = true ; break
            end
        end


        # if !compared && u.y[2] <= ptl.y[2] && u.y[1] >= ptl.y[1]
        #     if EPB
        #         node.localNadirPts = Vector{Vector{Float64}}()
        #     end
        #     return false
        # end
        
        # case 4 : condition dominance violated, then stock the non-dominated local nadir pts to prepare EPB
        if compared && !existence 
            fathomed = false
            if EPB
                if !isRoot(node) && (u.y in node.pred.localNadirPts || u.y == node.pred.nadirPt || u.y == node.nadirPt)    # the current local nadir pt is already branched 
                    node.localNadirPts = Vector{Vector{Float64}}() ; return fathomed 
                    # nothing
                else 
                    push!(node.localNadirPts, u.y)
                end 
            else
                return fathomed
            end
            
        end
 
        if !compared && u.y[2] ≥ ptr.y[2]
            if  ptr.y[1] <= u.y[1] <= ptl.y[1]
                fathomed = false
                if EPB
                    if !isRoot(node) && (u.y in node.pred.localNadirPts || u.y == node.pred.nadirPt || u.y == node.nadirPt)    # the current local nadir pt is already branched 
                        node.localNadirPts = Vector{Vector{Float64}}() ; return fathomed 
                        # nothing
                    else 
                        push!(node.localNadirPts, u.y)
                    end 
                else
                    return fathomed
                end
            elseif u.y[1] > ptl.y[1] 
                if EPB node.localNadirPts = Vector{Vector{Float64}}() end               # no need to (extended) pareto branching
                return false
            end
        end
    end

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