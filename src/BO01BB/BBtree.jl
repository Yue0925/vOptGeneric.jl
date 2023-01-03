using JuMP

include("struct.jl")
include("cutPool.jl")

"""
Definition of the node object in B&B tree.
"""
mutable struct Node
    num::Int64                    
    depth::Int64                # depth in tree
    pred::Node                  # predecessor
    succs::Vector{Node}         # successors
    var::Int64                  # index of the chosen variable to be split
    var_bound::Int64            # variable bound
    localNadirPts::Vector{Vector{Float64}}      # stok a list of non-dominated nadir pts for EPB 
    EPB::Bool                   # if this node is (extended) pareto branched 
    nadirPt::Vector{Float64}    # if EPB, indicate the pt branched from  
    duplicationBound::Float64   # an additional bound on objective z₁ (maybe redundant) to avoid duplicated search area during EPB 
    RBS::RelaxedBoundSet        # local relaxed bound set    
    activated::Bool             # if the node is active
    pruned::Bool                # if the node is pruned
    prunedType::PrunedType      # if the node is fathomed, restore pruned type
    deleted::Bool               # if the node is supposed to be deleted
    con_cuts::Vector{ConstraintRef}             
    cutpool::CutPool
    assignment::Dict{Int64, Int64}

    Node() = new()

    """
    A complete node constructor.
    """
    function Node(num::Int64, depth::Int64; pred::Node=Node(), succs::Vector{Node}=Vector{Node}(),
        var::Int64=0, var_bound::Int64=0, EPB::Bool=false, nadirPt::Vector{Float64}=Vector{Float64}(), duplicationBound::Float64=Inf)
        n = new()
        n.num = num
        n.depth = depth
    
        n.pred = pred
        n.succs = succs
        n.var = var
        n.var_bound = var_bound

        n.localNadirPts = Vector{Vector{Float64}}()
        n.EPB = EPB
        n.nadirPt = nadirPt
        n.duplicationBound = duplicationBound
    
        n.RBS = RelaxedBoundSet()
        n.activated = true
        n.pruned = false
        n.prunedType = NONE
        n.deleted = false
        # n.objs = Vector{JuMP.GenericAffExpr}()
        # n.cuts_ref = Vector{CutScore}()
        n.con_cuts = Vector{ConstraintRef}()
        n.cutpool = CutPool()
        n.assignment = Dict{Int64, Int64}()
    
        f(t) = nothing # @async println("Finalizing node $(t.num).")
        finalizer(f, n)
    
        # return n
    end
    
end



function Base.:show(io::IO, n::Node)
    println(io, "\n\n # ----------- node $(n.num) : \n", 
    "depth = $(n.depth) \n",
    "pred = $(n.pred.num) \n",
    # "succs = $(n.succs) \n",
    "var[ $(n.var) ] = $(n.var_bound) \n",
    "EPB = $(n.EPB) \n",
    "localNadirPts = $(n.localNadirPts) \n",
    "nadirPt = $(n.nadirPt) \n",
    # "LBS = $(n.RBS.natural_order_vect) \n",
    "activated = $(n.activated) \n",
    "pruned = $(n.pruned) \n",
    "prunedType = $(n.prunedType)"
    )
    print(io, "succs = [ ")
    for s in n.succs print(io, "$(s.num), ") end
    println(io, " ]")

    print(io, "LBS = [ ")
    for s in n.RBS.natural_order_vect.sols
        print(io, "$(s.y) , ")
    end
    println(io, "] ")
end


"""
Return `true` if the given node is the root of a branch-and-bound tree.
"""
function isRoot(node::Node)
    return node.depth == 0 # !isdefined(node, :pred) 
end

"""
Return `true` if the given `node` has activated/non-explored successor(s).
"""
function hasNonExploredChild(node::Node)
    for c in node.succs
        if c.activated return true end 
    end
    return false
end

"""
Delete the given node in B&B tree. (should be a private function)
"""
function Base.delete!(node::Node)
    node.deleted = true ; node = nothing               # remove from the memory
end

"""
Prune the given node in a B&B tree and delete all successors of the pruned node.
"""
function prune!(node::Node, reason::PrunedType)
    node.pruned = true
    node.prunedType = reason
    to_delete = node.succs[:]
    node.succs = Vector{Node}()

    while length(to_delete) > 0
        n = pop!(to_delete)
        to_delete = vcat(to_delete, n.succs[:])
        delete!(n)
    end
end

"""
From the actual node, go up to the root to get the partial assignment of variables.
"""
function getPartialAssign(actual::Node)::Dict{Int64, Int64}
    assignment = Dict{Int64, Int64}() # var index => bound value
    if isRoot(actual) # the actual node is the root 
        return assignment
    end
    predecessor = actual.pred
    if !actual.EPB assignment[actual.var] = actual.var_bound end 

    while !isRoot(predecessor)     
        actual = predecessor ; predecessor = actual.pred
        if !actual.EPB assignment[actual.var] = actual.var_bound end 
    end
    return assignment
end


"""
Going through all the predecessors until the root, add variables or objective bounds branched in the predecessors.

Return a list of objective bounds (symbolic constraint).
"""
function setVarObjBounds(actual::Node, pb::BO01Problem)
    if isRoot(actual) # the actual node is the root 
        return 
    end
    predecessor = actual.pred ; con_cuts = []

    # set actual objective/variable bound
    if actual.EPB
        append!(con_cuts, setObjBound(pb, actual.nadirPt, actual.duplicationBound))
    else
        setVarBound(pb, actual.var, actual.var_bound)
    end

    # set actual objective/variable bounds in predecessors
    while !isRoot(predecessor)     
        actual = predecessor ; predecessor = actual.pred
        if actual.EPB
            append!(con_cuts, setObjBound(pb, actual.nadirPt, actual.duplicationBound))
        else
            setVarBound(pb, actual.var, actual.var_bound)
        end
    end
    return con_cuts
end


"""
Remove variables or objective bounds set in the predecessors.
"""
function removeVarObjBounds(actual::Node, pb::BO01Problem, objcons)
    if isRoot(actual) # the actual node is the root 
        return 
    end
    predecessor = actual.pred
    if !actual.EPB removeVarBound(pb, actual.var, actual.var_bound) end

    while !isRoot(predecessor)     
        actual = predecessor ; predecessor = actual.pred
        if !actual.EPB removeVarBound(pb, actual.var, actual.var_bound) end
    end

    for con in objcons
        if JuMP.is_valid( pb.m, con)
            JuMP.delete( pb.m, con) ; JuMP.unregister( pb.m, :con) # remove the symbolic reference
        end
    end
end


"""
Given a patrial assignment of variables, remove the fixed bounding.
"""
function removeVarBound(pb::BO01Problem, var::Int64, bound::Int64)
    if bound == 0
        JuMP.set_upper_bound(pb.varArray[var], 1)
    elseif bound == 1
        JuMP.set_lower_bound(pb.varArray[var], 0)
    end
end


"""
Given a partial assignment on variables values, add the corresponding bounds.
"""
function setVarBound(pb::BO01Problem, var::Int64, bound::Int64)
    if bound == 0
        JuMP.set_upper_bound(pb.varArray[var], 0)
    elseif bound == 1
        JuMP.set_lower_bound(pb.varArray[var], 1)
    end
end

"""
Given a un-dominated nadir point, add the correspond EPB objective bound.
"""
function setObjBound(pb::BO01Problem, nadirPt::Vector{Float64}, duplication_bound::Float64)
    cons_obj = []
    # for i=1:2 
        con = JuMP.@constraint(pb.m, pb.c[1, 1] + pb.c[1, 2:end]'*pb.varArray ≤ nadirPt[1]) ; push!(cons_obj, con)
        con = JuMP.@constraint(pb.m, pb.c[2, 1] + pb.c[2, 2:end]'*pb.varArray ≤ nadirPt[2]) ; push!(cons_obj, con)
        if duplication_bound != Inf con = JuMP.@constraint(pb.m, pb.c[1, 1] + pb.c[1, 2:end]'*pb.varArray ≥ duplication_bound) ; push!(cons_obj, con) end 
    # end
    return cons_obj
end