"""
    This file contains objects types used for BO01BB algorithm.
    Consider the following multi-objective program :

    min     f1(⋅), f2(⋅)
    s.t.      Ax ≤ b
              x∈{0,1}
"""

@enum PrunedType NONE INFEASIBILITY OPTIMALITY DOMINANCE

TOL = 1e-4
include("cutPool.jl")

"""
Storage the total number of cuts applied etc...
"""
mutable struct CutsInfo
    ite_total::Int64
    cuts_applied::Int64
    sp_cuts::Int64
    mp_cuts::Int64
    cuts_total::Int64
    times_calling_dicho::Float64
    times_calling_separators::Float64
    times_oper_cutPool::Float64
    times_total_for_cuts::Float64
    times_add_retrieve_cuts::Float64
end

function CutsInfo()
    return CutsInfo(0, 0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0)
end


function Base.:show(io::IO, cinfo::CutsInfo)
    println(io, "\n\n # ----------- info about cuts : \n", 
    "ite_total = $(cinfo.ite_total) \n",
    "cuts_applied = $(cinfo.cuts_applied) \n",
    "sp_cuts = $(cinfo.sp_cuts) \n",
    "mp_cuts = $(cinfo.mp_cuts) \n", 
    "cuts_total = $(cinfo.cuts_total) \n", 
    "times_calling_dicho = $(cinfo.times_calling_dicho) \n",
    "times_calling_separators = $(cinfo.times_calling_separators) \n", 
    "times_oper_cutPool = $(cinfo.times_oper_cutPool) \n", 
    "times_total_for_cuts = $(cinfo.times_total_for_cuts) \n",
    "times_add_retrieve_cuts = $(cinfo.times_add_retrieve_cuts) \n"
    )
end

"""
Storage the statistics information of the BO01BB algorithm.
"""
mutable struct StatInfo
    total_times::Float64
    nb_nodes::Int64
    nb_nodes_pruned::Int64
    Gap::Float64
    relaxation_time::Float64
    test_dom_time::Float64
    update_incumb_time::Float64
    tree_size::Float64
    cp_activated::Bool
    cuts_infos::CutsInfo
    nb_nodes_EPB::Int64
    nb_nodes_VB::Int64
    root_relax::Bool
    TO::Bool
    rootLBS::Int64
    LBSexhaustive::Bool
    # status::MOI.TerminationStatusCode 
end

function StatInfo()
    return StatInfo(0.0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, false, CutsInfo(), 0, 0, false, false, 0, true)
end

function Base.:show(io::IO, info::StatInfo)
    println(io, " # informations of B&B algorithm : \n",
        "total_times_used = $(info.total_times) \n",
        "total_nodes = $(info.nb_nodes) \n",
        "pruned_nodes = $(info.nb_nodes_pruned) \n",
        "GAP = $(info.Gap) \n",
        "relaxation_time = $(info.relaxation_time) \n",
        "test_dominance_time = $(info.test_dom_time) \n",
        "update_incumbent_time = $(info.update_incumb_time) \n",
        "tree_size = $(info.tree_size) \n",
        "nb_nodes_EPB = $(info.nb_nodes_EPB) \n",
        "nb_nodes_VB = $(info.nb_nodes_VB) \n",
        "root_relax = $(info.root_relax) \n",
        "TO = $(info.TO) \n", 
        "rootLBS = $(info.rootLBS) \n", 
        "LBSexhaustive = $(info.LBSexhaustive)"
    )
    if info.cp_activated println(io, info.cuts_infos) end
end


"""
Some parameters used in the B&B algorithm.
"""
mutable struct BBparam
    time_limit::Float64     # time limit for B&B algo
    traverse::Symbol        # traverse strategy such as dfs, bfs...
    branching::Symbol       # branching strategy
    cp_activated::Bool     # if apply cuts at each node
    EPB::Bool               # if consider the EP branching 
    root_relax::Bool        # use the relaxation value at root in an IP solver 
end

function BBparam()
    return BBparam(300, :bfs, :arbitrary, false, false, false)
end


"""
Storage the components definning a bi-objective 0-1 linear program.
"""
mutable struct BO01Problem
    varArray::Array{JuMP.VariableRef}
    varArray_copied::Array{JuMP.VariableRef}
    m::JuMP.Model
    lp_copied::JuMP.Model
    param::BBparam
    info::StatInfo
    A::Matrix{Float64}
    b::Vector{Float64}
    c::Matrix{Float64}
    SP_CG_model_defined::Bool
    SP_CG_model::JuMP.Model
    SP_CG_α::Vector{JuMP.VariableRef}
    MP_CG_model_defined::Bool
    MP_CG_model::JuMP.Model
    MP_CG_α::Vector{JuMP.VariableRef}
end


"""
A point `y` in the criteria space may be projected from a set of equivalent solutions `x`.
"""
mutable struct Solution
    xEquiv::Vector{Vector{Float64}}            # a set of equivalent solutions x defining the same y
    y::Vector{Float64}
    is_binary::Bool
    λ::Vector{Float64}
    ct::Float64
end

function Solution()
    return Solution(Vector{Vector{Float64}}(), Vector{Float64}(), false, Vector{Float64}(), Inf)
end

function Solution(x::Vector{Float64}, y::Vector{Float64}, λ::Vector{Float64}=Vector{Float64}(), ct::Float64=Inf)
    is_binary = length(x) >0 ? true : false 
    for i = 1:length(x)
        if !(abs(x[i]-0.0) ≤ TOL || abs(x[i]-1.0) ≤ TOL)
            is_binary = false; break
        end
    end    
    return Solution([ x ], y , is_binary, λ, ct)
end


function roundSol(pb::BO01Problem, sol::Solution)
    if !sol.is_binary return end

    x_ = round.(sol.xEquiv[1]) ; 
    for x_i in x_
        if x_i < 0.0 && x_i > 1.0
            println("!!! error rounding x = $(sol.xEquiv[1]) \n !!! x_rounded = $x_")
        end
    end

    for i in size(pb.A, 1)
        if pb.A[i, :]'* x_ > pb.b[i] return end 
    end

    # todo : verify deepcopy ?
    sol.xEquiv[1] = x_; sol.is_binary = isBinary(x_)
    sol.y = [x_'*pb.c[1, 2:end] + pb.c[1, 1] ,x_'*pb.c[2, 2:end] + pb.c[2, 1]  ]
end


"""
Given a vector `x`, return true if `x` is approximately binary.
"""
function isBinary(x::Vector{Float64})
    if length(x) == 0 return false end 
    for i in 1:length(x)
        if !(abs(x[i]-0.0)<=TOL || abs(x[i]-1.0) <=TOL)
            return false
        end
    end
    return true
end

"""
Add an equivalent solution associated to point `y`. 
"""
function addEquivX(sol::Solution, x::Vector{Float64})
    # @assert length(x) > 0 "x cannot be empty"
    if length(x) == 0 return end 

    # check if x is approximately binary
    if isBinary(x)
        sol.is_binary = true
        # sol.y = round.(sol.y, digits = 4)
        push!(sol.xEquiv, x) # round.(x, digits = 4)
    else
        push!(sol.xEquiv, x)
    end
    
end

function addEquivX(sol::Solution, vecX::Vector{Vector{Float64}})
    for x in vecX
        for added_x in sol.xEquiv
            if x == added_x return end
        end
        addEquivX(sol, x)
    end
end

"""
Overload operators for the dominance order between two solutions.
"""
function Base.:show(io::IO, s::Solution)
    println(io, "Solution( \n |xEquiv| = ", length(s.xEquiv), 
    "\t is_binary ? ", s.is_binary,
    "\n y = ", s.y,
    "\n λ = ", s.λ, "\t ct = ", s.ct, 
    "\n xEquiv = ", s.xEquiv, 
    " )")
end


function Base.:<=(a::Solution, b::Solution)
    @assert length(a.y) > 0
    @assert length(b.y) > 0
    return a.y[1] ≤ b.y[1] && a.y[2] ≤ b.y[2]
end

function Base.:<(a::Solution, b::Solution)
    @assert length(a.y) > 0
    @assert length(b.y) > 0
    return a.y[1] < b.y[1]  && a.y[2] < b.y[2]
end

function Base.:>=(a::Solution, b::Solution)
    @assert length(a.y) > 0
    @assert length(b.y) > 0
    return a.y[1] ≥ b.y[1] && a.y[2] ≥ b.y[2]
end

function Base.:>(a::Solution, b::Solution)
    @assert length(a.y) > 0
    @assert length(b.y) > 0
    return a.y[1]>b.y[1] && a.y[2]>b.y[2]
end

function Base.:(==)(a::Solution, b::Solution)
    @assert length(a.y) > 0
    @assert length(b.y) > 0
    return abs(a.y[1] -b.y[1]) ≤ TOL && abs(a.y[2] - b.y[2]) ≤ TOL
end

function Base.:(!=)(a::Solution, b::Solution)
    @assert length(a.y) > 0
    @assert length(b.y) > 0
    return abs(a.y[1] - b.y[1]) > TOL || abs(a.y[2] - b.y[2]) > TOL
end

function Base.isequal(a::Solution, b::Solution) # for hash table
    @assert length(a.y) > 0
    @assert length(b.y) > 0
    return a == b
end


"""
Return `true` if solution `a` (weakly) dominates sobution `b`; `false` otherwise.
"""
function dominate(a::Solution, b::Solution)
    return a ≤ b && a != b
end

"""
Return `true` if solution `a` strictly dominates sobution `b`; `false` otherwise.
"""
function strictDominate(a::Solution, b::Solution)
    return a ≤ b && a.y[1] ≠ b.y[1] && a.y[2] ≠ b.y[2] 
end

"""
straight line       a x + b y = c 
            <=>   -Δz2 z1 + Δz1 z2 = Δz1 c 

a = -Δz2 = -λ1   b = Δz1 = -λ2    ct = Δz1 c 
        # todo (attention) : λ always absolute value 
"""
function updateCT(s::Solution)
    if length(s.λ) ≠ 2 || s.λ == [0.0, 0.0] return end  # s.ct ≠ Inf || 

    s.ct = -s.λ[1] * s.y[1] - s.λ[2] * s.y[2]
end


"""
A vector of solutions in the natrual order (from left to right in bi-objective space). 
"""
mutable struct NaturalOrderVector
    sols::Vector{Solution}
end

function NaturalOrderVector()
    return NaturalOrderVector(Vector{Solution}())
end

function Base.length(v::NaturalOrderVector)
    return length(v.sols)
end

function Base.:show(io::IO, nov::NaturalOrderVector)
    println(io, "NaturalOrderVector[")
    for sol in nov.sols
        print(io, sol)
    end
    println("]")
end


"""
Push a solution into a vector of natrual ordered solutions, return `true` if it is successfully added;
or `false`, if it is weakly dominated by one (or more) solution(s) in the vector. 

In case of successfully added and `filtered=true` (by defaut false), delete the old solutions that are weakly dominated by the new one.

Return 

    -   the position `idx` successfully inserted,
# todo : return 
    -   in case of same point with different ct , update ct and return `idx`
    -   -1 : the same point with exactly ct
            OR  the inconming point is dominated 
"""
function Base.push!(natural_sols::NaturalOrderVector, sol::Solution; filtered::Bool=false)::Tuple{Int,Bool}
    # sol.y = round.(sol.y, digits = 4) ; 
    idx = -1

    # add s directly if sols is empty
    if length(natural_sols) == 0
        push!(natural_sols.sols, sol) ; return 1, true
    end

    # a binary/dichotomy search finds the location to insert 
    l = 1; r = length(natural_sols); m = 0
    while l ≤ r
        m = Int(floor((l+r)/2))

        if abs(sol.y[2] - natural_sols.sols[m].y[2]) ≤ TOL
            # in case of the approximately equality on the first objective, compare the second obj
            if abs(sol.y[1] - natural_sols.sols[m].y[1]) ≤ TOL
                
                # if natural_sols.sols[m].is_binary
                #     if sol.is_binary
                #         addEquivX(natural_sols.sols[m], sol.xEquiv)
                #         natural_sols.sols[m].λ = deepcopy(sol.λ) ; updateCT(natural_sols.sols[m])
                #     end
                # else
                #     # replaced by the new solution directly 
                #     if length(sol.xEquiv[1])>0 
                #         natural_sols.sols[m].xEquiv = deepcopy(sol.xEquiv)
                #         natural_sols.sols[m].y = deepcopy(sol.y)
                #         natural_sols.sols[m].is_binary = sol.is_binary
                #         natural_sols.sols[m].λ = deepcopy(sol.λ)
                #         natural_sols.sols[m].ct = sol.ct
                #     end
                # end                

                # ------------------------------------------------
                # # todo : different λ same pente
                if length(sol.xEquiv[1])>0 && length(sol.λ) == 2 && length(natural_sols.sols[m].λ) == 2 

                    if abs(sol.λ[1]- natural_sols.sols[m].λ[1])>TOL || abs(sol.λ[2]- natural_sols.sols[m].λ[2])>TOL

                        natural_sols.sols[m].λ = deepcopy(sol.λ)
                        updateCT(natural_sols.sols[m])
                    end
                end
                addEquivX(natural_sols.sols[m], sol.xEquiv) 
                # -------------------------------------------------

                return m, false
            elseif sol.y[1] > natural_sols.sols[m].y[1]
                l = m+1
            elseif sol.y[1] < natural_sols.sols[m].y[1]
                r  = m-1
            end
        elseif sol.y[2] < natural_sols.sols[m].y[2]
            l = m+1
        elseif sol.y[2] > natural_sols.sols[m].y[2]
            r = m-1
        end
    end

    if r==0 # insert at the top
        natural_sols.sols = vcat([sol], natural_sols.sols)
        m = 1
    elseif l==length(natural_sols)+1 # insert at the bottom
        push!(natural_sols.sols, sol)
        m = l
    else # insert at m
        m = m > Int(floor((l+r)/2)) ? m : m+1
        natural_sols.sols = vcat(vcat(natural_sols.sols[1:m-1], sol), natural_sols.sols[m:end])
    end
    idx = m 

    # find points weakly dominated by the new point and delete it/them
    if filtered
        i = 1
        while i < length(natural_sols.sols)
            j = i+1
            while j<= length(natural_sols.sols)
                if dominate(natural_sols.sols[i], natural_sols.sols[j])
                    deleteat!(natural_sols.sols, j)
                    if j == m idx= -1 elseif m > j idx -= 1 end 
                    
                elseif dominate(natural_sols.sols[j], natural_sols.sols[i])
                    deleteat!(natural_sols.sols, i)
                    if i == m idx= -1 elseif m > i idx -= 1 end 
                    j -= 1 ; break
                else
                    j += 1
                end
            end
            if j > length(natural_sols.sols) i += 1 end 
        end
    end
    return idx< 1 ? (idx, false) : (idx, true)
end


"""
The relaxed bound set consists of segments and natural ordered solutions. 

Since all solutions are natural ordered, segments is defined by a dictionary where each point solution s is associated
with a boolean if s and the next/adjacent point in right forms a segment.
"""
mutable struct RelaxedBoundSet
    natural_order_vect::NaturalOrderVector
end

function RelaxedBoundSet()
    return RelaxedBoundSet(NaturalOrderVector()) 
end


"""
The incumbent set consists in feasible solutions.
"""
mutable struct IncumbentSet
    natural_order_vect::NaturalOrderVector
end

function IncumbentSet()
    return IncumbentSet(NaturalOrderVector())
end
