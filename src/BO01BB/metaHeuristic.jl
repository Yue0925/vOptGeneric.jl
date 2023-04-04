include("BBtree.jl")
include("../algorithms.jl")
include("struct.jl")
include("separators.jl")
include("cutPool.jl")

using JuMP 
using Random

const IterLimit = 5
const EPS = 0.001






"""
Given a fractional solution `l`, return a rounded jumped interger solution `s`.
"""
function rounding_jumping(l::Solution, pb::BO01Problem, assign::Dict{Int64, Int64}; proba::Float64=0.5)::Solution
    x̄ = [-1 for _ in 1:length(l.xEquiv[1])]
    rhs = deepcopy(pb.b)

    # don't violate the vars assignment
    for (j, bnd) in assign
        if bnd == 0
            x̄[j] = 0
        else
            x̄[j] = 1
            for i in 1:size(pb.A, 1)
                rhs[i] -= pb.A[i, j]
            end
        end
    end

    # ----------------------------------
    # rounding to 1 with a probability
    # ----------------------------------
    idx = shuffle!( [j for j in 1:length(l.xEquiv[1]) if x̄[j] == -1 ] )
    for j in idx
        if abs(l.xEquiv[1][j] - 0) ≤ EPS
            x̄[j] = 0
        elseif abs(l.xEquiv[1][j] - 1) ≤ EPS && rand() ≥ proba
            for i in 1:size(pb.A, 1)
                if pb.b[i] *(rhs[i] - pb.A[i, j]) < 0
                    continue
                end
            end

            x̄[j] = 1
            for i in 1:size(pb.A, 1)
                rhs[i] -= pb.A[i, j]
            end
        end
    end

    # ---------------------------
    # jumping 
    # ---------------------------
    while -1 in x̄
        # {j => [viol(xj=0), viol(xj=1)]}
        viol_reduct = Dict{Int64, Vector{Float64}}( j => [sum(abs.(rhs)), sum(abs.(rhs .- pb.A[:, j]))] for j in 1:length(l.xEquiv[1]) if x̄[j]==-1)
        # {j => [min viol, obj coeff]}
        preference = Dict{Int64, Vector{Float64}}(j => [minimum(viol_reduct[j]), l.λ'*pb.c[:, j+1]] for j in 1:length(l.xEquiv[1]) if x̄[j]==-1)
        # sorted idx
        idx_sorted = sort!([j for j in 1:length(l.xEquiv[1]) if x̄[j]==-1], by= v ->(preference[v][2], preference[v][1]) )

        j = idx_sorted[1]
        if preference[j][1] == viol_reduct[j][1]
            x̄[j] = 0
        else
            x̄[j] = 1
            for i in 1:size(pb.A, 1)
                rhs[i] -= pb.A[i, j]
            end
        end
    end

    y = [x̄'* pb.c[1, 2:end] + pb.c[1, 1], x̄'* pb.c[2, 2:end] + pb.c[2, 1]]
    return Solution([x̄], y, true, Vector{Float64}())
end


"""
Return `true` if `s` is a feasible solution.
"""
function isFeasible(s::Solution, pb::BO01Problem)::Bool
    for i in 1:size(pb.A, 1)
        if s.xEquiv[1]'*pb.A[i, :] > pb.b[i]
            return false
        end
    end
    return true 
end


"""
Find all local nadir points dominated by the given lower bound. 
"""
function nadirPtsZone(incumbent::IncumbentSet, l::Solution)
    nadir_pts = NaturalOrderVector()

    if length(incumbent.natural_order_vect) == 1 
        if dominate(l, incumbent.natural_order_vect.sols[1])
            return incumbent.natural_order_vect 
        end
    end 

    for i = 1:length(incumbent.natural_order_vect)-1
        u = Solution(
            Vector{Vector{Float64}}(),
            [incumbent.natural_order_vect.sols[i].y[1],
            incumbent.natural_order_vect.sols[i+1].y[2]
            ],
            true, Vector{Float64}() 
        )

        if dominate(l, u) push!(nadir_pts, u, filtered=true) end 
    end

    return nadir_pts
end


"""
Flip some vars with high difference between fractional sol and integer sol, and closet to 0.5.
"""
function flip(s̃::Solution, s̄::Solution, assign::Dict{Int64, Int64}, pb::BO01Problem)::Solution
    n = length(s̄.xEquiv[1])
    forbid = keys(assign) ; idx = [j for j in 1:n if !(j in forbid)]

    T = length(idx)/4 ; nb_flipped = rand(round(Int64, T/2):round(Int64, T*3/4))
    x = deepcopy(s̄.xEquiv[1][:])

    score = Dict{Int64, Vector{Float64}}(
            j => [-abs(s̃.xEquiv[1][j] - s̄.xEquiv[1][j]), abs(s̃.xEquiv[1][j] - 0.5) ] for j in idx
    )
    idx_sorted = sort!(idx, by=v -> (score[v][1], score[v][2]) )

    while nb_flipped > 0
        nb_flipped -= 1
        j = popfirst!(idx_sorted)
        if s̄.xEquiv[1][j] == 1
            x[j] = 0
        else
            x[j] = 1
        end
    end

    y = [x'* pb.c[1, 2:end] + pb.c[1, 1], x'* pb.c[2, 2:end] + pb.c[2, 1]]
    return Solution([x], y, true, Vector{Float64}())
end


"""
Return `true` if vector `H` contains solution `s`.
"""
function contains(H::Vector{Solution}, s::Solution)::Bool
    for s_bis in H
        if s_bis == s return true end 
    end
    return false
end


function Δ_opt(pb::BO01Problem, s::Solution, nadir_pts::NaturalOrderVector)::Solution
    # ---------------
    col = length(pb.varArray_copied) ; row = length(pb.b)
    idx0 = [j for j in 1:col if s.xEquiv[1][j] == 0] ; idx1 = [j for j in 1:col if s.xEquiv[1][j] == 1]

    JuMP.set_objective(pb.lp_copied, MOI.MIN_SENSE,
         sum(pb.varArray_copied[j] for j in idx0 ; init=0) + sum(1-pb.varArray_copied[j] for j in idx1 ; init=0)
    )

    # pareto bound 
    ctr_symbol = []
    for u in nadir_pts.sols
        ctr = JuMP.@constraint(pb.lp_copied, pb.varArray_copied'* pb.c[1, 2:end] + pb.c[1, 1] ≤ u.y[1]) ; push!(ctr_symbol, ctr)
        ctr = JuMP.@constraint(pb.lp_copied, pb.varArray_copied'* pb.c[2, 2:end] + pb.c[2, 1] ≤ u.y[2]) ; push!(ctr_symbol, ctr)
    end

    JuMP.optimize!(pb.lp_copied, ignore_optimize_hook=true) ; status = JuMP.termination_status(pb.lp_copied)

    x = JuMP.value.(pb.varArray_copied)
    y = [x'* pb.c[1, 2:end] + pb.c[1, 1], x'* pb.c[2, 2:end] + pb.c[2, 1]]

    # retrieve new ctr added 
    if length(ctr_symbol) > 0
        for ctr in ctr_symbol
            if JuMP.is_valid(pb.lp_copied, ctr)
                JuMP.delete(pb.lp_copied, ctr) ; JuMP.unregister(pb.lp_copied, :ctr)
            end
        end
    end

    return Solution(x, y, Vector{Float64}())

end


"""
Return non-dominated new feasible solutions. 
"""
function feasPumingJumping(node::Node, pb::BO01Problem, incumbent::IncumbentSet; verbose::Bool=false)
    LBS = node.RBS.natural_order_vect.sols
    U_newfea = NaturalOrderVector()

    # ∀ l lower bound
    for l in LBS
        if l.is_binary continue end 

        s̄ = rounding_jumping(l, pb, node.assignment)
        if isFeasible(s̄, pb) push!(U_newfea, s̄, filtered=true) end

        # # pumping only at root  
        # if node.depth >0 continue end

        # H = Vector{Solution}()
        # push!(H, s̄) ; zone = nadirPtsZone(incumbent, l)

        # iter = 0
        # while iter < IterLimit
        #     iter += 1

        #     s̃ = Δ_opt(pb, s̄, zone) ; s̃.λ = l.λ

        #     if s̃.is_binary 
        #         push!(U_newfea, s̃, filtered=true) ; break  
        #     end 

        #     s̄ = rounding_jumping(s̃, pb, node.assignment)

        #     if contains(H, s̄)
        #         s̄ = flip(s̃, s̄, node.assignment, pb)
        #         verbose && println("flip ...")
        #     end
        #     push!(H, s̄)
        #     if isFeasible(s̄, pb) push!(U_newfea, s̄, filtered=true) end 
        # end
    end

    return U_newfea
end

