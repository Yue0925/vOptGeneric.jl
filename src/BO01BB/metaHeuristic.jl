include("cuttingPlanes.jl")

using Random

const IterLimit = 5
const EPS = 0.001

# todo : new var U_new new feas encountered ; check all bounds cuts for m/m_copied 

"""
Given a fractional solution `l`, return a rounded jumped interger solution `s`.
"""
function rounding_jumping(l::Solution, pb::BO01Problem; proba::Float64=0.5)::Solution
    x = [-1 for _ in 1:length(l.x)]
    rhs[:] = pb.b[:]

    # ----------------------------------
    # rounding to 1 with a probability
    # ----------------------------------
    idx = shuffle!(collect(1:length(l.x)))
    for j in idx
        if abs(l.x[j] - 0) ≤ EPS
            x[j] = 0
        elseif abs(l.x[j] - 1) ≤ EPS && rand() ≥ proba
            x[j] = 1
            for i in 1:size(pb.A, 1)
                rhs[i] -= pb.A[i, j]
            end
        end
    end

    # ---------------------------
    # jumping 
    # ---------------------------
    while -1 in x
        # {j => [viol(xj=0), viol(xj=1)]}
        viol_reduct = Dict{Int64, Vector{Float64}}( j => [sum(abs.(rhs)), sum(abs.(rhs .- pb.A[:, j]))] for j in 1:length(l.x) if x[j]==-1)
        # {j => [min viol, obj coeff]}
        preference = Dict{Int64, Vector{Float64}}(j => [min(viol_reduct[j]), l.λ'*pb.c[:, j+1]] for j in 1:length(l.x) if x[j]==-1)
        # sorted idx
        idx_sorted = sort!([j for j in 1:length(l.x) if x[j]==-1], by= v ->(preference[v][1], preference[v][2]) )

        j = idx_sorted[1]
        if preference[j][1] == viol_reduct[j][1]
            x[j] = 0
        else
            x[j] = 1
            for i in 1:size(pb.A, 1)
                rhs[i] -= pb.A[i, j]
            end
        end
    end

    y = [x'* pb.c[1, 2:end] + pb.c[1, 1], x'* pb.c[2, 2:end] + pb.c[2, 1]]
    return Solution([x], y, true, Vector{Float64}())
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
function flip(s̃::Solution, s̄::Solution)::Solution
    n = length(s̄.xEquiv[1])
    T = n/4 ; nb_flipped = rand(round(Int64, T/2):round(Int64, T*3/4))
    x[:] = s̄.xEquiv[1][:]

    score = Dict{Int64, Vector{Float64}}(
            j => [-abs(s̃.xEquiv[1][j] - s̄.xEquiv[1][j]), abs(s̃.xEquiv[1][j] - 0.5) ] for j in 1:n
    )
    idx_sorted = sort!([j for j in 1:n], by=v -> (score[v][1], score[v][2]) )

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
    # best_bound = objective_bound(m)
    # ctr_bound = JuMP.@constraint(lp_copied, λ1*f1_copied + λ2*f2_copied >= best_bound)
    # JuMP.optimize!(lp_copied, ignore_optimize_hook=true)
    # x_star = JuMP.value.(varArray_copied)

    # if JuMP.is_valid(lp_copied, ctr_bound)
    #     JuMP.delete(lp_copied, ctr_bound) ; JuMP.unregister(lp_copied, :ctr_bound)
    # end

    # ---------------
    col = size(pb.A, 1) ; row = size(pb.A, 2)
    idx0 = [j for j in 1:col if s.xEquiv[1][j] == 0] ; idx1 = [j for j in 1:col if s.xEquiv[1][j] == 1]

    varArray_copied = JuMP.all_variables(pb.lp_copied)
    JuMP.set_objective(pb.lp_copied, MOI.MIN_SENSE,
         sum(varArray_copied[j] for j in idx0) + sum(1-varArray_copied[j] for j in idx1)
    )

    # pareto bound 
    ctr_symbol = []
    for u in nadir_pts.sols
        ctr = JuMP.@constraint(pb.lp_copied, varArray_copied'* pb.c[1, 2:end] + pb.c[1, 1] ≤ u.y[1]) ; push!(ctr_symbol, ctr)
        ctr = JuMP.@constraint(pb.lp_copied, varArray_copied'* pb.c[2, 2:end] + pb.c[2, 1] ≤ u.y[2]) ; push!(ctr_symbol, ctr)
    end

    # todo : LBS bounds which normally no influence 
    JuMP.optimize!(pb.lp_copied, ignore_optimize_hook=true)

    x = JuMP.value.(varArray_copied)
    # retrieve new ctr added 
    for ctr in ctr_symbol
        if JuMP.is_valid(pb.lp_copied, ctr)
            JuMP.delete(pb.lp_copied, ctr) ; JuMP.unregister(pb.lp_copied, :ctr)
        end
    end

end


function feasPumingJumping(node::Node, pb::BO01Problem, incumbent::IncumbentSet, round_results, verbose ; args...)
    LBS = node.RBS.natural_order_vect.sols

    # ∀ l lower bound
    for l in LBS
        if l.is_binaryis_binary continue end 

        H = Vector{Solution}()
        s̄ = rounding_jumping(l, pb)

        push!(H, s̄) ; zone = nadirPtsZone(incumbent, l)
        # todo : be attention if new s̄ is under the zone !

        iter = 0
        while !isFeasible(s̄, pb) && iter < IterLimit
            iter += 1

            #todo : solving lp -> s̃


            if s̃.is_binaryis_binary break end 

            s̄ = rounding_jumping(s̃, pb)

            if contains(H, s̄)
                s̄ = flip(s̃, s̄)
            end
            push!(H, s̄)
        end
    end
    #todo : update incumbent new generated , return 
end