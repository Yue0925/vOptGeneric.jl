include("BBtree.jl")
include("../algorithms.jl")
include("struct.jl")
include("separators.jl")
include("cutPool.jl")

using JuMP 
using Random

const IterLimit = 5
const EPS = 0.001

# todo : check all bounds cuts for m/m_copied 

"""
Given a fractional solution `l`, return a rounded jumped interger solution `s`.
"""
function rounding_jumping(l::Solution, pb::BO01Problem; proba::Float64=0.5)::Solution
    x̄ = [-1 for _ in 1:length(l.xEquiv[1])]
    rhs = deepcopy(pb.b)

    # ----------------------------------
    # rounding to 1 with a probability
    # ----------------------------------
    idx = shuffle!(collect(1:length(l.xEquiv[1])))
    for j in idx
        if abs(l.xEquiv[1][j] - 0) ≤ EPS
            x̄[j] = 0
        elseif abs(l.xEquiv[1][j] - 1) ≤ EPS && rand() ≥ proba
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
        idx_sorted = sort!([j for j in 1:length(l.xEquiv[1]) if x̄[j]==-1], by= v ->(preference[v][1], preference[v][2]) )

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
function flip(s̃::Solution, s̄::Solution)::Solution
    n = length(s̄.xEquiv[1])
    T = n/4 ; nb_flipped = rand(round(Int64, T/2):round(Int64, T*3/4))
    x = deepcopy(s̄.xEquiv[1][:])

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
    # ---------------
    col = length(pb.varArray) ; row = length(pb.b)
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

    # todo : (check) LBS bounds which normally no influence 
    JuMP.optimize!(pb.lp_copied, ignore_optimize_hook=true)

    x = JuMP.value.(varArray_copied)
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


function feasPumingJumping(node::Node, pb::BO01Problem, incumbent::IncumbentSet)
    LBS = node.RBS.natural_order_vect.sols
    U_newfea = NaturalOrderVector()

    # ∀ l lower bound
    for l in LBS
        if l.is_binary continue end # push!(U_newfea, l, filtered=true) ; 

        H = Vector{Solution}()
        s̄ = rounding_jumping(l, pb)

        push!(H, s̄) ; zone = nadirPtsZone(incumbent, l)
        # todo : be attention if new s̄ is under the zone !

        iter = 0
        while !isFeasible(s̄, pb) && iter < IterLimit
            iter += 1

            s̃ = Δ_opt(pb, s̄, zone)

            if s̃.is_binaryis_binary break end 

            s̄ = rounding_jumping(s̃, pb)

            if contains(H, s̄)
                s̄ = flip(s̃, s̄)
            end
            push!(H, s̄)
        end
        if isFeasible(s̄, pb) push!(U_newfea, s̄, filtered=true) end
    end

    return U_newfea
end


# -----------------------------
# solving mono MaxCut_10 ... 
# -----------------------------
# n = 55.0
# solved time 0.13

# -----------------------------
# solving MaxCut_10 by bc_rootRelaxCPEPB  ... 
# -----------------------------
# [ Info: 1 new feasible points found !
# ERROR: LoadError: Result index of attribute MathOptInterface.VariablePrimal(1) out of bounds. There are currently 0 solution(s) in the model.
# Stacktrace:
#  [1] check_result_index_bounds
#    @ ~/.julia/packages/MathOptInterface/a4tKm/src/attributes.jl:207 [inlined]
#  [2] get(model::CPLEX.Optimizer, attr::MathOptInterface.VariablePrimal, x::MathOptInterface.VariableIndex)
#    @ CPLEX ~/.julia/packages/CPLEX/X6Mno/src/MOI/MOI_wrapper.jl:2884
#  [3] get(b::MathOptInterface.Bridges.LazyBridgeOptimizer{CPLEX.Optimizer}, attr::MathOptInterface.VariablePrimal, index::MathOptInterface.VariableIndex)
#    @ MathOptInterface.Bridges ~/.julia/packages/MathOptInterface/a4tKm/src/Bridges/bridge_optimizer.jl:1116
#  [4] get(model::MathOptInterface.Utilities.CachingOptimizer{MathOptInterface.Bridges.LazyBridgeOptimizer{CPLEX.Optimizer}, MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}}, attr::MathOptInterface.VariablePrimal, index::MathOptInterface.VariableIndex)
#    @ MathOptInterface.Utilities ~/.julia/packages/MathOptInterface/a4tKm/src/Utilities/cachingoptimizer.jl:911
#  [5] _moi_get_result(::MathOptInterface.Utilities.CachingOptimizer{MathOptInterface.Bridges.LazyBridgeOptimizer{CPLEX.Optimizer}, MathOptInterface.Utilities.UniversalFallback{MathOptInterface.Utilities.Model{Float64}}}, ::MathOptInterface.VariablePrimal, ::Vararg{Any})
#    @ JuMP ~/.julia/packages/JuMP/0STkJ/src/JuMP.jl:1134
#  [6] get(model::Model, attr::MathOptInterface.VariablePrimal, v::VariableRef)
#    @ JuMP ~/.julia/packages/JuMP/0STkJ/src/JuMP.jl:1191
#  [7] value(v::VariableRef; result::Int64)
#    @ JuMP ~/.julia/packages/JuMP/0STkJ/src/variables.jl:1005
#  [8] value
#    @ ~/.julia/packages/JuMP/0STkJ/src/variables.jl:1004 [inlined]
#  [9] _broadcast_getindex_evalf
#    @ ./broadcast.jl:670 [inlined]
# [10] _broadcast_getindex
#    @ ./broadcast.jl:643 [inlined]
# [11] getindex
#    @ ./broadcast.jl:597 [inlined]
# [12] macro expansion
#    @ ./broadcast.jl:961 [inlined]
# [13] macro expansion
#    @ ./simdloop.jl:77 [inlined]
# [14] copyto!
#    @ ./broadcast.jl:960 [inlined]
# [15] copyto!
#    @ ./broadcast.jl:913 [inlined]
# [16] copy
#    @ ./broadcast.jl:885 [inlined]
# [17] materialize
#    @ ./broadcast.jl:860 [inlined]
# [18] Δ_opt(pb::Main.vOptGeneric.BO01Problem, s::Main.vOptGeneric.Solution, nadir_pts::Main.vOptGeneric.NaturalOrderVector)
#    @ Main.vOptGeneric ~/coding/vOptGeneric.jl/src/BO01BB/metaHeuristic.jl:164
# [19] feasPumingJumping(node::Main.vOptGeneric.Node, pb::Main.vOptGeneric.BO01Problem, incumbent::Main.vOptGeneric.IncumbentSet)
#    @ Main.vOptGeneric ~/coding/vOptGeneric.jl/src/BO01BB/metaHeuristic.jl:199
# [20] compute_LBS(node::Main.vOptGeneric.Node, pb::Main.vOptGeneric.BO01Problem, incumbent::Main.vOptGeneric.IncumbentSet, round_results::Bool, verbose::Bool; args::Base.Pairs{Symbol, Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}}, Tuple{Symbol}, NamedTuple{(:args,), Tuple{Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}}}}})
#    @ Main.vOptGeneric ~/coding/vOptGeneric.jl/src/BO01BB/cuttingPlanes.jl:58
# [21] LPRelaxByDicho(node::Main.vOptGeneric.Node, pb::Main.vOptGeneric.BO01Problem, incumbent::Main.vOptGeneric.IncumbentSet, round_results::Bool, verbose::Bool; args::Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}})
#    @ Main.vOptGeneric ~/coding/vOptGeneric.jl/src/BO01BB/fathoming.jl:82
# [22] LPRelaxByDicho
#    @ ~/coding/vOptGeneric.jl/src/BO01BB/fathoming.jl:76 [inlined]
# [23] macro expansion
#    @ ~/.julia/packages/TimerOutputs/4yHI4/src/TimerOutput.jl:237 [inlined]
# [24] iterative_procedure(todo::DataStructures.Queue{Base.RefValue{Main.vOptGeneric.Node}}, node::Main.vOptGeneric.Node, pb::Main.vOptGeneric.BO01Problem, incumbent::Main.vOptGeneric.IncumbentSet, worst_nadir_pt::Vector{Float64}, round_results::Bool, verbose::Bool; args::Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}})
#    @ Main.vOptGeneric ~/coding/vOptGeneric.jl/src/BO01BB/BO01BB.jl:221
# [25] iterative_procedure(todo::DataStructures.Queue{Base.RefValue{Main.vOptGeneric.Node}}, node::Main.vOptGeneric.Node, pb::Main.vOptGeneric.BO01Problem, incumbent::Main.vOptGeneric.IncumbentSet, worst_nadir_pt::Vector{Float64}, round_results::Bool, verbose::Bool)
#    @ Main.vOptGeneric ~/coding/vOptGeneric.jl/src/BO01BB/BO01BB.jl:136
# [26] macro expansion
#    @ ~/.julia/packages/TimerOutputs/4yHI4/src/TimerOutput.jl:237 [inlined]
# [27] macro expansion
#    @ ~/coding/vOptGeneric.jl/src/BO01BB/BO01BB.jl:373 [inlined]
# [28] macro expansion
#    @ ~/.julia/packages/TimerOutputs/4yHI4/src/TimerOutput.jl:237 [inlined]
# [29] solve_branchboundcut(m::Model, cp::Bool, root_relax::Bool, EPB::Bool, round_results::Bool, verbose::Bool; args::Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}})
#    @ Main.vOptGeneric ~/coding/vOptGeneric.jl/src/BO01BB/BO01BB.jl:366
# [30] solve_branchboundcut
#    @ ~/coding/vOptGeneric.jl/src/BO01BB/BO01BB.jl:305 [inlined]
# [31] vSolve(m::Model; relax::Bool, method::Symbol, step::Float64, round_results::Bool, verbose::Bool, args::Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}})
#    @ Main.vOptGeneric ~/coding/vOptGeneric.jl/src/vOptGeneric.jl:70
# [32] vopt_solve(method::Symbol, outputName::String; step::Float64)
#    @ Main ~/coding/vOptGeneric.jl/experiments/instances/MaxCutalea/vOptMaxCut.jl:79
# [33] vopt_solve(method::Symbol, outputName::String)
#    @ Main ~/coding/vOptGeneric.jl/experiments/instances/MaxCutalea/vOptMaxCut.jl:33
# [34] solve(fname::String, method::String)
#    @ Main ~/coding/vOptGeneric.jl/experiments/instances/MaxCutalea/vOptMaxCut.jl:143
# [35] top-level scope
#    @ ~/coding/vOptGeneric.jl/experiments/instances/MaxCutalea/vOptMaxCut.jl:148
# in expression starting at /home/yue/coding/vOptGeneric.jl/experiments/instances/MaxCutalea/vOptMaxCut.jl:148
# ./instances/MaxCut_10 ...  bc_rootRelaxCP
# ^Z
# [10]+  Stopped                 ./run.sh
