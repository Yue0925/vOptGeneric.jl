# MIT License
# Copyright (c) 2017: Xavier Gandibleux, Anthony Przybylski, Gauthier Soleilhac, and contributors.
__precompile__()
module vOptGeneric
import JuMP, MathOptInterface, Combinatorics
import Base.push!
const MOI = MathOptInterface

include("vOptData.jl")

export vModel,
    getvOptData,
    printX_E,
    getY_N,
    getX_E,
    vSolve,
    @addobjective,
    parseMOP,
    writeMOP

include("MOP.jl")
include("algorithms.jl")
include("BO01BB/BO01BB.jl")


function vModel(optimizer_factory; args...)
    m = JuMP.Model(optimizer_factory ; args...)
    m.optimize_hook = vSolve
    # m.printhook = printhook
    m.ext[:vOpt] = vOptData()
    return m
end

function vModel(; args...)
    m = JuMP.Model(; args...)
    m.optimize_hook = vSolve
    # m.printhook = printhook
    m.ext[:vOpt] = vOptData()
    return m
end

function vSolve(m::JuMP.Model ; relax=false, method=nothing, step = 1., round_results = false, verbose = true, args...)

    vd = getvOptData(m)

    if relax != false
        @warn "linear relaxation not yet implemented"
    end
    
    if method == :epsilon
        solve_eps(m, step, round_results, verbose ; relaxation=relax, args...)
    elseif method == :dicho || method == :dichotomy
        solve_dicho(m, round_results, verbose ; relaxation=relax, args...)
    elseif method == :Chalmet || method == :chalmet
        solve_Chalmet(m, step, verbose ; relaxation=relax, args...)
    elseif method == :lex || method == :lexico
        solve_lexico(m, verbose ; relaxation=relax, args...)
    elseif method == :bb || method == :branchbound
        return solve_branchboundcut(m, false, false, round_results, verbose ; args...)
    elseif method == :bb_EPB || method == :branchboundEPB
        return solve_branchboundcut(m, false, true, round_results, verbose ; args...)
    elseif method == :bc || method == :branchcut
        return solve_branchboundcut(m, true, false, round_results, verbose ; args...)
    elseif method == :bc_EPB || method == :branchcutEPB
        return solve_branchboundcut(m, true, true, round_results, verbose ; args...)
    else
        @warn("use solve(m, method = :(epsilon | dichotomy | chalmet | lexico | branchbound) )")
    end

end

# function printhook(io::IO, m::Model)
#     vd = getvOptData(m)
#     for i = 1:length(vd.objs)
#         println(vd.objSenses[i] == :Min ? "Min " : "Max ", vd.objs[i])
#     end
#     str = JuMP.model_str(JuMP.REPLMode, m)
#     index = findfirst("Subject to", str)
#     index !== nothing && print(str[first(index):end])
# end


macro addobjective(m, args...)
    if length(args) != 2
        error("in @addobjective: needs three arguments: model, 
                objective sense (Max or Min) and linear expression.")
    end
    m = esc(m)
    args[1] != :Min && args[1] != :Max && error("in @addobjective: expected Max or Min for objective sense, got $(args[1])")
    sense = args[1] == :Min ? MOI.MIN_SENSE : MOI.MAX_SENSE
    expr = esc(args[2])
    return quote
        f = JuMP.AffExpr() + JuMP.@expression($m, $expr)
        !isa(f, JuMP.GenericAffExpr) && error("in @addobjective : vOptGeneric only supports linear objectives")
        vd = $m.ext[:vOpt]
        push!(vd.objSenses, $(esc(sense)))
        push!(vd.objs, f)
        f
    end
end

function printX_E(m::JuMP.Model)
    vd = getvOptData(m)
    for i = 1:length(vd.Y_N)
        print(vd.Y_N[i]," : ")
        for var in JuMP.all_variables(m)
            val = JuMP.value(var, i)
            if val != 0
                if JuMP.is_binary(var) || JuMP.is_integer(var)
                    print(JuMP.name(var), "=", round(Int, val), " ")
                else
                    print(JuMP.name(var)," =", val," ")
                end
            end
        end 
        println()
   end
end

function JuMP.value(v::JuMP.VariableRef, i::Int)
    vd = getvOptData(v.model)
    return vd.X_E[i][v.index.value]
end

function getY_N(m::JuMP.Model)
    return getvOptData(m).Y_N
end

function getX_E(m::JuMP.Model)
    return getvOptData(m).X_E
end

# function  getLogObjs(m::JuMP.Model)
#     return getvOptData(m).logObjs
# end


end
