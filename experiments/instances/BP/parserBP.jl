using JuMP, CPLEX

mutable struct BP 
    n :: Int64 
    m :: Int64 
    A :: Matrix{Float64}
    b :: Vector{Float64}
    c :: Vector{Int64}
    objSense :: MOI.OptimizationSense
    name :: String 
end

function BP(n :: Int64, m :: Int64)
    return BP(n , m , zeros(m, n), zeros(m), zeros(Int64, n), MOI.MIN_SENSE, "")
end

"""
Coefficients generator proposed by Pedersen et al, 2008.
"""
function generateC2(c1::Vector{Int64})::Vector{Int64}
    mini = minimum(c1) ; maxi = maximum(c1)
    c2 = zeros(Int64, length(c1))

    for i in 1:length(c1)
        if c1[i] < (maxi-mini)/2
            c2[i] = rand((maxi-c1[i]+mini):(maxi))
        elseif c1[i] >= (maxi-mini)/2
            c2[i] = rand((mini):(maxi-c1[i]+mini))
        else
            error("Coefficient error c1[$i] = $(c1[i]), maximum = $maxi, minimum = $mini !")
        end
    end
    return c2
end



function readModel(fname::String)::BP
    f = open(fname)
    name = split(split(readline(f), " ")[end], ".")[1]
    close(f)

    # loading model ...
    model = read_from_file(fname) ; set_optimizer(model, CPLEX.Optimizer)
    JuMP.set_silent(model)

    # reading model ...
    cstrEqualTo, cstrGreaterThan, cstrLessThan, cstrInterval = [
        JuMP.all_constraints(model, JuMP.GenericAffExpr{Float64,JuMP.VariableRef}, set_type) 
            for set_type in (MOI.EqualTo{Float64}, MOI.GreaterThan{Float64}, MOI.LessThan{Float64}, MOI.Interval{Float64})
    ]
    numRows = sum(length, (cstrEqualTo, cstrGreaterThan, cstrLessThan)) ; numRows += 2*length(cstrInterval)
    varArray = JuMP.all_variables(model); numVars = length(varArray)

    inst = BP(numVars, numRows) ; inst.name=name

    if numVars >= 1000 return inst end 

    optimize!(model) ; solved_time = round(solve_time(model), digits = 2)
    println("mono solved time $(solved_time) \n\n" )

    if solved_time > 300.0 return inst end 

    varIndex = Dict(varArray[i] => i for i=1:length(varArray))

    # objective
    inst.objSense = objective_sense(model) 
    obj = objective_function(model)
    for (var, coeff) in obj.terms
        inst.c[varIndex[var]] = round(Int64, coeff)
    end

    cstr_index = 0
    
    # Ax=b
    for cstrRef in cstrEqualTo
        cstr_index += 1

        con = JuMP.constraint_object(cstrRef)
        rhs = con.set.value

        terms = JuMP.linear_terms(con.func)
        for (coeff, term) in terms
            inst.A[cstr_index, varIndex[term]] = coeff
        end
        inst.b[cstr_index] = rhs
    end

    # Ax≥b
    for cstrRef in cstrGreaterThan
        cstr_index += 1

        con = JuMP.constraint_object(cstrRef)
        rhs = con.set.lower

        terms = JuMP.linear_terms(con.func)
        for (coeff, term) in terms
            inst.A[cstr_index, varIndex[term]] = coeff * -1
        end
        inst.b[cstr_index] = rhs * -1
    end

    # Ax≤b
    for cstrRef in cstrLessThan
        cstr_index += 1

        con = JuMP.constraint_object(cstrRef)
        rhs = con.set.upper

        terms = JuMP.linear_terms(con.func)
        for (coeff, term) in terms
            inst.A[cstr_index, varIndex[term]] = coeff
        end
        inst.b[cstr_index] = rhs
    end
    

    # lb ≤ AX ≤ ub
    for cstrRef in cstrInterval
        cstr_index += 1

        con = JuMP.constraint_object(cstrRef)
        lb = con.set.lower
        ub = con.set.upper

        terms = JuMP.linear_terms(con.func)
        for (coeff, term) in terms
            inst.A[cstr_index, varIndex[term]] = coeff ; inst.A[cstr_index+1, varIndex[term]] = coeff * -1
        end
        inst.b[cstr_index] = ub ; inst.b[cstr_index+1] = lb * -1
        cstr_index += 1
    end

    # # todo write second objective
    # c2 = generateC2(inst.c)
    # folder = "./objective"
    # if !isdir(folder)
    #     mkdir(folder)
    # end

    # outputName = folder * "/" * inst.name
    # fout = open(outputName, "w")
    # println(fout, "c2 = $c2 ")
    # close(fout)
    return inst 
end

# readModel(ARGS[1])

# readModel("./mipLib/cov1075.mps")