include("parserHardKP.jl")

using JuMP, CPLEX

# include("../../../src/vOptGeneric.jl")
# using .vOptGeneric 


# function vopt_solve(inst::HKP, method, fname, outputName; step=0.5)
#     # ---- setting the model
#     model = vModel( CPLEX.Optimizer )
#     JuMP.set_silent(model)

#     @variable(model, x[1:inst.n], Bin)
#     @constraint(model, x'* inst.A ≤ inst.b)
#     @addobjective(model, Max, x'* inst.c)
#     #todo : construct C2 
#     # @addobjective(model, Min, sum(inst.C2[j] * x[j] for j=1:inst.n))


#     if method == :bb
#         infos = vSolve( model, method=:bb, verbose=false )
#         println(infos)
#     elseif method == :bc 
#         infos = vSolve( model, method=:bc, verbose=false )
#         println(infos)
#     elseif method == :dicho 
#         start = time()
#         vSolve( model, method=:dicho, verbose=false )
#         total_time = round(time() - start, digits = 2)
#     elseif method==:epsilon 
#         start = time()
#         vSolve( model, method=:epsilon, step=step, verbose=false )
#         total_time = round(time() - start, digits = 2)
#         println("total_time = $total_time ")

#     elseif method == :bc_rootRelax 
#         infos = vSolve( model, method=:bc_rootRelax, verbose=false )
#         println(infos)
#     elseif method == :bc_rootRelaxCP  
#         infos = vSolve( model, method=:bc_rootRelaxCP , verbose=false )
#         println(infos)    
#     elseif method == :bb_EPB 
#         infos = vSolve( model, method=:bb_EPB, verbose=false )
#         println(infos)
#     elseif method == :bc_rootRelaxEPB
#         infos = vSolve( model, method=:bc_rootRelaxEPB, verbose=false )
#         println(infos)
#     elseif method == :bc_rootRelaxCPEPB
#         infos = vSolve( model, method=:bc_rootRelaxCPEPB, verbose=false )
#         println(infos)
#     elseif method == :bc_EPB
#         infos = vSolve( model, method=:bc_EPB, verbose=false )
#         println(infos)
#     else
#         @error "Unknown method parameter $(method) !"
#     end

#     # ---- Querying the results
#     Y_N = getY_N( model )
#     println("length Y_N = ", length(Y_N))

#     X_E = getX_E( model )
#     # println("length X_E = ", length(X_E))


#     # (method != :dicho && method != :epsilon) ? writeResults(inst.n, inst.m, fname, outputName, method, Y_N, X_E; infos) :
#     #     writeResults(inst.n, inst.m, fname, outputName, method, Y_N, X_E; total_time)

# end

function solve(fname::String)
    instances = readHKP(fname) 

    for inst in instances
        if inst.n >= 1000 continue end 
        println("\n -----------------------------")
        println(" building $(inst.name) ... ")
        println(" -----------------------------")
 
        model = Model(CPLEX.Optimizer) ; JuMP.set_silent(model)
        @variable(model, x[1:inst.n], Bin )
        @objective(model, Max, x'* inst.c)

        @constraint(model, x'* inst.A ≤ inst.b)

        # optimize
        optimize!(model) 
        println(" n = $(inst.n) , m = 1")
        println("solved time $(round(solve_time(model), digits = 2))" )
    end
end

solve(ARGS[1])
