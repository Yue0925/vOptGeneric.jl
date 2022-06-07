# Bi-objective linear assignment problem (bilap)
#
# Example 9.38 (from Ulungu and Teghem, 1994), page 255 of
# Multicriteria Optimization (2nd edt), M. Ehrgott, Springer 2005.


# ---- Packages to use
using JuMP, GLPK

include("../../../src/vOptGeneric.jl")
include("../../../src/BO01BB/displayGraphic.jl")
using .vOptGeneric


function writeResults(vars::Int64, constr::Int64, fname::String, outputName::String, method, Y_N; total_time=nothing, infos=nothing)

        fout = open(outputName, "w")
        println(fout, "vars = $vars ; constr = $constr ")
      
        if method == :bb
          println(fout, infos)
        else
          println(fout, "total_times_used = $total_time")
        end
        println(fout, "size_Y_N = ", length(Y_N))
        println(fout, "Y_N = ", Y_N)
      
        close(fout)
      
        displayGraphics(fname,Y_N, outputName)
end

function vSolveBOLAP(method::Symbol, fname::String; step=0.5)

        # ---- Values of the instance to solve
        C1 = [  5  1  4  7 ;  # coefficients's vector of the objective 1
        6  2  2  6 ;
        2  8  4  4 ;
        3  5  7  1   ]

        C2 = [  3  6  4  2 ;  # coefficients's vector of the objective 2
        1  3  8  3 ;
        5  2  2  3 ;
        4  2  3  5   ]

        n  = size(C2,1)       # number of lines/columns


        # ---- setting the model
        bilap = vModel( GLPK.Optimizer )

        @variable( bilap, x[1:n,1:n], Bin )

        @addobjective( bilap , Min, sum( C1[i,j]*x[i,j] for i=1:n,j=1:n ) )
        @addobjective( bilap , Min, sum( C2[i,j]*x[i,j] for i=1:n,j=1:n ) )

        @constraint( bilap , cols[i=1:n], sum(x[i,j] for j=1:n) == 1 )
        @constraint( bilap , rows[j=1:n], sum(x[i,j] for i=1:n) == 1 )

        if method == :bb
                infos = vSolve( bilap, method=:bb, verbose=false )
        elseif method == :dicho 
                start = time()
                vSolve( bilap, method=:dicho, verbose=false )
                total_time = round(time() - start, digits = 2)
        elseif method==:epsilon 
                start = time()
                vSolve( bilap, method=:epsilon, step=step, verbose=false )
                total_time = round(time() - start, digits = 2)
        else
                @error "Unknown method parameter $(method) !"
        end

        # ---- Querying the results
        Y_N = getY_N( bilap )

        (method == :bb) ? 
                writeResults(n, 2*n, "LinearAssignmentUlungu1994", fname, method, Y_N; infos) : 
                writeResults(n, 2*n, "LinearAssignmentUlungu1994", fname, method, Y_N; total_time)

end


function main()
        folder = "../../results/smallExamples/"
        for method in [ :bb] # :dicho, :epsilon,
                result_dir = method≠:bb ? folder * "/" * string(method) : folder * "/" * string(method) * "/default"
                if !isdir(result_dir)
                        mkdir(result_dir)
                end
                fname = result_dir * "/" * "LinearAssignmentUlungu1994"

                vSolveBOLAP(method, fname)
        end
end


function run_epsilon_ctr(epsilon::String)
        step = parse(Float64, epsilon)
        folder = "../../results/smallExamples/"
        method = :epsilon
        result_dir = folder * "/" * string(method) * "/" * string(method) * "_" * string(step)
        if !isdir(result_dir)
                mkdir(result_dir)
        end
        fname = result_dir * "/" * "LinearAssignmentUlungu1994"

        vSolveBOLAP(method, fname; step=step)  
end

# run_epsilon_ctr(ARGS[1])

main()