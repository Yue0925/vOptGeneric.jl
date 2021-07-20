# Bi-objective linear assignment problem (bilap)
#
# Example 9.38 (from Ulungu and Teghem, 1994), page 255 of
# Multicriteria Optimization (2nd edt), M. Ehrgott, Springer 2005.


# ---- Packages to use
using vOptGeneric, JuMP, GLPK


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


# ---- Invoking the solver (Chalmet method)
vSolve( bilap, method=:chalmet, step=0.5 )


# ---- Querying the results
Y_N = getY_N( bilap )


# ---- Displaying the results
printX_E( bilap )
