# ==============================================================================
# Xavier Gandibleux - November 2021
#   Implemented in Julia 1.6

# ==============================================================================
using PyPlot

# ------------------------------------------------------------------------------
# kung Algorithm for extracting SN from a static set S of points
# function kungAlgorithm(S1, S2)
#     S = []
#     for i=1:length(S1)
#         push!(S, (S1[i] , S2[i]) )
#     end
#     sort!(S, by = x -> x[1])
#     SN=[] ; push!(SN, S[1]) ; minS2 = S[1][2]
#     for i=2:length(S1)
#         if S[i][2] < minS2
#             push!(SN, S[i]) ; minS2 = S[i][2]
#         end
#     end
#     return SN
# end

# ------------------------------------------------------------------------------
# compute the corner points of a set (S1,S2) and its lower envelop
function computeCornerPointsLowerEnvelop(S1, S2)
    # when all objectives have to be maximized
    Env1=[]; Env2=[]
    for i in 1:length(S1)-1
        push!(Env1, S1[i]); push!(Env2, S2[i])
        push!(Env1, S1[i+1]); push!(Env2, S2[i])
    end
    push!(Env1, S1[end]);push!(Env2, S2[end])
    return Env1,Env2
end

# ------------------------------------------------------------------------------
# display different results
function displayGraphics(fname,YN, output::String, LBS, λ)
    println("displayGraphics")
    PlotOrthonormedAxis = true  # Axis orthonormed or not
    DisplayYN   = true          # Non-dominated points corresponding to efficient solutions
    DisplayUBS  = false         # Points belonging to the Upper Bound Set
    DisplayLBS  = true         # Points belonging to the Lower Bound Set
    DisplayInt  = false         # Points corresponding to integer solutions
    DisplayProj = false         # Points corresponding to projected solutions
    DisplayFea  = false         # Points corresponding to feasible solutions
    DisplayPer  = false         # Points corresponding to perturbated solutions

    YN_1=[];YN_2=[]
    for i in 1:length(YN)
        push!(YN_2, YN[i][2])
        push!(YN_1, YN[i][1])
    end

    xL = []; yL = []
    for i in 1:length(LBS)
        push!(xL, LBS[i][1])
        push!(yL, LBS[i][2])
    end

    ixL = []; iyL = []
    if length(LBS) == 3
        push!(ixL, LBS[1][1]) ; push!(iyL, LBS[1][2])
        k = (λ[2][2]/λ[2][1]) ; b = LBS[2][1] + k * LBS[2][2]
        push!(ixL, (b-k*LBS[1][2])) ; push!(iyL, LBS[1][2])
        push!(ixL, LBS[3][1]) ; push!(iyL, (b-LBS[3][1])/k)
        push!(ixL, LBS[3][1]) ; push!(iyL, LBS[3][2])

    else
        for i in 1:length(LBS)
            push!(ixL, LBS[i][1])
            push!(iyL, LBS[i][2])
        end
    end

    # --------------------------------------------------------------------------
    # Setup
    figure("Project MOMH 2021") # ,figsize=(6.5,5) 
    # if PlotOrthonormedAxis
    #     vmin = 0.99 * min(minimum(YN_1),minimum(YN_2))
    #     vmax = 1.01 * max(maximum(YN_1),maximum(YN_2))
    #     xlim(vmin,vmax)
    #     ylim(vmin,vmax)
    # end
    xlabel(L"z^2(x)")
    ylabel(L"z^1(x)")
    PyPlot.title("Bi-01BKP | $fname")

    # --------------------------------------------------------------------------
    # Display Non-Dominated points
    if DisplayYN
        # display only the points corresponding to non-dominated points
        scatter(YN_2, YN_1, color="red", marker="x", label = L"y \in UBS")
        # display segments joining adjacent non-dominated points
        # plot(YN_2, YN_1, color="red", linewidth=0.75, marker="+", markersize=1.0, linestyle=":")
        # display segments joining non-dominated points and their corners points
        Env1,Env2 = computeCornerPointsLowerEnvelop(YN_2, YN_1)
        plot(Env1,Env2, color="red", linewidth=0.75, marker="x", markersize=1.0, linestyle=":")
    end

    # --------------------------------------------------------------------------
    # Display a Lower bound set (dual, by excess)
    if DisplayLBS
        plot(iyL, ixL, color="green", linewidth=0.75, linestyle=":")
        scatter(yL, xL, color="green", marker="x", label = L"y \in LBS")
    end

    # --------------------------------------------------------------------------
    # Display integer points (feasible and non-feasible in GravityMachine)
    if DisplayInt
        scatter(XInt,YInt, color="orange", marker="s", label = L"y"*" rounded")
    end

    # --------------------------------------------------------------------------
    # Display projected points (points Δ(x,x̃) in GravityMachine)
    if DisplayProj
        scatter(XProj,YProj, color="red", marker="x", label = L"y"*" projected")
    end

    # --------------------------------------------------------------------------
    # Display feasible points
    if DisplayFea
        scatter(XFeas,YFeas, color="green", marker="o", label = L"y \in F")
    end

    # --------------------------------------------------------------------------
    # Display perturbed points (after a cycle in GravityMachine)
    if DisplayPer
        scatter(XPert,YPert, color="magenta", marker="s", label ="pertub")
    end

    # --------------------------------------------------------------------------
    legend(bbox_to_anchor=[1,1], loc=0, borderaxespad=0, fontsize = "x-small")
    savefig(output * ".png")
    PyPlot.close()
end

# ==============================================================================

LBS = [ [-108470.0, -288290.0] , [-111490.0, -269114.0] , [-112846.0, -207674.0] , ] 
λ = [ [0.0, 1.0] , [80616.0, 4376.0] , [1.0, 0.0] , ] 

UBS = [ [-84825.0, -332671.0] , [-86960.0, -329588.0] , [-87523.0, -326528.0] , [-87565.0, -324666.0] , [-88510.0, -322731.0] , [-89193.0, -318703.0] , [-89913.0, -318003.0] , [-89965.0, -317897.0] , [-90010.0, -317086.0] , [-94793.0, -315399.0] , [-96053.0, -310686.0] , [-107650.0, -309761.0] , [-108075.0, -303050.0] , [-108550.0, -302600.0] , [-109270.0, -301900.0] , [-110700.0, -297259.0] , [-111090.0, -290685.0] , [-115048.0, -285634.0] , [-115643.0, -278995.0] , [-118852.0, -275280.0] , [-120392.0, -260067.0] , [-121112.0, -254598.0] , [-121875.0, -227865.0] , [-122582.0, -212225.0] , [-124813.0, -191151.0] , [-125824.0, -183274.0] , [-128520.0, -168403.0] , [-129592.0, -162196.0] , [-129921.0, -139262.0] , [-131758.0, -135393.0] , [-132183.0, -128682.0] , ]

displayGraphics("node_1144", UBS, "node_1144", LBS, λ)
