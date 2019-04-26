# MIT License
# Copyright (c) 2017: Xavier Gandibleux, Anthony Przybylski, Gauthier Soleilhac, and contributors.
function solve_lexico(m::Model, verbose; args...)
    #Retrieve objectives and their senses from vOptData
    vd = getvOptData(m)
    empty!(vd.Y_N) ; empty!(vd.X_E)
    objs = vd.objs
    objSenses = vd.objSenses
    nbObj = length(objs)

    #Check that problem is feasible and that no objective is unbounded
    for i = 1:nbObj
        JuMP.set_objective(m, objSenses[i], objs[i])
        @suppress optimize!(m, ignore_optimize_hook=true )
        # status != :Optimal && return status TODO
    end

    #Create the constraints on the objectives
    #JuMP doesn't support changing constraint coefficients so we use >= and <= constraints,
    #and typemin/typemax to deactivate the constraints
    cstr_obj = ScalarConstraint[]
    for i = 1:nbObj
        if objSenses[i] == JuMP.MOI.MAX_SENSE
            cstr = @build_constraint(objs[i].aff >= -1e16)
            push!(cstr_obj, cstr)
            JuMP.add_constraint(m, cstr)
        else
            cstr = @build_constraint(objs[i].aff <= -1e16)
            push!(cstr_obj, cstr)
            JuMP.add_constraint(m, cstr)
        end
    end

    for p in permutations(1:nbObj, nbObj)
        verbose && println("solving for objectives $p")
        solve_permutation(m, p, cstr_obj ; args...)
    end
   
    return :Optimal
end

function solve_permutation(m::Model, p, cstr_obj ; args...)

    vd = getvOptData(m)
    objs = vd.objs
    objSenses = vd.objSenses

    #Set the first objective of the permutation as an objective in the JuMP model
    JuMP.set_objective(m, objSenses[p[1]], objs[p[1]])

    #Solve with that objective
    status = @suppress optimize!(m, ignore_optimize_hook=true)
    status != :Optimal && return status

    for i = 2:length(p)
        fVal = m.objVal #get the value for the last objective solved
        slack = objSenses[p[i-1]] == JuMP.MOI.MAX_SENSE ? -1e-6 : 1e-6
        JuMP.setRHS(cstr_obj[p[i-1]], fVal - objs[p[i-1]].aff.constant + slack) #set the constraint for the last objective solved
        JuMP.set_objective(m, objSenses[p[i]], objs[p[i]]) #set the i-th objective of the permutation in the JuMP model
        @suppress optimize!(m, ignore_optimize_hook=true ) #and solve
    end    

    varArray = JuMP.all_variables(m)
    #Store results in vOptData
    push!(vd.Y_N, map(JuMP.getvalue, objs))
    push!(vd.X_E, JuMP.getvalue.(varArray))

    for i = 1:length(p)
        if objSenses[i] == JuMP.MOI.MAX_SENSE
            JuMP.setRHS(cstr_obj[i], -1e16)
        else
            JuMP.setRHS(cstr_obj[i], 1e16)
        end
    end

    nothing
end

function solve_eps(m::Model, ϵ::Float64, round_results, verbose ; args...)
    #Retrieve objectives and their senses from vOptData
    vd = getvOptData(m)
    empty!(vd.Y_N) ; empty!(vd.X_E)
    f1,f2 = vd.objs[1],vd.objs[2]
    f1Sense,f2Sense = vd.objSenses[1],vd.objSenses[2]
    
    #Set the first objective as an objective in the JuMP model
    JuMP.set_objective(m, f1Sense, f1)
    
    #Solve with that objective
    optimize!(m, ignore_optimize_hook=true)
    status = termination_status(m)
    @show status
    #If a solution exists
    if status == JuMP.MOI.OPTIMAL

        varArray = JuMP.all_variables(m)

        #Get the score on the objectives
        f1Val = value(f1)
        f2Val = value(f2)

        #Store results in vOptData
        push!(vd.X_E, value.(varArray))
        push!(vd.Y_N, [f1Val, f2Val])

        R1 = f1Sense==JuMP.MOI.MIN_SENSE ? (<=) : (>=)
        R2 = f2Sense==JuMP.MOI.MIN_SENSE ? (<=) : (>=)
        weak_dom(a, b) = R1(a[1], b[1]) && R2(a[2], b[2])

        #Declare the epsilon-constraint (RHS will be set later)
        RHS = @variable(m)
        eps = f2Sense == JuMP.MOI.MIN_SENSE ? @constraint(m, f2 <= RHS) : @constraint(m, f2 >= RHS)
        first_optimize = true

        #While a solution exists
        while status == JuMP.MOI.OPTIMAL

            if !first_optimize
                Y_m1 = vd.Y_N[end]
                #If last solution found is dominated by this one
                if weak_dom((f1Val, f2Val), Y_m1)
                    #Remove last solution from Y_N and X_E
                    pop!(vd.Y_N) ; pop!(vd.X_E)
                end
                #Store results in vOptData
                push!(vd.X_E, value.(varArray))
                push!(vd.Y_N, [f1Val, f2Val])
            end

            verbose && print("z1 = ", f1Val, ", z2 = ", f2Val)

            #Set the RHS of the epsilon-constraint
            if f2Sense == JuMP.MOI.MIN_SENSE
                JuMP.fix(RHS, f2Val - f2.constant - ϵ)
                verbose && println(". Solving with f2 <= ", f2Val - ϵ)
            else
                JuMP.fix(RHS, f2Val - f2.constant + ϵ)
                verbose && println(". Solving with f2 >= ", f2Val + ϵ)
            end

            #And solve again
            first_optimize = false
            optimize!(m, ignore_optimize_hook=true)
            status = termination_status(m)
            #Get the score on the objectives
            f1Val = value(f1)
            f2Val = value(f2)

        end

        #Sort X_E and Y_N
        s = sortperm(vd.Y_N, by = e -> e[1])
        vd.Y_N = vd.Y_N[s]
        vd.X_E = vd.X_E[s]

        #Remove the epsilon constraint
        delete(m, eps) ; delete(m, RHS)
    else
        return status
    end
    return :Optimal
end

function solve_dicho(m::Model, round_results ; args...)
    vd = getvOptData(m)
    empty!(vd.Y_N) ; empty!(vd.X_E)
    f1,f2 = vd.objs[1],vd.objs[2]
    f1Sense,f2Sense = vd.objSenses[1],vd.objSenses[2]
    varArray = JuMP.all_variables(m)

    #Set the first objective as an objective in the JuMP model
    JuMP.set_objective(m, f1Sense, f1)
    
    #Solve with that objective
    status = optimize!(m, ignore_optimize_hook=true ;  args...)

    #If a solution exists
    if status == :Optimal

        yr_1 = m.objVal
        yr_2 = JuMP.getvalue(f2)

        #Store results in vOptData
        push!(vd.Y_N, round_results ? round.([yr_1, yr_2]) : [yr_1, yr_2])
        push!(vd.X_E, JuMP.getvalue.(varArray))


        #Set the second objective as an objective in the JuMP model
        JuMP.set_objective(m, f2sense, f2)

        #Solve with that objective
        status = optimize!(m, ignore_optimize_hook=true )

        if status == :Optimal

            ys_1 = JuMP.getvalue(f1)
            ys_2 = m.objVal

            if !isapprox(yr_1, ys_1, atol=1e-3) || !isapprox(yr_2, ys_2, atol=1e-3)
                push!(vd.Y_N, round_results ? round.([ys_1, ys_2]) : [ys_1, ys_2])
                push!(vd.X_E, JuMP.getvalue.(varArray))
                dichoRecursion(m, yr_1, yr_2, ys_1, ys_2, varArray, round_results ; args...)
            end
        
            #Sort X_E and Y_N
            s = sortperm(vd.Y_N, by = first)
            vd.Y_N = vd.Y_N[s]
            vd.X_E = vd.X_E[s]

            R1 = f1Sense==JuMP.MOI.MIN_SENSE ? (<=) : (>=)
            R2 = f2Sense==JuMP.MOI.MIN_SENSE ? (<=) : (>=)
            weak_dom(a, b) = R1(a[1], b[1]) && R2(a[2], b[2])

            #Filter X_E and Y_N :
            inds = Int[]
            for i = 1:length(vd.Y_N)-1
                if weak_dom(vd.Y_N[i], vd.Y_N[i+1])
                    push!(inds, i+1)
                elseif weak_dom(vd.Y_N[i+1], vd.Y_N[i])
                    push!(inds, i)
                end
            end
            deleteat!(vd.Y_N, inds)
            deleteat!(vd.X_E, inds)
        end
    end

    status
end

function dichoRecursion(m::Model, yr_1, yr_2, ys_1, ys_2, varArray, round_results ; args...)

    vd = getvOptData(m)
    f1,f2 = vd.objs[1],vd.objs[2]
    f1Sense,f2Sense = vd.objSenses[1],vd.objSenses[2]

    λ1 = abs(yr_2 - ys_2)
    λ2 = abs(ys_1 - yr_1)

    m.objSense = f1Sense

    if f1Sense==f2Sense
        lb = λ1*yr_1 + λ2*yr_2
        JuMP.set_objective_function(λ1*f1 + λ2*f2)
    else
        lb = λ1*yr_1 - λ2*yr_2
        JuMP.set_objective_function(λ1*f1 - λ2*f2)
    end

    optimize!(m, ignore_optimize_hook=true )

    yt_1 = JuMP.getvalue(f1)
    yt_2 = JuMP.getvalue(f2)


    if f1Sense==f2Sense
        val = λ1*yt_1 + λ2*yt_2
    else
        val = λ1*yt_1 - λ2*yt_2
    end

    if (f1Sense == JuMP.MOI.MIN_SENSE && val < lb - 1e-4) || val > lb + 1e-4
        push!(vd.Y_N, round_results ? round.([yt_1, yt_2]) : [yt_1, yt_2])
        push!(vd.X_E, JuMP.getvalue.(varArray))
        dichoRecursion(m, yr_1, yr_2, yt_1, yt_2, varArray, round_results ; args...)
        dichoRecursion(m, yt_1, yt_2, ys_1, ys_2, varArray, round_results ; args...)
    end

end


function solve_Chalmet(m::Model, step ; args...)
    vd = getvOptData(m)
    empty!(vd.Y_N) ; empty!(vd.X_E)
    f1,f2 = vd.objs[1],vd.objs[2]
    f1Sense,f2Sense = vd.objSenses[1],vd.objSenses[2]
    varArray = JuMP.all_variables(m)

    #Set the first objective as an objective in the JuMP model
    JuMP.set_objective(m, f1Sense, f1)

    #Solve with that objective
    status = optimize!(m, ignore_optimize_hook=true )

    #If a solution exists
    if status == :Optimal

        yr_1 = m.objVal
        yr_2 = JuMP.getvalue(f2)

        #Store results in vOptData
        push!(vd.Y_N, [yr_1, yr_2])
        push!(vd.X_E, JuMP.getvalue.(varArray))


        #Set the second objective as an objective in the JuMP model
        JuMP.set_objective(m, f2Sense, f2)


        #Solve with that objective
        status = optimize!(m, ignore_optimize_hook=true )

        if status == :Optimal

            ys_1 = JuMP.getvalue(f1)
            ys_2 = m.objVal

            if !isapprox(yr_1, ys_1, atol=1e-3) || !isapprox(yr_2, ys_2, atol=1e-3)
                push!(vd.Y_N, [ys_1, ys_2])
                push!(vd.X_E, JuMP.getvalue.(varArray))
                #Declare the constraints on z1 and z2 (RHS will be set later)
                if f1Sense == JuMP.MOI.MIN_SENSE
                    cstrz1 = @constraint(m, f1.aff <= 0.0)
                else
                    cstrz1 = @constraint(m, f1.aff >= 0.0)
                end
                if f2Sense == JuMP.MOI.MIN_SENSE
                    cstrz2 = @constraint(m, f2.aff <= 0.0)
                else
                    cstrz2 = @constraint(m, f2.aff >= 0.0)
                end
                ChalmetRecursion(m, yr_1, yr_2, ys_1, ys_2, varArray, cstrz1, cstrz2, step ; args...)

                #Sort X_E and Y_N
                s = sortperm(vd.Y_N, by = e -> e[1])
                vd.Y_N = vd.Y_N[s]
                vd.X_E = vd.X_E[s]

                #We can use weak dominance definition here
                dom_min_min(a,b) = a[1] <= b[1] && a[2] <= b[2]
                dom_max_max(a,b) = a[1] >= b[1] && a[2] >= b[2]
                dom_min_max(a,b) = a[1] <= b[1] && a[2] >= b[2]
                dom_max_min(a,b) = a[1] >= b[1] && a[2] <= b[2]

                if f1Sense==f2Sense==JuMP.MOI.MIN_SENSE
                    dominates = dom_min_min
                elseif f1Sense==f2Sense==JuMP.MOI.MAX_SENSE
                    dominates = dom_max_max
                elseif f1Sense==JuMP.MOI.MIN_SENSE
                    dominates = dom_min_max
                else
                    dominates = dom_max_min
                end

                #Filter X_E and Y_N :
                inds = Int[]
                for i = 1:length(vd.Y_N)-1
                    if dominates(vd.Y_N[i], vd.Y_N[i+1])
                        push!(inds, i+1)
                    elseif dominates(vd.Y_N[i+1], vd.Y_N[i])
                        push!(inds, i)
                    end
                end
                deleteat!(vd.Y_N, inds)
                deleteat!(vd.X_E, inds)


                if f1Sense == JuMP.MOI.MIN_SENSE
                    JuMP.setRHS(cstrz1, 1e16)
                else
                    JuMP.setRHS(cstrz1,-1e16)
                end

                if f2Sense == JuMP.MOI.MIN_SENSE
                    JuMP.setRHS(cstrz2, 1e16)
                else
                    JuMP.setRHS(cstrz2,-1e16)
                end
            end
        end
    end

    return status
end

function ChalmetRecursion(m::Model, yr_1, yr_2, ys_1, ys_2, varArray, cstrz1, cstrz2, ϵ ; args...)

    vd = getvOptData(m)
    f1,f2 = vd.objs[1],vd.objs[2]
    f1Sense,f2Sense = vd.objSenses[1],vd.objSenses[2]

    m.objSense = f1Sense
    if f1Sense==f2Sense
        JuMP.set_objective_function(f1 + f2)
    else
        JuMP.set_objective_function(f1 - f2)
    end

    lbz1 = f1Sense==JuMP.MOI.MAX_SENSE ? min(yr_1, ys_1) : max(yr_1, ys_1)
    lbz2 = f2Sense==JuMP.MOI.MAX_SENSE ? min(yr_2, ys_2) : max(yr_2, ys_2)

    if f1Sense == JuMP.MOI.MIN_SENSE
        JuMP.setRHS(cstrz1, lbz1 - f1.aff.constant - ϵ)
        print(". Solving with f1 <= ", lbz1 - ϵ)
    else
        JuMP.setRHS(cstrz1, lbz1 - f1.aff.constant + ϵ)
        print(". Solving with f1 >= ", lbz1 + ϵ)
    end

    if f2Sense == JuMP.MOI.MIN_SENSE
        JuMP.setRHS(cstrz2, lbz2 - f2.aff.constant - ϵ)
        println(" and f2 <= ", lbz2 - ϵ)
    else
        JuMP.setRHS(cstrz2, lbz2 - f2.aff.constant + ϵ)
        println(" and f2 >= ", lbz2 + ϵ)
    end

    status = @suppress optimize!(m, ignore_optimize_hook=true)

    if status == :Optimal

        yt_1 = JuMP.getvalue(f1)
        yt_2 = JuMP.getvalue(f2)

        push!(vd.Y_N, [yt_1, yt_2])
        push!(vd.X_E, JuMP.getvalue.(varArray))
        ChalmetRecursion(m, yr_1, yr_2, yt_1, yt_2, varArray, cstrz1, cstrz2, ϵ ; args...)
        ChalmetRecursion(m, yt_1, yt_2, ys_1, ys_2, varArray, cstrz1, cstrz2, ϵ ; args...)

    end
end