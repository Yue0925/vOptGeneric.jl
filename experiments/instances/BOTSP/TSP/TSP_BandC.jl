using JuMP
using CPLEX
using CPUTime
using Graphs

include("TSP_IO.jl")
include("TSP_BandC_SeparationAlgo.jl")

# MOI is a shortcut for MathematicalOptimizationInterface

function edge_value(x,i,j)
   if (i<j)
      return value(x[i,j])
   else
      return value(x[j,i])
   end
end

# Si x est entier et respecte les contraintes de degre sum_{v=1}^n x[u,v)]=1 pour tout u
# Renvoie la liste des sommets formant un cycle passant par u
function find_cycle_in_integer_x(x, u)
      S = Int64[]
      #push!(S,u)
      i=u
      prev=-1
      while true
         j=1
         while  ( j==i || j==prev || edge_value(x,i,j) < 0.999 ) 
            j=j+1
         end
         push!(S,j)
         prev=i
         i=j
	 	 #println(i)
         (i!=u) || break   # si i==u alors fin de boucle
      end
      return S
end



function BandC_TSP(G)

    c = calcul_dist(G)

    LP = Model(CPLEX.Optimizer)
  
    # Setting some stat variables
    nbViolatedMengerCut_fromIntegerSep = 0
    nbViolatedMengerCut_fromFractionalSep = 0
     
    # Setting the Model LP

    #variables
    @variable(LP, x[1:G.nb_points, 1:G.nb_points], Bin)  # We will only use x_{ij} with i<j
    # objective function
    @objective(LP, Min, sum((sum(x[i, j] * c[i, j]) for j = i+1:G.nb_points) for i = 1:G.nb_points ) )
    

    # Edge constraints
    for i in 1:G.nb_points
	  @constraint(LP, (sum(x[j, i] for j in 1:i-1)) + (sum(x[i, j] for j in i+1:G.nb_points))== 2)
    end   

 
   #println(LP)
   
   # Initialization of a graph to compute min cut for the fractional separation
   
  G_sep=complete_digraph(G.nb_points)

  #################
  # our function lazySep_ViolatedMengerCut
    function lazySep_ViolatedMengerCut(cb_data)
        # cb_data is the CPLEX value of our variables for the separation algorithm
        # In the case of a LazyCst, the value is integer, but sometimes, it is more 0.99999 than 1

        # Get the x value from cb_data and round it
        xsep =zeros(Int64,G.nb_points, G.nb_points); #Int64[G.nb_points;G.nb_points]
        
        for i in 1:G.nb_points
           for j in i+1:G.nb_points
             if (callback_value(cb_data, x[i,j])>0.999)
               xsep[i,j]=1
             end
             if (callback_value(cb_data, x[i,j])<0.0001) 
               xsep[i,j]=0
             end
           end
        end
        # for i in 1:G.nb_points
        #   print(xsep[i]," ")
        # end
        # println()
        
#        violated, W = ViolatedMengerCut_IntegerSeparation(G,xsep)
        
        start=rand(1:G.nb_points)
  
        W =find_cycle_in_integer_x(xsep, start)


        if size(W,1)!=G.nb_points    # size(W) renvoie sinon (taille,)
      
           #println(W)
           
           con = @build_constraint(sum(x[i,j] for i ∈ W for j ∈ i+1:G.nb_points if j ∉ W) 
                                   + sum(x[j,i] for i ∈ W for j ∈ 1:i-1 if j ∉ W)  
                                   >= 2)
           
          #println(con)
           
           MOI.submit(LP, MOI.LazyConstraint(cb_data), con) 
           nbViolatedMengerCut_fromIntegerSep=nbViolatedMengerCut_fromIntegerSep+1
           
        end
        
    end
  #
  #################


  #################
  # our function userSep_ViolatedMengerCut
    function userSep_ViolatedMengerCut(cb_data)
        # cb_data is the CPLEX value of our variables for the separation algorithm
        # In the case of a usercut, the value is fractional or integer (and can be -0.001)

        # Get the x value from cb_data 
        xsep =zeros(Float64,G.nb_points, G.nb_points); #Int64[G.nb_points;G.nb_points]           
        for i in 1:G.nb_points
           for j in 1:i-1
               xsep[i,j]=callback_value(cb_data, x[j,i])
           end
           for j in i+1:G.nb_points
               xsep[i,j]=callback_value(cb_data, x[i,j])
           end
        end
        
       Part,valuecut=mincut(G_sep,xsep)  # Part is a vector indicating 1 and 2 for each node to be in partition 1 or 2
       
       W=Int64[]
       for i in 1:G.nb_points
          if Part[i]==1
             push!(W,i)
          end
       end
       
       if (valuecut<2.0)
      #     println(W)
           
           con = @build_constraint(sum(x[i,j] for i ∈ W for j ∈ i+1:G.nb_points if j ∉ W) 
                                   + sum(x[j,i] for i ∈ W for j ∈ 1:i-1 if j ∉ W)  
                                   >= 2)
           
      #     println(con)
           
           MOI.submit(LP, MOI.UserCut(cb_data), con) 
           nbViolatedMengerCut_fromFractionalSep=nbViolatedMengerCut_fromFractionalSep+1
         
       end
             
    end
  #
  #################

  #################
  # our function primalHeuristicTSP
    function primalHeuristicTSP(cb_data)
    
#        L=[]
#        for i in 1:nv(G)
#           push!(L, ( i , callback_value(cb_data, x[i]) ) )
#        end
#        sort!(L, by = x -> x[2])

#        sol=[0 for i in 1:nv(G)]



#        MOI.submit(LP, MOI.HeuristicSolution(cb_data), x, sol)
    
    end
  #
  #################

  #################
  # Setting callback in CPLEX
    # our lazySep_ViolatedAcyclic function sets a LazyConstraintCallback of CPLEX
    MOI.set(LP, MOI.LazyConstraintCallback(), lazySep_ViolatedMengerCut) 
    
    # our userSep_ViolatedAcyclic function sets a LazyConstraintCallback of CPLEX   
    MOI.set(LP, MOI.UserCutCallback(), userSep_ViolatedMengerCut)
    
    # our primal heuristic to "round up" a primal fractional solution
    #MOI.set(LP, MOI.HeuristicCallback(), primalHeuristicTSP)
  #
  #################


   println("Résolution du PLNE par le solveur")
   optimize!(LP)
   println("Fin de la résolution du PLNE par le solveur")
   	
   #println(solution_summary(m, verbose=true))

   status = termination_status(LP)

   # un petit affichage sympathique
   if status == JuMP.MathOptInterface.OPTIMAL
      println("Valeur optimale = ", objective_value(LP))
      println("Solution primale optimale :")
      
      #for i= 1:G.nb_points
      #   for j= i+1:G.nb_points
      #      println("x(",i,",",j,")=",value(x[i,j]))
      #   end
      #end
     
      S= find_cycle_in_integer_x(x, 1)
      push!(S,first(S))
      
      println("Temps de résolution :", solve_time(LP))
      println("Number of generated Menger Cut constraints  : ", nbViolatedMengerCut_fromIntegerSep+nbViolatedMengerCut_fromFractionalSep)
      println("   from IntegerSep : ", nbViolatedMengerCut_fromIntegerSep)
      println("   from FractionalSep :", nbViolatedMengerCut_fromFractionalSep)

      return S
    else
      println("Problème lors de la résolution")
    end
     
  
end




##########################
### Launch compact PLNE and visualization for the stable set problem
##########################

function Find_TSP_BandC(filename)

    G=Read_undirected_TSP(filename) # on load notre graphe

	@time @CPUtime S=BandC_TSP(G) # on le résout
 
	
	#println(Sol)

	println("generate pdf")
	
	filename = replace(filename, ".tsp" => "_exact_BandC")

	WritePdf_visualization_solution_ordre(G,S,filename) # on fait le pdf



end
