using Plots
using Shuffle

include("TSP_IO.jl")


# fonction qui construit une solution grâce à la méthode du plus proche voisin 
function plus_proche_voisin(G) 

	X = copy(G.X)
	Y = copy(G.Y)
	
	# permer de se souvenir du nom du point dont on parle 
	I = Vector(1:(G.nb_points -1))
	I = I .+ 1

	# arbitrairement on part du point 1 
	cx = pop!(X)
	cy = pop!(Y)
	
	res = Int64[]
	push!(res, 1)
	
	# pour chaque point courant
	for i in 1:G.nb_points -1
		
		ppv = 1
		pcd = (X[1] - cx)^2 + (Y[1] - cy)^2
		
		# on cherche son plus proche voisin
		for j in 2:G.nb_points -i
		
			dist = (X[j] - cx)^2 + (Y[j] - cy)^2
			
			if(pcd > dist)
				
				pcd = dist
				ppv = j
			
			end
		
		end
		
		# puis on ajoute ce fameux plus proche voisin
		cx = X[ppv]
		cy = Y[ppv]
		
		push!(res, I[ppv])
		deleteat!(X, ppv)
		deleteat!(Y, ppv)
		deleteat!(I, ppv)
	
	end
	
	val=Compute_value_TSP(G, res)
	println("Solution Heuristique S=",res)
	println("Valeur: ",val)
	println()
	
	return res
	
end



# permet de faire une execution et de générer un pdf
function Find_TSP_ppv(filename)
	
	I = Read_undirected_TSP(filename)
	
	s = plus_proche_voisin(I)
	
	filename = replace(filename, ".tsp" => "_ppv")
	
	WritePdf_visualization_solution_ordre(I, s, filename)
	
	
end

