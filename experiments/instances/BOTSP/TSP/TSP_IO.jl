# Code largement inspiré du travail de Florian Belhassen-Dubois
using Plots

# notre structure graphe
struct NuagePoints
	nb_points
	X
	Y
end

# on va lire un fichier tsp et construire puis renvoyer l'instance du TSP correspondant
function Read_undirected_TSP(filename)

	I = NuagePoints(0, [], [])

	open(filename) do f
			
		node_coord_section = 0
		nbp = 0
		X = Array{Float64}(undef, 0)
		Y = Array{Float64}(undef, 0)
			
		# on lit une par une nos ligne du fichier	
		for (i,line) in enumerate(eachline(f))
			
			# on sépare cette ligne en mots
			x = split(line," ") # For each line of the file, splitted using space as separator
             
			# on supprime les mots vides, en effet une espace suivi d'un autre espace renvoi le mot vide
			deleteat!(x, findall(e->e=="", x))

			
			if(node_coord_section == 0)       # If it's not a node section
					
				# si on est dans les coordonnées on va s'arrêter et remplir notre instance sinon il nous reste des labels à lire
				if(x[1] == "NODE_COORD_SECTION")
					node_coord_section = 1
				# si on est dans le label dimension, on le récupère
				elseif(x[1] == "DIMENSION")
					nbp = parse(Int, x[3])
				end
			
			# on est enfin dans nos coordonnées ! On les lit et on rempli notre instance avec
			elseif(node_coord_section == 1 && x[1] != "EOF")
				
				push!(X, parse(Float64, x[2]))
				push!(Y, parse(Float64, x[3]))
				
			else
				
				node_coord_section = 2
				
			end
		end
		
		# on construit notre nuage de points
		I = NuagePoints(nbp, X, Y)
		
	end
	return I
end

# on va visualiser nos points dans un fichier pdf dont le nom est passé en paramètre
function WritePdf_visualization_TSP(I,filename)

	filename_splitted_in_two_parts = split(filename,".") # split to remove the file extension
	filename_with_pdf_as_extension= filename_splitted_in_two_parts[1]*".pdf"
	# save to pdf
	
	# un simple plot en fonction des coordonnées 
	p = plot(I.X, I.Y, seriestype = :scatter)
	savefig(p, filename_with_pdf_as_extension)

end


# tout est dit dans le nom... distance euclidienne 
function dist(G, i, j) 
	
	return ((G.X[i] - G.X[j])^2 + (G.Y[i] - G.Y[j])^2)^(0.5)
	
end

# pour connaitre la distance de tous nos points 
function calcul_dist(G)

	c = Array{Float64}(undef, (G.nb_points, G.nb_points))
	
	for i in 1:G.nb_points
	
		for j in 1:G.nb_points
		
			c[i, j] = dist(G, i, j)
		
		end
	
	end
	
	return c
	
end

# calcule la somme des coûts de notre arête solution
function Compute_value_TSP(G, solution)
	
	X = copy(G.X)
	Y = copy(G.Y)
	
	X = resort(X, solution, G.nb_points)
	Y = resort(Y, solution, G.nb_points)
	
	push!(X, X[1])
	push!(Y, Y[1])
	
	res = ((X[1] - X[G.nb_points])^2 + (Y[1] - Y[G.nb_points])^2)^(0.5)
	for i = 1:(G.nb_points -1)	
		res = res + ((X[i] - X[(i +1)])^2 + (Y[i] - Y[(i +1)])^2)^(0.5)
	end
	
	return res

end


# permet de visualiser notre solution (un circuit / cycle) dans un fichier pdf dont le nom est spécifié en paramètres
# La solution est donnée par la liste ordonné des points à visiter commençant par 1
function WritePdf_visualization_solution_ordre(I, S, filename)

	filename_splitted_in_two_parts = split(filename,".") # split to remove the file extension
	filename_with_pdf_as_extension= filename_splitted_in_two_parts[1]*".pdf"
	# save to pdf
	
	tabX= Float16[]
	tabY= Float16[]
	
    for i in S
       push!(tabX, I.X[i])
       push!(tabY, I.Y[i])
    end
	
	# on ajoute le premier point pour le plot, c'est important sinon il manque l'arête entre 1 et n...
	push!(tabX, I.X[1])
	push!(tabY, I.Y[1])
	
	p = plot(I.X, I.Y, seriestype = :scatter,legend = false)
	plot!(p, tabX, tabY,legend = false)
	savefig(p, filename_with_pdf_as_extension)

end

