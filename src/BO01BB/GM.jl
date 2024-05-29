# include("../../../GravityMachine/src/GMmain.jl")

# ==============================================================================
# The gravity machine (Man of Steel) -> to terraform the world


const verbose = false
const graphic = false

verbose ? println("""\nAlgorithme "Gravity machine" --------------------------------\n""") : nothing

verbose ? println("-) Active les packages requis\n") : nothing

using JuMP, CPLEX, PyPlot, Printf, Random

verbose ? println("  Fait \n") : nothing

generateurVisualise = -1

# ==============================================================================
include("../../../GravityMachine/src/GMdatastructures.jl") # types, datastructures and global variables specially defined for GM
include("../../../GravityMachine/src/GMparsers.jl")        # parsers of instances and non-dominated points
include("../../../GravityMachine/src/GMgenerators.jl")     # compute the generators giving the L bound set
include("../../../GravityMachine/src/GMjumpModels.jl")     # JuMP models for computing relaxed optima of the SPA 
include("../../../GravityMachine/src/GMrounding.jl")       # Startegies for rounding a LP-solution to a 01-solution
include("../../../GravityMachine/src/GMprojection.jl")     # JuMP models for computing the projection on the polytope of the SPA
include("../../../GravityMachine/src/GMmopPrimitives.jl")  # usuals algorithms in multiobjective optimization
include("../../../GravityMachine/src/GMperturbation.jl")   # routines dealing with the perturbation of a solution when a cycle is detected
include("../../../GravityMachine/src/GMquality.jl")        # quality indicator of the bound set U generated


# ==============================================================================
# Ajout d'une solution relachee initiale a un generateur

function ajouterX0!(vg::Vector{tGenerateur}, k::Int64, s::tSolution)

    vg[k].sRel = deepcopy(s) # met en place le generateur \bar{x}^k
    vg[k].sPrj = deepcopy(s) # le generateur est la premiere projection \bar{x}^{k,0}
    return nothing
end


# ==============================================================================
# Ajout d'une solution entiere (arrondie ou perturbee) a un generateur

function ajouterXtilde!(vg::Vector{tGenerateur}, k::Int64, x::Vector{Float64}, y::Vector{Float64})

    vg[k].sInt.x = copy(x)
    vg[k].sInt.y = copy(y)
    return nothing
end


# ==============================================================================
# Ajout d'une solution fractionnaire (projetee) a un generateur

function ajouterXbar!(vg::Vector{tGenerateur}, k::Int64, x::Vector{Float64}, y::Vector{Float64})

    vg[k].sPrj.x = copy(x)
    vg[k].sPrj.y = copy(y)
    return nothing
end


# ==============================================================================
# Elabore 2 ensembles d'indices selon que xTilde[i] vaut 0 ou 1

function split01(xTilde::Array{Float64,1})

   indices0 = (Int64)[]
   indices1 = (Int64)[]

   for i=1:length(xTilde)
       if isapprox(xTilde[i], 0.0, atol=1e-3)
           push!(indices0,i)
       else
           push!(indices1,i)
       end
        # # todo : sensibility of assigning 1 
        # if isapprox(xTilde[i], 1.0, atol=1e-3)
        #     push!(indices1,i)
        # else
        #     push!(indices0,i)
        # end
    
    end

   return indices0, indices1
end


# ==============================================================================
# test si une solution est admissible en verifiant si sa relaxation lineaire
# conduit a une solution entiere

function estAdmissible(x::Vector{Float64})

    admissible = true
    i=1
    while admissible && i<=length(x)
        if round(x[i], digits=3)!=0.0 && round(x[i], digits=3)!=1.0
            admissible = false
        end
        i+=1
    end
    return admissible
end


# ==============================================================================
# calcule la performance z d'une solution x sur les 2 objectifs

function evaluerSolution(x::Vector{Float64}, c::Matrix{Float64})

    z1 = c[1, 1]; z2 = c[2, 1]
    for i in 1:length(x)
        z1 += x[i] * c[1, i+1]
        z2 += x[i] * c[2, i+1]
    end
    return z1,z2
end


# ==============================================================================
# nettoyage des valeurs des variables d'une solution x relachee sur [0,1]

function nettoyageSolution!(x::Vector{Float64})

    nbvar=length(x)
    for i in 1:nbvar
        if     round(x[i], digits=3) == 0.0
                   x[i] = 0.0
        elseif round(x[i], digits=3) == 1.0
                   x[i] = 1.0
        else
                   x[i] = round(x[i], digits=3)
        end
    end
end


# ==============================================================================
# predicat : verifie si une solution entiere est realisable
function isFeasible(vg::Vector{tGenerateur}, k::Int64)
    #verbose && vg[k].sFea == true ? println("   feasible") : nothing
    return (vg[k].sFea == true)
end


# ==============================================================================
# predicat : verifie si le nombre d'essai maximum a ete tente
function isFinished(trial::Int64, maxTrial::Int64)
#    verbose && trial > maxTrial ? println("   maxTrial") : nothing
    return (trial > maxTrial)
end


# ==============================================================================
# predicat : verifie si le budget de calcul maximum a ete consomme
function isTimeout(temps, maxTime)
#    verbose && time()- temps > maxTime ? println("   maxTime") : nothing
    return (time()- temps > maxTime)
end


# ==============================================================================
# elabore pC le pointeur du cone ouvert vers L

function elaborePointConeOuvertversL(vg::Vector{tGenerateur}, k::Int64, pB::tPoint, pA::tPoint)

    # recupere les coordonnees du point projete
    pC=tPoint(vg[k].sPrj.y[1], vg[k].sPrj.y[2])

    # etablit le point nadir pN au depart des points pA et pB adjacents au generateur k
    pN = tPoint( pA.x , pB.y )

#    print("Coordonnees du cone 2 : ")
#    @show pC, pN

    # retient pN si pC domine pN (afin d'ouvrir le cone)
    if pC.x < pN.x  &&  pC.y < pN.y
        # remplace pC par pN
        pC=tPoint( pA.x , pB.y )
    end

    return pC
end


# ==============================================================================
#= Retourne un booléen indiquant si un point se trouve dans un secteur défini dans
  le sens de rotation trigonométrique (repère X de gauche à droite, Y du haut vers
  le bas).
  https://www.stashofcode.fr/presence-dun-point-dans-un-secteur-angulaire/#more-328
  M    Point dont la position est à tester (point resultant a tester)
  O    Point sommet du secteur (point generateur)
  A    Point de départ du secteur (point adjacent inferieur)
  B    Point d'arrivée du secteur (point adjacent superieur)
  sortie : Booléen indiquant si le point est dans le secteur ou non.

  Exemple :

  B=point(2.0,1.0)
  O=point(2.5,2.5)
  A=point(5.0,5.0)

  M=point(5.0,4.0)
  inSector(M, O, A, B)
=#

function inSector(M, O, A, B)

    cpAB = (A.y - O.y) * (B.x - O.x) - (A.x - O.x) * (B.y - O.y)
    cpAM = (A.y - O.y) * (M.x - O.x) - (A.x - O.x) * (M.y - O.y)
    cpBM = (B.y - O.y) * (M.x - O.x) - (B.x - O.x) * (M.y - O.y)

    if (cpAB > 0)
        if ((cpAM > 0) && (cpBM < 0))
            return true
        else
            return false
        end
    else
        if (!((cpAM < 0) && (cpBM > 0)))
            return true
        else
            return false
        end
    end
end

function inCone(pOrg, pDeb, pFin, pCur)
    # pOrg : point origine du cone (la ou il est pointe)
    # pDeb : point depart du cone (point du rayon [pOrg,pDeb])
    # pFin : point final du cone (point du rayon [pOrg,pFin])
    # pCur : point courant a tester
    # retourne VRAI si pCur est dans le cone pDeb-pFin-pOrg, FAUX sinon

    cp_pDeb_pFin = (pDeb.x - pOrg.x) * (pFin.y - pOrg.y) - (pDeb.y - pOrg.y) * (pFin.x - pOrg.x)
    cp_pDeb_pCur = (pDeb.x - pOrg.x) * (pCur.y - pOrg.y) - (pDeb.y - pOrg.y) * (pCur.x - pOrg.x)
    cp_pFin_pCur = (pFin.x - pOrg.x) * (pCur.y - pOrg.y) - (pFin.y - pOrg.y) * (pCur.x - pOrg.x)

    if (cp_pDeb_pFin > 0)
        if ((cp_pDeb_pCur >= 0) && (cp_pFin_pCur <= 0))
            return true
        else
            return false
        end
    else
        if (!((cp_pDeb_pCur < 0) && (cp_pFin_pCur > 0)))
            return true
        else
            return false
        end
    end
end

function inCone1VersZ(pOrg, pDeb, pFin, pCur)
    return inCone(pOrg, pDeb, pFin, pCur)
end

function inCone2Vers0(pOrg, pDeb, pFin, pCur)
    return !inCone(pOrg, pDeb, pFin, pCur)
end


# ==============================================================================
# Selectionne les points pour le cone pointe sur le generateur k (pCour) et ouvert vers Y
function selectionPoints(vg::Vector{tGenerateur}, k::Int64)
    nbgen = size(vg,1)
    if k==1
        # premier generateur (point predecesseur fictif)
        pPrec = tPoint(vg[k].sRel.y[1], vg[k].sRel.y[2]+1.0)
        pCour = tPoint(vg[k].sRel.y[1], vg[k].sRel.y[2])
        pSuiv = tPoint(vg[k+1].sRel.y[1], vg[k+1].sRel.y[2])
    elseif k==nbgen
        # dernier generateur (point suivant fictif)
        pPrec = tPoint(vg[k-1].sRel.y[1], vg[k-1].sRel.y[2])
        pCour = tPoint(vg[k].sRel.y[1], vg[k].sRel.y[2])
        pSuiv = tPoint(vg[k].sRel.y[1]+1.0, vg[k].sRel.y[2])
    else
        # generateur non extreme
        pPrec = tPoint(vg[k-1].sRel.y[1], vg[k-1].sRel.y[2])
        pCour = tPoint(vg[k].sRel.y[1], vg[k].sRel.y[2])
        pSuiv = tPoint(vg[k+1].sRel.y[1], vg[k+1].sRel.y[2])
    end
#    print("Coordonnees du cone 1 : ")
#    @show pPrec, pCour, pSuiv
    return pPrec, pCour, pSuiv
end


# ==============================================================================
# Calcule la direction d'interet du nadir vers le milieu de segment reliant deux points generateurs
function calculerDirections(L::Vector{tSolution}, vg::Vector{tGenerateur})
   # function calculerDirections(L, vg::Vector{tGenerateur})

    nbgen = size(vg,1)
    for k in 2:nbgen

        n1 = L[end].y[1]
        n2 = L[1].y[2]

        x1,y1 = vg[k-1].sRel.y[1], vg[k-1].sRel.y[2]
        x2,y2 = vg[k].sRel.y[1], vg[k].sRel.y[2]
        xm=(x1+x2)/2.0
        ym=(y1+y2)/2.0
        Δx = abs(n1-xm)
        Δy = abs(n2-ym)
        λ1 =  1 - Δx / (Δx+Δy)
        λ2 =  1 - Δy / (Δx+Δy)
        verbose ? @printf("  x1= %7.2f   y1= %7.2f \n",x1,y1) : nothing
        verbose ? @printf("  x2= %7.2f   y2= %7.2f \n",x2,y2) : nothing
        verbose ? @printf("  Δx= %7.2f    Δy= %7.2f \n",Δx,Δy) : nothing
        verbose ? @printf("  λ1= %6.5f    λ2= %6.5f \n",λ1,λ2) : nothing
        # plot(n1, n2, xm, ym, linestyle="-", color="blue", marker="+")
        # annotate("",
        #          xy=[xm;ym],# Arrow tip
        #          xytext=[n1;n2], # Text offset from tip
        #          arrowprops=Dict("arrowstyle"=>"->"))
        # println("")
    end

end


# ==============================================================================
# Calcule la direction d'interet du nadir vers un point generateur
function calculerDirections2(L::Vector{tSolution}, vg::Vector{tGenerateur})
    #function calculerDirections2(L, vg::Vector{tGenerateur})

    nbgen = size(vg,1)
    λ1=Vector{Float64}(undef, nbgen)
    λ2=Vector{Float64}(undef, nbgen)
    for k in 1:nbgen

        n1 = L[end].y[1]
        n2 = L[1].y[2]

        xm=vg[k].sRel.y[1]
        ym=vg[k].sRel.y[2]
        Δx = abs(n1-xm)
        Δy = abs(n2-ym)
        λ1[k] =  1 - Δx / (Δx+Δy)
        λ2[k] =  1 - Δy / (Δx+Δy)
        verbose ? @printf("  k= %3d   ",k) : nothing
        verbose ? @printf("  xm= %7.2f   ym= %7.2f ",xm,ym) : nothing
        verbose ? @printf("  Δx= %8.2f    Δy= %8.2f ",Δx,Δy) : nothing
        verbose ? @printf("  λ1= %6.5f    λ2= %6.5f \n",λ1[k],λ2[k]) : nothing
        # if generateurVisualise == -1 
        #     # affichage pour tous les generateurs
        #     plot(n1, n2, xm, ym, linestyle="-", color="blue", marker="+")
        #     annotate("",
        #              xy=[xm;ym],# Arrow tip
        #              xytext=[n1;n2], # Text offset from tip
        #              arrowprops=Dict("arrowstyle"=>"->"))        
        # elseif generateurVisualise == k
        #     # affichage seulement pour le generateur k
        #     plot(n1, n2, xm, ym, linestyle="-", color="blue", marker="+")
        #     annotate("",
        #              xy=[xm;ym],# Arrow tip
        #              xytext=[n1;n2], # Text offset from tip
        #              arrowprops=Dict("arrowstyle"=>"->"))        
        # end 
        #println("")
    end
    return λ1, λ2
end
 

function GM( model::JuMP.Model,  x::Array{JuMP.VariableRef}, 
    c::Matrix{Float64},
    tailleSampling::Int64,
    maxTrial::Int64,
    maxTime::Int64
  )

@assert tailleSampling>=3 "Erreur : Au moins 3 sont requis"

verbose ? @printf("0) instance et parametres \n\n") : nothing
verbose ? println("  tailleSampling = $tailleSampling | maxTrial = $maxTrial | maxTime = $maxTime\n\n") : nothing

nbvar = length(x)
nbobj = 2

# structure pour les points qui apparaitront dans l'affichage graphique
d = tListDisplay([],[], [],[], [],[], [],[], [],[], [],[], [],[])

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
verbose ? @printf("1) calcule les etendues de valeurs sur les 2 objectifs\n\n") : nothing

#todo (model): calcule la valeur optimale relachee de f1 seule et le point (z1,z2) correspondant
f1RL, xf1RL = computeLinearRelax(model, x, c, -1, 1) # opt fct 1
minf1RL, maxf2RL = evaluerSolution(xf1RL, c)

# calcule la valeur optimale relachee de f2 seule et le point (z1,z2) correspondant
f2RL, xf2RL = computeLinearRelax(model, x, c, -1, 2) # opt fct 2
maxf1RL, minf2RL = evaluerSolution(xf2RL, c)

verbose ? @printf("  f1_min=%8.2f ↔ f1_max=%8.2f (Δ=%.2f) \n",minf1RL, maxf1RL, maxf1RL-minf1RL) : nothing
verbose ? @printf("  f2_min=%8.2f ↔ f2_max=%8.2f (Δ=%.2f) \n\n",minf2RL, maxf2RL, maxf2RL-minf2RL) : nothing


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
verbose ? @printf("2) calcule les generateurs par e-contrainte alternant minimiser z1 et z2\n\n") : nothing

nbgen, L = calculGenerateurs(model, x, c, tailleSampling, minf1RL, maxf2RL, maxf1RL, minf2RL, d)

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# allocation de memoire pour la structure de donnees -----------------------

vg = allocateDatastructure(nbgen, nbvar, nbobj)

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
verbose ? @printf("3) place L dans structure et verifie l'admissibilite de chaque generateur\n\n") : nothing

for k=1:nbgen

verbose ? @printf("  %2d  : [ %8.2f , %8.2f ] ", k, L[k].y[1], L[k].y[2]) : nothing

# copie de l'ensemble bornant inferieur dans la stru de donnees iterative ---
ajouterX0!(vg, k, L[k])

# test d'admissibilite et marquage de la solution le cas echeant -------
if estAdmissible(vg[k].sRel.x)
   ajouterXtilde!(vg, k, vg[k].sRel.x, L[k].y)
   vg[k].sFea   = true
   verbose ? @printf("→ Admissible \n") : nothing
   # archive le point obtenu pour les besoins d'affichage    
   if generateurVisualise == -1 
       # archivage pour tous les generateurs
       push!(d.XFeas,vg[k].sInt.y[1])
       push!(d.YFeas,vg[k].sInt.y[2])
   elseif generateurVisualise == k
       # archivage seulement pour le generateur k
       push!(d.XFeas,vg[k].sInt.y[1])
       push!(d.YFeas,vg[k].sInt.y[2])
   end 
else
   vg[k].sFea   = false
   verbose ? @printf("→ x          \n") : nothing
end

end
verbose ? println("") : nothing

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# calcule les directions (λ1,λ2) pour chaque generateur a utiliser lors des projections
λ1,λ2 = calculerDirections2(L,vg)

# ==========================================================================

verbose ? @printf("4) terraformation generateur par generateur \n\n") : nothing

for k in [i for i in 1:nbgen if !isFeasible(vg,i)]
temps = time()
trial = 0
H =(Vector{Float64})[]

#perturbSolution30!(vg,k,c1,c2,d)

# rounding solution : met a jour sInt dans vg --------------------------
#roundingSolution!(vg,k,c1,c2,d)  # un cone
#roundingSolutionnew24!(vg,k,c1,c2,d) # deux cones
roundingSolutionNew23!(vg,k,c,d) # un cone et LS sur generateur

push!(H,[vg[k].sInt.y[1],vg[k].sInt.y[2]])
verbose ? println("   t=",trial,"  |  Tps=", round(time()- temps, digits=4)) : nothing

while !(t1=isFeasible(vg,k)) && !(t2=isFinished(trial, maxTrial)) && !(t3=isTimeout(temps, maxTime))

   trial+=1

   # projecting solution : met a jour sPrj, sInt, sFea dans vg --------
   projectingSolution!(model, x, vg,k, c, λ1,λ2,d)
   verbose ? println("   t=",trial,"  |  Tps=", round(time()- temps, digits=4)) : nothing

   if !isFeasible(vg,k)

       # rounding solution : met a jour sInt dans vg --------------------------
       #roundingSolution!(vg,k,c1,c2,d)
       #roundingSolutionnew24!(vg,k,c1,c2,d)
       roundingSolutionNew23!(vg,k,c,d)
       verbose ? println("   t=",trial,"  |  Tps=", round(time()- temps, digits=4)) : nothing

       # test detection cycle sur solutions entieres ------------------
       cycle = [vg[k].sInt.y[1],vg[k].sInt.y[2]] in H
       if (cycle == true)
           verbose ? println("CYCLE!!!!!!!!!!!!!!!") : nothing
           # perturb solution
           perturbSolution30!(vg,k,c,d)
       end
       push!(H,[vg[k].sInt.y[1],vg[k].sInt.y[2]])

   end
end
if t1
   verbose ? println("   feasible \n") : nothing
elseif t2
   verbose ? println("   maxTrial \n") : nothing
elseif t3
   verbose ? println("   maxTime \n") : nothing
end


end

verbose ? println("") : nothing

# ==========================================================================

verbose ? @printf("5) Extraction des resultats\n\n") : nothing


for k=1:nbgen
verbose ? @printf("  %2d  : [ %8.2f , %8.2f ] ", k, vg[k].sInt.y[1],vg[k].sInt.y[2]) : nothing
# test d'admissibilite et marquage de la solution le cas echeant -------
if vg[k].sFea
   verbose ? @printf("→ Admissible \n") : nothing
else
   verbose ? @printf("→ x          \n") : nothing
end
end


    # U = []
    # for k=1:nbgen
    #     verbose ? @printf("  %2d  : [ %8.2f , %8.2f ] ", k, vg[k].sInt.y[1],vg[k].sInt.y[2]) : nothing
    #     # test d'admissibilite et marquage de la solution le cas echeant -------
    #     if vg[k].sFea
    #         push!(U, vg[k].sInt.y)
    #     end
    # end

    return vg, nbgen
end
