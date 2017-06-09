# **vOptSolver**: solver of Multiobjective linear optimization problems

### About
vOptSolver is currently supported by the ANR/DFG-14-CE35-0034-01 research project [(link)](https://voptproject.wordpress.com/). 
The version 0.2 integrates exact algorithms for computing a complete set of non-dominated points for structured and non-structured optimization problems with two [and three] objectives ([ ] = forthcoming).

### Content

- [Presentation](./README.md#presentation)
- [Features](./README.md#features)
- [Instructions](./README.md#instructions)
- [References] (./README.md#references)


### News
03-Jun-2017: The next release of vOptSolver (version 0.2) is scheduled for June 2017.

### Feedback
All bugs, feature requests, pull requests, feedback, etc., are welcome. 

### Coordinator
Prof. Dr. Xavier Gandibleux, University of Nantes - France [(contact)](http://www.univ-nantes.fr/gandibleux-x)

### How To Cite
X. Gandibleux, G. Soleilhac, A. Przybylski, S. Ruzika. 
vOptSolver: an open source software environment for multiobjective mathematical optimization.
*IFORS2017: 21st Conference of the International Federation of Operational Research Societies*. 
July 17-21, 2017 Quebec City, Canada.

## Presentation

### Aims
- Solver of multiobjective linear optimization problems for scientifics and practionners
- Easy to formulate a problem, to provide data, to solve a problem, to collect the outputs, to analyze the solutions
- Natural and intuitive use for mathematicians, informaticians, engineers

### Purposes
- Solving needs: methods and algorithms for performing numerical experiments
- Research needs: support and primitives for the development of new algorithms
- Pedagogic needs: environment for practicing of theories and algorithms

### Characteristics
- Efficient, flexible, evolutive solver
- Free, open source, multi-platform, reusing existing specifications
- Easy installation, no need of being expert in computer science

### Background
- Julia programming language [(link)](http://julialang.org/)
- JuMP algebraic language [(link)](https://jump.readthedocs.io/en/latest/)
- Usual free (GLPK) and commercial (CPLEX, GUROBI) MILP solvers

## Features

### Problems / Definition
- Non-structured problems / algebraic language: 
    -  LP: Linear Program
    -  MILP: Mixed Integer Linear Program
    -  ILP: Integer Linear Program 
- Structured problems / Application Programming Interface (API): 
    -  OSP: One machine Scheduling Problem
    -  LAP: Linear Assignment Problem 
    -  [UKP, MKP, UDFLP, SSCFLP, UMFLP, CFLP, PATHS]

### Algorithms
The solving algorithms included compute an exact complete set of non-dominated points
- Generic algorithm for structured or non-structured discrete problems: 
    - Haimes1971: epsilon-constraint method / 2ILP 
    - Aneja1979: Aneja & Nair method / 2ILP
    - Chalmet1986: Chalmet et al. method / 2ILP
    - Vincent2013: branch-and-cut / 2MILP
- Specific algorithm for structured (MOCO/MOMILP) problem: 
    - Wassenhove1980: 2OSP1980 (Julia)
    - Przybylski2008: 2LAP2008 (C)
    - [Jorge2010: 2UKP2010 (C++); Gandibleux2012: 2UDFLP2012 (C++); Delmee2017: 2UMFLP2016 (C++); Gandibleux2006: PATHS (C)]
- Algorithms and datastructures for handling non-dominated points: 
    - [AVL2016/2ILP (Julia)]

### Inputs
- Non-structured problems: 
    - direct with the provided languages (Julia, JuMP)
    - standard MOP format (ILP, MILP, LP)
    - specific problem format (MILP)
- Structured problems: 
    -  direct with the language (Julia), 
    -  specific problem format (2LAP, 2UKP, 2UFLP)

### Outputs
- Non-structured problems: 
    - standard 2MOP format (ILP, MILP, LP)
- Structured problems: 
    - specific problem format (2LAP, 2UKP, 2UFLP)

## Instructions 

### Information
- Julia is available on macOS, linux, windows for a local use or online on [JuliaBox](https://juliabox.com/) for a distant use
- vOptSolver is free, open source under [GPL] licence, tested with Julia 0.5 on macOS and linux-Ubuntu

### Installation Instructions
For a local use, a working version of:
- Julia must be ready; instructions for the installation are available [here](https://julialang.org/downloads/)
- your favorite MILP solver must be ready (mandatory for dealing with non-structured problems; GLPK is suggested); 
  instructions for the installation are available [here](http://jump.readthedocs.io/en/latest/installation.html)

Before your first local or distant use, run Julia and when the terminal is ready, add as follow the two mandatory packages to your Julia distribution: 


```
julia> Pkg.add("JuMP")
julia> Pkg.clone("vOptSolver")
```

That's all folk! 

### Usage Instructions

Install it, and then have fun with vOptSolver. 

Resources availables:

- [Examples](./examples/): 
    - problems ready to be solved
- [Documentation](./doc/): 
    - description of the Application Programming Interface 

## References

-   [Haimes1971] Y.V. Haimes, L.S. Lasdon, D.A. Wismer: 
    On a bicriterion formation of the problems of integrated system identification and system optimization.
    *IEEE Transactions on Systems, Man and Cybernetics*, Volume SMC-1, Issue 3, Pages 296-297, July 1971.
    
-   [Aneja1979] Y. P. Aneja and K. P. K. Nair: 
    Bicriteria Transportation Problem.
    *Management Science*, 25:1, 73-78 1979. 

-   [Wassenhove1980] L. N. Van Wassenhove, L. F. Gelders: 
    Solving a bicriterion scheduling problem.
    *European Journal of Operational Research*, Volume 4, Issue 1, Pages 42-48, 1980.

-   [Chalmet1986] L.G. Chalmet, L. Lemonidis, D.J. Elzinga: 
    An algorithm for the bi-criterion integer programming problem.
    *European Journal of Operational Research*, Volume 25, Issue 2, Pages 292-300, 1986.

-   [Gandibleux2006] X. Gandibleux, F. Beugnies, S. Randriamasy:  
    Martins' algorithm revisited for multi-objective shortest path problems with a MaxMin cost function. 
    *4OR: A Quarterly Journal of Operations Research*, Springer Verlag, 4 (1), pp.47-59, 2006.

-   [Przybylski2008] A. Przybylski, X. Gandibleux, M. Ehrgott: 
    Two phase algorithms for the bi-objective assignment problem.
    *European Journal of Operational Research*, Volume 185, Issue 2, Pages 509-533, 2008.

-   [Jorge2010] J. Jorge: 
    *Nouvelles propositions pour la résolution exacte du sac à dos multi-objectif unidimensionnel en variables binaires.* 
    PhD Thesis (in French), Université de Nantes - France, 2010.

-   [Gandibleux2012] X. Gandibleux, A. Przybylski , S. Bourougaa, A. Derrien, A. Grimault: 
    Computing the Efficient Frontier for the 0/1 Biobjective Uncapacitated Facility Location Problem 
    *CORS/MOPGP’2012 (10th international conference on Multiple Objective Programming and Goal Programming).* June 11-13, 2012, Niagara Falls, Canada.

-   [Vincent2013] Th. Vincent:
    *Caractérisation des solutions efficaces et algorithmes d'énumération exacts pour l'optimisation multiobjectif en variables mixtes binaires.* 
    PhD Thesis (in French), Université de Nantes - France, 2013.

-   [Delmee2017] Q. Delmée, X. Gandibleux and A. Przybylski: 
    Résolution exacte du problème de localisation de services bi-objectif sans contrainte de capacité en variables mixtes.
    *ROADEF2017 (18ème édition du congrès annuel de la Société Française de Recherche Opérationnelle et d'Aide à la Décision).* 22-24 février 2017, Metz, France.



---

Terms and acronyms used
- LP: Linear Program
- MILP: Mixed Integer Linear Program
- ILP: Integer Linear program
- MOCO: MultiObjective Combinatorial Optimization
- MOMILP: MultiObjective Mixed Integer Linear Program
- MOP: MultiObjective Program
- OSP: One machine Scheduling Problem
- LAP: Linear Assignment Problem
- UKP: Unidimensional 01 Knapsack Problem
- MKP: Multidimensional 01 Knapsack Problem
- UFLP: Uncapacitated Facility Location Problem
- UDFLP: Discrete Uncapacitated  Facility Location Problem
- SSCFLP: Single Source Capacitated Facility Location Problem
- UMFLP:  Mixed variables Uncapacitated Facility Location Problem
- CFLP: Capacitated Facility Location Problem
- PATHS: shortest paths problem
- Julia: name of the programming language
- JuMP: stands for Julia for Mathematical Optimization, a modeling language for mathematical optimization embedded in Julia
- AVL tree is a self-balancing binary search tree
- API: stands for Application Programming Interface
- GPL: stands for GNU General Public License
- GLPK: stands for GNU Linear Programming Kit, an open source solver
- CPLEX: a commercial solver
- GUROBI: a commercial solver


