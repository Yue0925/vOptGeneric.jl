# vOptGeneric: part of vOptSolver for non-structured problems

[![Build Status](https://travis-ci.org/vOptSolver/vOptGeneric.jl.svg?branch=master)](https://travis-ci.org/vOptSolver/vOptGeneric.jl)
[![codecov.io](http://codecov.io/github/vOptSolver/vOptGeneric.jl/coverage.svg?branch=master)](http://codecov.io/github/vOptSolver/vOptGeneric.jl?branch=master)
[![Coverage Status](https://coveralls.io/repos/github/vOptSolver/vOptGeneric.jl/badge.svg?branch=master)](https://coveralls.io/github/vOptSolver/vOptGeneric.jl?branch=master)

**vOptSolver** is a solver of multiobjective linear optimization problems (MOMIP, MOLP, MOIP, MOCO).
This repository concerns **vOptGeneric**, the part of vOptSolver devoted to **multiobjective non-structured problems**. With vOptGeneric, the problem is expressed using JuMP algebraic language extended to multiple objectives. vOptGeneric runs on macOS, linux-ubuntu, and windows.

We suppose you are familiar with vOptSolver; if not, read first this [presentation](https://voptsolver.github.io/vOptSolver/).


## Instructions 
For a local use, a working version of:
- Julia must be ready; instructions for the installation are available [here](https://julialang.org/downloads/)
- your favorite MIP solver must be ready (GLPK is suggested); 
  instructions for the installation are available [here](https://github.com/jump-dev/JuMP.jl)
  
### Run Julia

On linux:

- open a console on your computer or in the cloud
- when the prompt is ready, type in the console `julia`

On macOS:

- locate the application `julia` and 
- click on the icon, the julia console comes to the screen

### Installation Instructions

Before your first use, 
1. run Julia and when the terminal is ready with the prompt `julia` on screen, 
2. add as follow the mandatory packages to your Julia distribution: 

```
julia> using Pkg
julia> Pkg.add("vOptGeneric")
julia> Pkg.add("JuMP")
julia> Pkg.add("GLPK")
```

That's all folk; at this point, vOptGeneric is properly installed.

### Usage Instructions

When vOptGeneric is properly installed,

1. run Julia and when the terminal is ready with the prompt `julia` on screen, 
2. invoke vOptGeneric, JuMP and the MILP solver to activate in typing in the console:
```
julia> using vOptGeneric
julia> using JuMP
julia> using GLPK
```
vOptGeneric is ready. See examples for further informations and have fun with the solver! 

## Problems available

| Problem | Description                          | Output    | Method                       | Parameter (if required)  | Name          |
|:--------|:-------------------------------------|:---------:| ---------------------------: | ------------| :--------|
| 2-IP    | bi-objective Integer Linear Program  | Y_N     | **:epsilon**                 | step = *realValue*       | ϵ-constraint  | 
| 2-IP    | bi-objective Integer Linear Program  | Y_N     | **:chalmet** or **:Chalmet** | step = *realValue*       | Chalmet       |
| 2-IP    | bi-objective Integer Linear Program  | Y_{SN}  | **:dicho** or **:dichotomy** | (none)                   | Aneja & Nair  |
| p-MIP  | multi-objective Mixed Integer Linear Program | Y_{lex} | **:lex** or **:lexico**      | (none)                   | Lexicographic |


## Examples
The folder `examples` provides (1) source code of problems ready to be solved and (2) selected datafiles into different formats.

## Limitations
No special limitation; the solving strength of vOptGeneric is currently provided by the MIP solver (GLPK, Clp/Cbc, Cbc, GUROBI, etc.) invoked.
