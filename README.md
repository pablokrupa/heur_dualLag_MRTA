# A heuristic dual Lagrangian solver for multi-agent task allocation problems

This Matlab package provides a heuristic dual Lagrangian solver for Multi-Robot Task Allocation (MRTA) problems.

The dual Lagrangian problem decomposes the central Mixed-Integer Linear Problem (MILP) of the MRTA problem into several smaller-scale MILP problems; one for each agent.
These problems can then be solved in parallel in each iteration of the dual Lagrangian method, with dual variables steering the method towards convergence of the coupling variables between agents.
However, due to the non-convex nature of MILP problems, convergence towards consensus in the coupling constraints is not guaranteed.
To solve this, we incorporate a heuristic update of the dual variables to improve convergence towards consensus.

## Installing the package

This Matlab package was developed and tested using Matlab R2024b.
The package requires a working installation of the MOSEK optimization toolbox for Matlab.
Installation instructions can be found at <https://docs.mosek.com/latest/toolbox/install-interface.html>.
This Matlab package was developed and tested using MOSEK version 10.2.6.

The package requires no installation beyond adding it to Matlab's path.

## Numerical examples included in the package

The package contains scripts with names `results_n_xx.m`, where `xx` is some integer.
These script execute examples solving random MRTA problems of increasing dimension, where `xx` indicates the number of locations of the MRTA problem.
The results of these scrips are saved in `/results`, and can be analyzed by executing `analyze_rand_batch.m` (assigning to `get_n` to the value of `xx` to be analyzed). 

