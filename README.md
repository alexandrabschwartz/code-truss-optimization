# code-truss-optimization
Some functions to optimize truss structures using vanishing constraints for stress constraints and nonslender bars. The included functions are described below.

I did test this code, but implementing the constraints correctly is a bit tricy. So expect bugs especially for 3-dimensional trusses. If you find errors or have suggestions for improvements, please let me know.

## optimizeTRUSS

This function takes a truss structure and optional some options for the solver as input arguments and then uses MPVC solvers to optimize the truss. Within the truss structure, you can include
* different load cases
* upper bounds on the bar diameters
* lower bounds on the diamter of realized bars (nonslender bars)
* bounds on the stresses on the bars in the loadCases
* bounds on the compliance of the truss in the loadCases
* bounds on the displacements of nodes inthe loadCases

## setupTRUSS_missingData

This function takes a truss structure as input and checks the truss data for completeness and consistence and tries to fill in missing data using default values. E.g. if you did not proviude bounds on the displacements, it inserts inf.

## setupTRUSS_computePotentialBars

This function  takes a truss structure as input and computes all potential bars such that at elast one end node is free and there is no other node on the bar.

# plotTRUSS_initialStructure, plotTRUSS_unloadedStructure, plotTRUSS_loadedStructure

These functions take an (optimized) truss as input and plot the intitial groundstructure, the optimized unloaded truss and the optimized truss in each load case.

# exampleTRUSS_fivebar, exampleTRUSS_tenbar

Small examples of trusses with 5/ 10 potential bars.

# exampleTRUSS_cantilever, exampleTRUSS_hook

Cantilever and hook-shaped trusses of variable dimension.


