<img src="https://isrlab.github.io/images/julia%20robust%20control%20library.png" alt="Julia RCS Logo" width="350"/>

This is a robust control toolbox in Julia. The goal is to provide necessary functions for designing $\mathcal{H}_2/\mathcal{H}_\infty$ controllers and estimators with model uncertainty.

The synthesis is formulated as a convex optimization problem, which is solved using JuMP/Convex. The default solver is SCS. If Mosek is available, it can be used. Future version will also interface with SLICOT.

Some functions have been implemented. See testCode.jl under src/ to see basic analysis usage. 

This toolbox is not ready for use yet. Please stay tuned.
