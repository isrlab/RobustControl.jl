module RobustControl
using LinearAlgebra, Plots, Printf, Convex, LaTeXStrings

# include all the files
include("StateSpace.jl");
include("Utilities.jl");
include("Response.jl");
include("Plotting.jl");
include("Interconnections.jl");
include("ModelReduction.jl");
include("Realization.jl");
end # module
