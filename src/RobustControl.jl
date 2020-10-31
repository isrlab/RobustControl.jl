module RobustControl
using LinearAlgebra, PyPlot

# include all the files
include("Interconnections.jl");
include("ModelReduction.jl");
include("Plotting.jl");
include("Realization.jl");
include("Response.jl");
include("StateSpace.jl");
include("Utilities.jl");
end # module
