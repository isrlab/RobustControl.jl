module RobustControl
using LinearAlgebra, PyPlot, Printf

# include all the files
include("StateSpace.jl");
include("Utilities.jl");
include("Response.jl");
include("Plotting.jl");
include("Interconnections.jl");
include("ModelReduction.jl");
include("Realization.jl");
end # module
