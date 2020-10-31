# include all the files
include("Interconnections.jl";
include("ModelReduction.jl");
include("Plotting.jl");
include("Realization.jl");
include("Response.jl");
include("StateSpace.jl");
include("Utilities.jl");

module RobustControl
println("This is RobustControl Toolbox.")
end # module
