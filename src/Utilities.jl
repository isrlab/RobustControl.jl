# Utilities for Robust Control Toolbox
using LinearAlgebra

function eye(n)::Matrix{AbstractFloat}
	return Matrix{AbstractFloat}(I, n, n)
end

function linspace(x1::AbstractFloat,x2::AbstractFloat,N::Integer)::Vector{AbstractFloat}
   return LinRange(x1,x2,N);
end

function logspace(x1::AbstractFloat,x2::AbstractFloat,N::Integer)::Vector{AbstractFloat}
    if min(x1,x2) <= 0
        error("logspace is undefined for negative values");
        return([]);
    else
        x = linspace(log(x1),log(x2),N);
        return(exp.(x));
    end
end
