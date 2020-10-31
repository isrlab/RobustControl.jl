# Utilities for Robust Control Toolbox
using LinearAlgebra

function eye(n)::Matrix{Float64}
	return Matrix{Float64}(I, n, n)
end

function linspace(x1::Real,x2::Real,N::Integer)::Vector{Float64}
    if N>1
        dx = (x2-x1)/(N-1);
        return(x1:dx:x2);
    else
        return([]);
    end
end

function logspace(x1::Real,x2::Real,N::Integer)::Vector{Float64}
    if min(x1,x2) <= 0
        error("logspace is undefined for negative values");
        return([]);
    else
        x = linspace(log(x1),log(x2),N);
        return(exp.(x));
    end
end
