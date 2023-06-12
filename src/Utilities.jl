# Utilities for Robust Control Toolbox
using LinearAlgebra


function linspace(x1::T,x2::T,N::Integer)::Vector{T} where T<: Real
   return LinRange(x1,x2,N);
end

function logspace(x1::T,x2::T,N::Integer) where T<: Real
    if min(x1,x2) <= 0
        error("logspace is undefined for negative values");
        return([]);
    else
        x = linspace(log(x1),log(x2),N);
        return(exp.(x));
    end
end
