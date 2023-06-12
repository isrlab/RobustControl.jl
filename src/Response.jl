# Code for time and frequency response of LTI systems
using LinearAlgebra

# Frequency response at user specified frequencies
function FrequencyResponse(sys::StateSpace,no::Integer,ni::Integer, Om::Array{AbstractFloat,1})::Vector{ComplexF64}
    CheckIODimensions(sys,no,ni);

    B = sys.B[:,ni];
    C = reshape(sys.C[no,:],1,length(sys.C[no,:])); # Must ensure row vector
    D = sys.D[no,ni];
    I = im*eye(size(sys.A)[1]);
    TF = Vector{ComplexF64}(undef,length(Om));
    for (i,om) in enumerate(Om)
        G = C*inv(om*I-sys.A)*B; # This is vector of one element
        TF[i] = G[1]+D; # that is why we need G[1], otherwise the addition is across different types.
    end
    return(TF);
end

# Frequency response with no user specified om. This makes an adhoc grid between [1E-3 100*wn.max], with 300 points.
function FrequencyResponse(sys::StateSpace,no::Integer,ni::Integer)::Tuple{Array{Complex{AbstractFloat},1},Array{AbstractFloat,1}}

    CheckIODimensions(sys,no,ni);

    l,ev = eigen(sys.A);
    lmax = maximum(abs.(l)); # Find highest natural frequency
    Om = logspace(1E-3,lmax*100,300); # Making adhoc decision about 100 times highest natural frequency.

    B = sys.B[:,ni];
    C = reshape(sys.C[no,:],1,length(sys.C[no,:])); # Must ensure row vector
    D = sys.D[no,ni];
    I = im*eye(size(sys.A)[1]);
    TF = Vector{ComplexF64}(undef,length(Om));
    for (i,om) in enumerate(Om)
        G = C*inv(om*I-sys.A)*B; # This is vector of one element ...
        TF[i] = G[1]+D; # ... that is why we need G[1], otherwise the addition is across different types.
    end
    return(TF,Om);
end

function StepResponse(sys::StateSpace,no::Integer,ni::Integer)::Tuple{Array{AbstractFloat,1}, Array{AbstractFloat,1}}
    l,ev = eigen(sys.A);
    lmin = minimum(abs.(l)); # Find the slowest mode
    lmax = maximum(abs.(l)); # Find the fastest mode
    tmax = 2.5*(2*pi)/lmin; # Quite adhoc -- 10 times the slowest time-constant.
    dt =  (2*pi)/lmax/50; # Quite adhoc -- Sampling is 50 times the fastest mode.
    t = Vector(0:dt:tmax);
    u = ones(length(t));
    y = ForcedResponse(sys,no,ni,t,u);
    return(t,y);
end

function StepResponse(sys::StateSpace,no::Integer,ni::Integer,t::Vector{AbstractFloat})::Vector{AbstractFloat}
    u = ones(length(t));
    y = ForcedResponse(sys,no,ni,t,u);
    return(y);
end

function ForcedResponse(sys::StateSpace,no::Integer,ni::Integer,t::Vector{AbstractFloat},u::Vector{AbstractFloat},x0::Vector{AbstractFloat}=Array{AbstractFloat}(undef,0))::Vector{AbstractFloat}
    CheckIODimensions(sys,no,ni);

    # Obtain a discrete-time system.
    dt = t[2]-t[1]; # Assumed to be uniform dt -- probably should put a check in here.
    ns = size(sys.A)[2];

    # We are only looking at SISO time response
    B = sys.B[:,ni];
    C = reshape(sys.C[no,:],1,length(sys.C[no,:])); # Must ensure row vector
    D = sys.D[no,ni];

    M = dt*[sys.A B; zeros(1,ns+1)]
    S = exp(M);
    F = S[1:ns,1:ns];
    G = S[1:ns,ns+1:ns+1];

    # Check initial condition
    if isempty(x0)
        x = zeros(ns);
    else
        x = x0;
    end

    # Time march
    y = similar(t);
    for i = 1:length(t)
        y[i] = (C*x)[1] + D*u[i];
        x = F*x + G*u[i];
    end
    return(y);
end

function ImpulseResponse(sys::StateSpace,no::Integer,ni::Integer)::Tuple{Array{AbstractFloat,1}, Array{AbstractFloat,1}}
    l,ev = eigen(sys.A);
    lmin = minimum(abs.(l)); # Find the slowest mode
    lmax = maximum(abs.(l)); # Find the fastest mode
    tmax = 2.5*(2*pi)/lmin; # Quite adhoc -- 2.5 times the slowest time-constant.
    dt =  (2*pi)/lmax/50; # Quite adhoc -- Sampling is 50 times the fastest mode.
    t = Vector(0:dt:tmax);
    y = ImpulseResponse(sys,no,ni,t);
    return(t,y);
end

function ImpulseResponse(sys::StateSpace,no::Integer,ni::Integer,t::Vector{AbstractFloat})::Vector{AbstractFloat}
    y = similar(t);
    B = sys.B[:,ni];
    C = reshape(sys.C[no,:],1,length(sys.C[no,:])); # Must ensure row vector
    for i=1:length(t)
        y[i] = (C*exp(sys.A*t[i])*B)[1];
    end
    return(y);
end

function InitialResponse(sys::StateSpace,x0::Vector{AbstractFloat})::Tuple{Vector{AbstractFloat}, Matrix{AbstractFloat}}
    l,ev = eigen(sys.A);
    lmin = minimum(abs.(l)); # Find the slowest mode
    lmax = maximum(abs.(l)); # Find the fastest mode
    tmax = 2.5*(2*pi)/lmin; # Quite adhoc -- 2.5 times the slowest time-constant.
    dt =  (2*pi)/lmax/50; # Quite adhoc -- Sampling is 50 times the fastest mode.
    t = Vector(0:dt:tmax);
    y = InitialResponse(sys,x0,t);
    return(t,y);
end

function InitialResponse(sys::StateSpace,x0::Vector{AbstractFloat},t::Vector{AbstractFloat})::Matrix{AbstractFloat}
    y = Matrix{AbstractFloat}(undef,length(t),length(x0));
    for i=1:length(t)
        y[i,:] = reshape(exp(sys.A*t[i])*x0,1,length(x0));
    end
    return(y);
end
