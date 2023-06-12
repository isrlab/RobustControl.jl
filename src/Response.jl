# Code for time and frequency response of LTI systems
using LinearAlgebra

# Frequency response at user specified frequencies
function FrequencyResponse(sys::StateSpace,no::Integer,ni::Integer, Om::Vector{T}) where T <: Number
    CheckIODimensions(sys,no,ni);

    B = sys.B[:,ni];
    C = reshape(sys.C[no,:],1,length(sys.C[no,:])); # Must ensure row vector
    D = sys.D[no,ni];
    
    G = om->C*inv(om*im*I(size(sys.A)[1])-sys.A)*B .+ D;
    return vcat(G.(Om)...)
end

# Frequency response with no user specified om. This makes an adhoc grid between [1E-3 100*wn.max], with 300 points.
function FrequencyResponse(sys::StateSpace,no::Integer,ni::Integer) 

    CheckIODimensions(sys,no,ni);

    l,ev = eigen(sys.A);
    lmax = maximum(abs.(l)); # Find highest natural frequency
    Om = logspace(1E-3,lmax*100,300); # Making adhoc decision about 100 times highest natural frequency.

    return (FrequencyResponse(sys,no,ni,Om),Om)
end

function StepResponse(sys::StateSpace,no::Integer,ni::Integer,x0)
    l,ev = eigen(sys.A);
    lmin = minimum(abs.(l)); # Find the slowest mode
    lmax = maximum(abs.(l)); # Find the fastest mode
    tmax = 2.5*(2*pi)/lmin; # Quite adhoc -- 10 times the slowest time-constant.
    dt =  (2*pi)/lmax/50; # Quite adhoc -- Sampling is 50 times the fastest mode.
    t = 0:dt:tmax;
    u = ones(length(t));
    y = ForcedResponse(sys,no,ni,t,u,x0);
    return(t,y);
end

function StepResponse(sys::StateSpace,no::Integer,ni::Integer,t,x0)
    u = ones(length(t));
    y = ForcedResponse(sys,no,ni,t,u,x0);
    return(y);
end

function ForcedResponse(sys::StateSpace,no::Integer,ni::Integer,t,u,x0)
    CheckIODimensions(sys,no,ni);

    # Obtain a discrete-time system.
    dt = t[2]-t[1]; # Assumed to be uniform dt -- probably should put a check in here.
    ns = size(sys.A)[2];

    # We are only looking at SISO time response
    B = sys.B[:,ni];
    C = reshape(sys.C[no,:],1,length(sys.C[no,:])); # Must ensure row vector
    D = sys.D[no,ni];
    x = x0;

    # Time march
    y = similar(t);
    for i in eachindex(t)
        y[i] = (C*x)[1] + D*u[i];
        if i < length(t)
            dt = t[i+1] - t[i]
            M = dt*[sys.A B; zeros(1,ns+1)]
            S = exp(M);
            F = S[1:ns,1:ns];
            G = S[1:ns,ns+1:ns+1];
            x = F*x + G*u[i];
        end
    end
    return(y);
end

function ImpulseResponse(sys::StateSpace,no::Integer,ni::Integer)
    l,ev = eigen(sys.A);
    lmin = minimum(abs.(l)); # Find the slowest mode
    lmax = maximum(abs.(l)); # Find the fastest mode
    tmax = 2.5*(2*pi)/lmin; # Quite adhoc -- 2.5 times the slowest time-constant.
    dt =  (2*pi)/lmax/50; # Quite adhoc -- Sampling is 50 times the fastest mode.
    t = 0:dt:tmax;
    y = ImpulseResponse(sys,no,ni,t);
    return t,y
end

function ImpulseResponse(sys::StateSpace,no::Integer,ni::Integer,t)
    B = sys.B[:,ni];
    C = reshape(sys.C[no,:],1,length(sys.C[no,:])); # Must ensure row vector
    y = t->C*exp(sys.A*t)*B;
    return vcat(y.(t)...)
end

function InitialResponse(sys::StateSpace,x0)
    l,ev = eigen(sys.A);
    lmin = minimum(abs.(l)); # Find the slowest mode
    lmax = maximum(abs.(l)); # Find the fastest mode
    tmax = 2.5*(2*pi)/lmin; # Quite adhoc -- 2.5 times the slowest time-constant.
    dt =  (2*pi)/lmax/50; # Quite adhoc -- Sampling is 50 times the fastest mode.
    t = Vector(0:dt:tmax);
    y = InitialResponse(sys,x0,t);
    return(t,y);
end

function InitialResponse(sys::StateSpace,x0,t)
    y = t-> exp(sys.A*t)*x0;
    return vcat(y.(t)'...)
end
