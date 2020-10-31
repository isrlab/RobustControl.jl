mutable struct StateSpace
    A::Matrix{Float64}
    B::Matrix{Float64}
    C::Matrix{Float64}
    D::Matrix{Float64}
    ns::Integer
    nu::Integer
    ny::Integer

    function StateSpace(A::Array, B::Array, C::Array, D::Array=Array{Float64,2}(undef,0,2))::StateSpace
        this = new();

        this.ns = size(A)[1];
        this.nu = size(B)[2];
        this.ny = size(C)[1];

        # Check dimensions
        if (size(A)[1] != size(A)[2]) error("Matrix A must be square."); end
        if (size(B)[1] != size(A)[1]) error("Matrices A and B must have same number of rows."); end
        if (size(C)[2] != size(A)[2]) error("Matrices A and C must have same number of cols."); end

        this.A = float(A);
        this.B = float(B);
        this.C = float(C);

        if isempty(D)
            this.D = zeros(this.ny,this.nu);
        else
            if (size(D)[1] != size(C)[1]) error("Matrices C and D must have same number of rows."); end
            if (size(D)[2] != size(B)[2]) error("Matrices B and D must have same number of cols."); end
            this.D = float(D);
        end
        return(this);
    end
end

function CheckIODimensions(sys::StateSpace,no::Integer,ni::Integer)
	if (ni<1 || ni>sys.nu) error("Invalid input specified.") end
    if (no<1 || no>sys.ny) error("Invalid output specified.") end
end

function Poles(S::StateSpace)
    lam = eigvals(S.A);
    return lam;
end

function NaturalFrequencies(S::StateSpace)::Vector{Float64}
    lam = Poles(S);
    ri = real.(lam);
    ii = imag.(lam);
    om = sqrt.(ri.^2 + ii.^2);
    return(om);
end

function Damping(S::StateSpace)::Vector{Float64}
    lam = Poles(S);
    ri = real.(lam);
    ii = imag.(lam);
    om = sqrt.(ri.^2 + ii.^2);
    d = abs.(ri./om);
    return(d);
end

function rifd(S::StateSpace,printFlag::Bool=false)
    lam = Poles(S);
    rr = Vector{Float64}(undef,length(lam));
    ii = Vector{Float64}(undef,length(lam));
    om = Vector{Float64}(undef,length(lam));
    d = Vector{Float64}(undef,length(lam));
    if printFlag
        @printf("   Real \t   Imag \t  Omega \t  Damping\n")
    end
    for (i,l) in enumerate(lam)
        rr[i] = real(l)
        ii[i] = imag(l)
        om[i] = sqrt(rr[i]^2+ii[i]^2)
        d[i]= abs(rr[i]/om[i])
        if printFlag
            @printf("%+.6f \t %+.6f \t %.6f \t %.6f\n",rr[i],ii[i],om[i],d[i])
        end
    end
    return(rr,ii,om,d)
end

function TransmissionZeros(S::StateSpace,no::Integer,ni::Integer)::Vector{Float64}

    CheckIODimensions(S,no,ni);

    B = S.B[:,ni];
    C = reshape(S.C[no,:],1,length(S.C[no,:])); # Must ensure row vector
    D = S.D[no,ni];
    L = [S.A B; C D];
    II = Matrix{Float64}(I, size(S.A)[1], size(S.A)[2]);
    M = [II 0*B;0*C 0*D];
    e,v = eigen(L,M);
    z = filter(x-> x!=Inf,e)
    return z;
end

function DCGain(S::StateSpace)::Matrix{Float64}
    lam = Poles(S);
    PolesAtOrigin = filter(x->x==eps(Float64),abs.(lam));
    if isempty(PolesAtOrigin)
        return(S.D - S.C*inv(S.A)*S.B);
    else
        return(ones(S.ny,S.nu)*Inf);
    end
end

function ControllabilityMatrix(S::StateSpace)::Matrix{Float64}
    R = Matrix{Float64}(undef,S.ns,S.ns*S.nu);
    for i=0:(S.ns-1)
        i1 = i*S.nu+1;
        i2 = (i+1)*S.nu;
        R[:,i1:i2] = S.A^i*S.B;
    end
    return(R);
end

function ObservabilityMatrix(S::StateSpace)::Matrix{Float64}
    R = Matrix{Float64}(undef,S.ns*S.ny,S.ns);
    for i=0:(S.ns-1)
        i1 = i*S.ny+1;
        i2 = (i+1)*S.ny;
        R[i1:i2,:] = S.C*S.A^i;
    end
    return(R);
end

function IsControllable(S::StateSpace)::Bool
    if (rank(ControllabilityMatrix(S))==S.ns)
        return true;
    else
        return false;
    end
end

function IsObservable(S::StateSpace)::Bool
    if (rank(ObservabilityMatrix(S))==S.ns)
        return true;
    else
        return false;
    end
end
