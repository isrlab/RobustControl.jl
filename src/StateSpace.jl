using Printf

mutable struct StateSpace
    A
    B
    C
    D
    ns::Integer
    nu::Integer
    ny::Integer

    function StateSpace(A, B, C, D)::StateSpace
        this = new();

        this.ns = size(A)[1];
        this.nu = size(B)[2];
        this.ny = size(C)[1];

        # Check dimensions
        if (size(A)[1] != size(A)[2]) error("Matrix A must be square."); end
        if (size(B)[1] != size(A)[1]) error("Matrices A and B must have same number of rows."); end
        if (size(C)[2] != size(A)[2]) error("Matrices A and C must have same number of cols."); end

        this.A = A;
        this.B = B;
        this.C = C;

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

StateSpace(A, B, C) = StateSpace(A,B,C,[]);


function CheckIODimensions(sys::StateSpace,no::Integer,ni::Integer)
	if (ni<1 || ni>sys.nu) error("Invalid input specified.") end
    if (no<1 || no>sys.ny) error("Invalid output specified.") end
end

function Poles(S::StateSpace)
    lam = eigvals(S.A);
    return lam;
end

function NaturalFrequencies(S::StateSpace)
    lam = Poles(S);
    ri = real.(lam);
    ii = imag.(lam);
    om = sqrt.(ri.^2 + ii.^2);
    return(om);
end

function Damping(S::StateSpace)
    lam = Poles(S);
    ri = real.(lam);
    ii = imag.(lam);
    om = sqrt.(ri.^2 + ii.^2);
    d = abs.(ri./om);
    return(d);
end

function rifd(S::StateSpace,printFlag::Bool=true)
    lam = Poles(S);
    rr = Vector{AbstractFloat}(undef,length(lam));
    ii = Vector{AbstractFloat}(undef,length(lam));
    om = Vector{AbstractFloat}(undef,length(lam));
    d = Vector{AbstractFloat}(undef,length(lam));
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

function TransmissionZeros(S::StateSpace,no::Integer,ni::Integer)

    CheckIODimensions(S,no,ni);

    nx = size(S.A)[1];
    B = S.B[:,ni];
    C = reshape(S.C[no,:],1,length(S.C[no,:])); # Must ensure row vector
    D = S.D[no,ni];
    L = [S.A B; C D];
    II =I(nx)
    M = [II 0*B;0*C 0*D];
    e,v = eigen(L,M);
    z = filter(x-> x!=Inf,e)
    return z;
end

function DCGain(S::StateSpace)
    lam = Poles(S);
    PolesAtOrigin = filter(x->x==eps(AbstractFloat),abs.(lam));
    if isempty(PolesAtOrigin)
        return(S.D - S.C*inv(S.A)*S.B);
    else
        return(ones(S.ny,S.nu)*Inf);
    end
end

function ControllabilityMatrix(S::StateSpace)
    R = Matrix{AbstractFloat}(undef,S.ns,S.ns*S.nu);
    for i=0:(S.ns-1)
        i1 = i*S.nu+1;
        i2 = (i+1)*S.nu;
        R[:,i1:i2] = S.A^i*S.B;
    end
    return(R);
end

function ObservabilityMatrix(S::StateSpace)
    R = Matrix{AbstractFloat}(undef,S.ns*S.ny,S.ns);
    for i=0:(S.ns-1)
        i1 = i*S.ny+1;
        i2 = (i+1)*S.ny;
        R[i1:i2,:] = S.C*S.A^i;
    end
    return(R);
end

function IsControllable(S::StateSpace)
    if (rank(ControllabilityMatrix(S))==S.ns)
        return true;
    else
        return false;
    end
end

function IsObservable(S::StateSpace)
    if (rank(ObservabilityMatrix(S))==S.ns)
        return true;
    else
        return false;
    end
end
