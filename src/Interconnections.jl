# Interconnections of linear systems

# include("StateSpace.jl");
# include("Utilities.jl");

function Negate(S::StateSpace)::StateSpace
    return StateSpace(S.A,S.B,-S.C,-S.D);
end

function Subsystem(S::StateSpace, inputs::Vector{Int64}, outputs::Vector{Int64})::StateSpace
    # Check IO dimensions.
    if (!isempty(filter(x->(x>S.nu),inputs))) error("Input index exceeds input dimension."); end
    if (!isempty(filter(x->(x>S.ny),outputs))) error("Output index exceeds input dimension."); end

    B = S.B[:,inputs];
    C = S.C[outputs,:];
    D = S.D[outputs,inputs];
    
    return StateSpace(S.A,B,C,D);
end

function Parallel(S1::StateSpace,S2::StateSpace)::StateSpace
    # Check io dimensions.
    if (S1.nu!=S2.nu || S1.ny!=S2.ny)
        error("Input and output dimensions must match.")
    end

    A = [S1.A zeros(S1.ns,S2.ns);zeros(S2.ns,S1.ns) S2.A];
    B = [S1.B;S2.B];
    C = [S1.C S2.C];
    D = (S1.D+S2.D);
    return StateSpace(A,B,C,D);
end

function Series(S1::StateSpace,S2::StateSpace)::StateSpace
    # Check io dimensions.
    if (S1.ny!=S2.nu)
        error("Output of S1 must match input of S2.");
    end

    A = [S1.A zeros(S1.ns,S2.ns);S2.B*S1.C S2.A];
    B = [S1.B;S2.B*S1.D];
    C = [S2.D*S1.C S2.C];
    D = S2.D*S1.D;

    return StateSpace(A,B,C,D);
end

function Feedback(S1::StateSpace,S2::StateSpace)::StateSpace
    # Check io dimensions -- or provide io options for feedback.
    if (S1.ny != S2.nu || S1.nu != S2.ny)
        error("Incompatible I/O dimensions.")
    end

    A1 = S1.A; B1 = S1.B; C1 = S1.C; D1 = S1.D;
    A2 = S2.A; B2 = S2.B; C2 = S2.C; D2 = S2.D;

    I = eye(S2.ny);
    X = I+D2*D1;

    if (det(X) == 0)
        error("Feedback interconnection is illposed.");
    end

    invX = inv(X);
    Cb = [D2*C1 C2];

    AA1 = [A1 zeros(S1.ns,S2.ns)] - B1*inv(X)Cb;
    AA2 = [B2*C1 A2];
    AA = [AA1;AA2];

    BB1 = B1*invX;
    BB2 = B2*D1;
    BB = [BB1;BB2];

    CC = [C1 zeros(S1.ny,S2.ns)] - D1*invX*Cb;
    DD = D1*inv(X);

    return StateSpace(AA,BB,CC,DD);
end

# Need to implement partial interconnections.
# Eventually, these will replace the above code.

function Parallel(S1::StateSpace,S2::StateSpace,
    in1::Vector{Int64}, in2::Vector{Int64},
    out1::Vector{Int64}, out2::Vector{Int64})::StateSpace
    #                               +------+
    #                 v1 ---------->|      |----------> z1
    #                               |  S1  |
    #                        u1 +-->|      |---+ y1
    #                           |   +------+   |
    #                  u ------>+              O------> y
    #                           |   +------+   |
    #                        u2 +-->|      |---+ y2
    #                               |  S2  |
    #                 v2 ---------->|      |----------> z2
    #                               +------+
    # New system maps (v1;u;v2) -> (z1;y;z2)
    #   -- in1 defines channel u1
    #   -- in2 defines channel u2
    #   -- out1 defines channel y1
    #   -- out2 defines channel y2


    # Check input/output indices are within signal dimensions
    if (!isempty(filter(x->(x>S1.nu),in1))) error("S1: Input index exceeds input dimension."); end
    if (!isempty(filter(x->(x>S1.ny),out1))) error("S1: Output index exceeds output dimension."); end
    if (!isempty(filter(x->(x>S2.nu),in2))) error("S2: Input index exceeds input dimension."); end
    if (!isempty(filter(x->(x>S2.ny),out2))) error("S2: Output index exceeds output dimension."); end

    # Check compatibility of interconnection
    if (length(in1)!= length(in2)) error("Size of in1 and in2 must be same."); end
    if (length(out1)!= length(out2)) error("Size of out1 and out2 must be same."); end

    # Get indices for partitions of system matrices
    in1c = setdiff(1:S1.nu,in1); out1c = setdiff(1:S1.ny,out1); # These can be empty
    in2c = setdiff(1:S2.nu,in2); out2c = setdiff(1:S2.ny,out2); # These can be empty

    nv1 = length(in1c);
    nv2 = length(in2c);
    nz1 = length(out1c);
    nz2 = length(out2c);
    ns1 = S1.ns;
    ns2 = S2.ns;

    # Partitions of S1
    B1 = S1.B[:,in1c]; # defines v1 channel
    B2 = S1.B[:,in1]; # defines u1 channel
    C1 = S1.C[out1c,:];
    C2 = S1.C[out1,:];
    D11 = S1.D[out1c,in1c];
    D12 = S1.D[out1c,in1];
    D21 = S1.D[out1,in1c];
    D22 = S1.D[out1,in1];

    # Partitions of S2
    F1 = S2.B[:,in2]; # defines u2 channel
    F2 = S2.B[:,in2c];  # defines v2 channel
    G1 = S2.C[out2,:];
    G2 = S2.C[out2c,:];
    H11 = S2.D[out2,in2];
    H12 = S2.D[out2,in2c];
    H21 = S2.D[out2c,in2];
    H22 = S2.D[out2c,in2c];

    # Construct the new system matrices
    AA = [S1.A zeros(ns1,ns2);
    zeros(ns2,ns1) S2.A];

    # If nz1, nz2, nv1, nv2 are zero, this will result in empty matrices in the partitions.
    # Creating matrices with empty submatrices is a problem.
    # Using hcat and vcat to create the matrices, eliminates this problem.
    BB1 = hcat(B1,B2,zeros(ns1,nv2));
    BB2 = hcat(zeros(ns2,nv1),F1,F2);
    BB  = vcat(BB1,BB2);

    CC1 = hcat(C1,zeros(nz1,ns2));
    CC2 = hcat(C2,G1);
    CC3 = hcat(zeros(nz2,ns1),G2);
    CC  = vcat(CC1,CC2,CC3);

    DD1 = hcat(D11, D12, zeros(nz1,nv2));
    DD2 = hcat(D21, (D22+H11), H12);
    DD3 = hcat(zeros(nz2,nv1), H21, H22);
    DD = vcat(DD1,DD2,DD3);

    return StateSpace(AA,BB,CC,DD);
end


function Series(S1::StateSpace,S2::StateSpace, out1::Vector{Int64}, in2::Vector{Int64})::StateSpace
    #                                       +------+
    #                                v2 --->|      |
    #                       +------+        |  S2  |-----> y2
    #                       |      |------->|      |
    #              u1 ----->|      |y1   u2 +------+
    #                       |  S1  |
    #                       |      |---> z1
    #                       +------+
    # New system maps (v2;u1) -> (y2;z1)
    # -- out1 defines y1
    # -- in2 defines u2

    # Check input/output indices are within signal dimensions
    if (!isempty(filter(x->(x>S1.ny),out1))) error("S1: Output index exceeds output dimension."); end
    if (!isempty(filter(x->(x>S2.nu),in2)))  error("S2: Input index exceeds input dimension."); end

    # Check compatibility of interconnection
    if (length(out1)!= length(in2)) error("Incompatible I/O dimension for series connection."); end

    # Get indices for partitions of system matrices
    out1c = setdiff(1:S1.ny,out1); # Defines z, These can be empty
    in2c = setdiff(1:S2.nu,in2);   # These can be empty

    nv2 = length(in2c);
    nz1 = length(out1c);

    ns1 = S1.ns;
    ns2 = S2.ns;

    # Partitions of S1 -- follows the diagram
    B = S1.B;
    C1 = S1.C[out1,:]; # x1->y1
    C2 = S1.C[out1c,:]; # x1->z1
    D1 = S1.D[out1,:]; # u->y1
    D2 = S1.D[out1c,:]; # u->z1

    # Partitions of S2
    F1 = S2.B[:,in2c]; # defines v2 channel
    F2 = S2.B[:,in2];  # defines u2 channel
    G = S2.C; # y2
    H1 = S2.D[:,in2c]; # v2->y2
    H2 = S2.D[:,in2];  # u2->y2

    # Construct the new system matrices
    AA = [S1.A zeros(ns1,ns2);
          F2*C1 S2.A];

    # If nz1, nv2 are zero, this will result in empty matrices in the partitions.
    # Creating matrices with empty submatrices is a problem.
    # Using hcat and vcat to create the matrices, eliminates this problem.
    BB1 = hcat(zeros(ns1,nv2),B);
    BB2 = hcat(F1,F2*D1);
    BB  = vcat(BB1,BB2);

    CC1 = hcat(H2*C1,G);
    CC2 = hcat(C2,zeros(nz1,ns2));
    CC  = vcat(CC1,CC2);

    DD1 = hcat(H1,H2*D1);
    DD2 = hcat(zeros(nz1,nv2),D2);
    DD = vcat(DD1,DD2);

    return StateSpace(AA,BB,CC,DD);
end

function Feedback(S1::StateSpace,S2::StateSpace, feedin::Vector{Int64}, feedout::Vector{Int64},FeedbackSign::Int64=-1)::StateSpace

    #                            +------+
    #                v --------->|      |--------> z
    #                            |  S1  |
    #                u --->O---->|      |----+---> y
    #                      |  u1 +------+    |
    #                      |                 |
    #                      +-----[  S2  ]<---+
    # Maps  (v;u) -> (z;y)
    #   -- feedin defines u1
    #   -- feedout defines y

    # Check input/output indices are within signal dimensions
    if (!isempty(filter(x->(x>S1.nu),feedin))) error("S1: Feedback input index exceeds input dimension."); end
    if (!isempty(filter(x->(x>S1.ny),feedout))) error("S1: Feedback output index exceeds output dimension."); end

    # Check compatibility of interconnection
    if (length(feedin)!= S2.ny || length(feedout)!=S2.nu)
        error("Incompatible I/O dimensions for feedback.");
    end

    # Get indices for partitions of system matrices
    in1c = setdiff(1:S1.nu,feedin); # defines v, can be empty
    out1c = setdiff(1:S1.ny,feedout); # defines z, can be empty

    nv = length(in1c);
    nz = length(out1c);
    ns1 = S1.ns;
    ns2 = S2.ns;
    ny = length(feedout);
    nu = length(feedin);

    # Partitions of S1
    A = S1.A;
    B1 = S1.B[:,in1c]; # defines v channel
    B2 = S1.B[:,feedin]; # defines u1 channel
    C1 = S1.C[out1c,:]; # defines z channel
    C2 = S1.C[feedout,:]; # defines y channel
    D11 = S1.D[out1c,in1c]; # v->z
    D12 = S1.D[out1c,feedin]; # u1->z
    D21 = S1.D[feedout,in1c]; # v->y
    D22 = S1.D[feedout,feedin]; # u1->y

    # System matrices of S2
    E = S2.A;
    F = S2.B;
    G = S2.C;
    H = S2.D;

    I = eye(nu);
    K = (I-FeedbackSign*H*D22);
    if (det(K) == 0)
        error("Feedback interconnection is illposed.");
    end
    invK = inv(K);

    X = FeedbackSign*invK*hcat(H*C2,G);
    Y = invK*hcat(FeedbackSign*H*D21,I);

    # Construct the new system matrices
    AA1 = hcat(A,zeros(ns1,ns2)) + B2*X;
    AA2 = hcat(F*C2,E) + F*D22*X;
    AA = vcat(AA1,AA2);

    BB1 = hcat(B1,zeros(ns1,nu)) + B2*Y;
    BB2 = hcat(F*D21,zeros(ns2,nu)) + F*D22*Y;
    BB = vcat(BB1,BB2);

    CC1 = hcat(C1,zeros(nz,ns2)) + D12*X;
    CC2 = hcat(C2,zeros(ny,ns2)) + D22*X;
    CC = vcat(CC1,CC2);

    DD1 = hcat(D11,zeros(nz,nu)) + D12*Y;
    DD2 = hcat(D21,zeros(ny,nu)) + D22*Y;
    DD = vcat(DD1,DD2);

    return StateSpace(AA,BB,CC,DD);
end

# TO DO
# 1. Linear fractional transformation, and LFT algebra.
# LinearFractionalTransformation(...) -- uppper and lower.
