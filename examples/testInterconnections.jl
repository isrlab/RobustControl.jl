# Test various interconnections using MATLAB

using MATLAB, LinearAlgebra
include("StateSpace.jl");
include("Interconnections.jl");

# Test Parallel(...)
println("\nTesting Parallel(...)");
for i=1:10
    A1 = rand(3,3);
    B1 = rand(3,3);
    C1 = rand(3,3);
    D1 = rand(3,3);

    A2 = rand(3,3);
    B2 = rand(3,3);
    C2 = rand(3,3);
    D2 = rand(3,3);

    s1 = StateSpace(A1,B1,C1,D1);
    s2 = StateSpace(A2,B2,C2,D2);
    s3 = Parallel(s1,s2,Vector(1:2),Vector(1:2),Vector(1:3),Vector(1:3));

    # Get MATLAB's implementation.
    mat"s1 = ss($A1,$B1,$C1,$D1)";
    mat"s2 = ss($A2,$B2,$C2,$D2)";
    mat"s3 = parallel(s1,s2,[1:2],[1:2],[1:3],[1:3])";

    mat"Ap=s3.A; Bp=s3.B;Cp=s3.C;Dp=s3.D"
    @mget Ap Bp Cp Dp

    eA = norm(s3.A-Ap);
    eB = norm(s3.B-Bp)
    eC = norm(s3.C-Cp)
    eD = norm(s3.D-Dp)

    println("Norm of error: ", norm([eA,eB,eC,eD]));
end

# Test Feedback(...)
println("\n\nTesting Feedback(...)");
for i=1:10
    A1 = rand(3,3);
    B1 = rand(3,3);
    C1 = rand(3,3);
    D1 = rand(3,3);

    A2 = rand(3,3);
    B2 = rand(3,2);
    C2 = rand(2,3);
    D2 = rand(2,2);

    s1 = StateSpace(A1,B1,C1,D1);
    s2 = StateSpace(A2,B2,C2,D2);
    s3 = Feedback(s1,s2,Vector(2:3),Vector(2:3));

    # Get MATLAB's implementation.
    mat"s1 = ss($A1,$B1,$C1,$D1)";
    mat"s2 = ss($A2,$B2,$C2,$D2)";
    mat"s3 = feedback(s1,s2,[2:3],[2:3])";
    mat"Ap=s3.A; Bp=s3.B;Cp=s3.C;Dp=s3.D"
    @mget Ap Bp Cp Dp

    eA = norm(s3.A-Ap);
    eB = norm(s3.B-Bp)
    eC = norm(s3.C-Cp)
    eD = norm(s3.D-Dp)

    println("Norm of error: ", norm([eA,eB,eC,eD]));
end

# Test Series(...)
println("\n\nTesting Series(...)");
for i=1:10
    A1 = rand(3,3);
    B1 = rand(3,3);
    C1 = rand(3,3);
    D1 = rand(3,3);

    A2 = rand(3,3);
    B2 = rand(3,3);
    C2 = rand(2,3);
    D2 = rand(2,3);

    s1 = StateSpace(A1,B1,C1,D1);
    s2 = StateSpace(A2,B2,C2,D2);
    s3 = Series(s1,s2,Vector(2:3),Vector(1:2));
    s33 = Subsystem(s3,Vector(2:4),Vector(1:2)); # Matlab returns this.

    # Get MATLAB's implementation.
    mat"s1 = ss($A1,$B1,$C1,$D1)";
    mat"s2 = ss($A2,$B2,$C2,$D2)";
    mat"s3 = series(s1,s2,[2:3],[1:2])";
    mat"Ap=s3.A; Bp=s3.B;Cp=s3.C;Dp=s3.D"
    @mget Ap Bp Cp Dp

    # The ordering of states is different in Matlab.
    # We check sorted elements of matrices
    eA = sort(s33.A[:])-sort(Ap[:]);
    eB = sort(s33.B[:])-sort(Bp[:]);
    eC = sort(s33.C[:])-sort(Cp[:]);
    eD = sort(s33.D[:])-sort(Dp[:]);

    println("Norm of error: ",norm([eA;eB;eC;eD]))
end
