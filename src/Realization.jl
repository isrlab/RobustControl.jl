# Various realizations of state-space.

function TransformSystem(S::StateSpace, T::Matrix)
    if IsSquareMatrix(T)
        invT = inv(T)
        A = invT*S.A*T
        B = invT*S.B
        C = S.C*T
        D = S.D
        return StateSpace(A,B,C,D)
    else
        error("Transformation matrix is not square")
    end
end


function ControllableCompanionForm(S::StateSpace)::StateSpace
    if S.nu == 1
        T = ControllabilityMatrix(S)
        return TransformSystem(S,T)
    else
        error("Input dimension is more than one. ")
    end
end

function ObservableCanonicalForm(S::StateSpace)::StateSpace
    if S.ny == 1
        T = ObservabilityMatrix(S)
        return TransformSystem(S,T)
    else
        error("Output dimension is more than one.")
    end
end

function ModalForm(S::StateSpace)::StateSpace
    eigval,eigvec = eigen(S.A)
    return TransformSystem(S,eigvec)
end

function KalmanDecomposition(S::StateSpace)::StateSpace
    # Todo
end

function BalancedForm(S::StateSpace)::StateSpace # Hankel form.
    # Todo
end

function MinimalRealization(S::StateSpace)::StateSpace
    # Todo
end

function SlowFastDecomposition(S::StateSpace)::StateSpace
    # Todo
end

function StableUnstableDecomposition(S::StateSpace)::StateSpace
    # Todo
end
