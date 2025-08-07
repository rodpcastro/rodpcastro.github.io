using SpecialFunctions: expintx

function Gᵈᵉᵉᵖ(P, Q, K)
    # Wave Green function for deep waters. This function
    # also returns the gradient and Hessian matrix.
    ξ, ζ = P
    x, z = Q

    x1 = x - ξ
    z1 = z - ζ
    z3 = z + ζ

    R = abs(x1)
    v₁ = abs(z1)
    v₃ = abs(z3)

    # Auxilary variables d.
    d1 = R^2
    d2 = v₁^2
    d3 = v₃^2
    d4 = d1 + d2
    d5  = d4^2
    d6 = d1 + d3
    d7 = d6^2
    d8 = 2R

    r₁ = √d4
    r₃ = √d6

    X = K*R
    V₃ = K*v₃
    Z = V₃ - im*X

    # Auxilary variables k.
    k1 = 2K
    k2 = k1*K

    # Auxilary variables e.
    e1 = exp(-Z)
    e2 = expintx(-Z)
    e3 = e2 + 1/Z
    e4 = e3 + 1/Z^2
    e5 = 2π*e1
    e6 = K*e5
    e7 = K*e6

    # Auxilary variables s.
    sx1 = sign(x1)
    sx2 = sx1^2
    sz1 = sign(z1)
    sz2 = sz1^2
    sz3 = sign(z3)
    sz4 = sz3^2
    sw1 = sx1 * sz1
    sw2 = sx1 * sz3

    if X ≤ 1
        # Near field
        G = log(K*r₁) + log(K*r₃) - 2 * (real(e2) + log(abs(Z))) - im*e5
    else
        # Far field
        G = log(r₁/r₃) - 2*real(e2) - im*e5
    end

    Gx = sx1 * (R/d4 - R/d6 + k1*imag(e3) + e6)
    Gz = sz1 * v₁/d4 + sz3 * (-v₃/d6 + k1*real(e3) + im*e6)
    ∇G = [Gx, Gz]

    H1 = (d2 - d1)/d5
    H2 = (d1 - d3)/d7 + k2*real(e4) + im*e7
    Gxx = sx2 * (H1 + H2)
    Gzz = -sz2 * H1 - sz4 * H2
    Gxz = -sw1 * d8*v₁/d5 + sw2 * (d8*v₃/d7 - k2*imag(e4) - e7)
    ∇²G = [Gxx Gxz; Gxz Gzz]

    return (G, ∇G, ∇²G)
end