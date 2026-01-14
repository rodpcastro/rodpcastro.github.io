function T(n::Int, x::Float64)
    if n == 0
        return 1.0
    elseif n == 1
        return x
    end

    t0 = 1.0
    t1 = x
    for i in 2:n
        t2 = 2x * t1 - t0
        t0 = t1
        t1 = t2
    end
    return t1
end


function L₁_cheby(A::Float64, B::Float64, logH::Float64)

    if A < 0.0 || A > 0.5
      throw(DomainError(A, "A must range from 0.0 to 0.5"))
    end

    if B < 0.0 || B > 1.0 
      throw(DomainError(B, "B must range from 0.0 to 1.0"))
    end

    if abs(logH) > 2.0
      throw(DomainError(logH, "log(H) must range from -2.0 to 2.0"))
    elseif logH < 0.15
      coefs = coefs1
      logH1, logH2 = -2.0, 0.15
    elseif logH < 1.0
      coefs = coefs2
      logH1, logH2 = 0.15, 1.0
    else
      coefs = coefs3
      logH1, logH2 = 1.0, 2.0
    end

    A1, A2 = 0.0, 0.5
    B1, B2 = 0.0, 1.0

    x = (2A - A2 - A1) / (A2 - A1)
    y = (2B - B2 - B1) / (B2 - B1)
    z = (2logH - logH2 - logH1) / (logH2 - logH1)

    Tx = reshape([T(p, x) for p in 0:size(coefs, 1)-1], :, 1, 1)
    Ty = reshape([T(q, y) for q in 0:size(coefs, 2)-1], 1, :, 1)
    Tz = reshape([T(r, z) for r in 0:size(coefs, 3)-1], 1, 1, :)

    return sum(coefs .* Tx .* Ty .* Tz)
end