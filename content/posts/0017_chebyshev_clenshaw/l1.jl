using Integrals

rtol = 1e-16
atol = 1e-16
imax = 1e6

function L₁(x::AbstractVector{<:Real})
    A, B, logH = x
    H = 10^logH
    hd = H * log(10.0)
    
    f(u) = (u + H) / (u - H - (u + H)*exp(-2u))
    f3(u) = 2 / (u - H - (u + H)*exp(-2u))^2
    f33(u) = 4*(1 + exp(-2u)) / (u - H - (u + H)*exp(-2u))^3
    
    g(u) = exp(-u*(2+B)) + exp(-u*(2-B))
    g2(u) = exp(-u*(2-B)) - exp(-u*(2+B))
    g22(u) = exp(-u*(2+B)) + exp(-u*(2-B))
    
    h(u) = cos(u*A)
    h1(u) = -sin(u*A)
    h11(u) = -cos(u*A)
    
    p(u, _) = (f(u) * g(u) * h(u) + exp(-u)) / u
    
    p1(u, _) = f(u) * g(u) * h1(u)
    p2(u, _) = f(u) * g2(u) * h(u)
    p3(u, _) = f3(u) * g(u) * h(u)

    p11(u, _) = f(u) * g(u) * h11(u) * u
    p12(u, _) = f(u) * g2(u) * h1(u) * u
    p13(u, _) = f3(u) * g(u) * h1(u) * u

    p21(u, _) = f(u) * g2(u) * h1(u) * u
    p22(u, _) = f(u) * g22(u) * h(u) * u
    p23(u, _) = f3(u) * g2(u) * h(u) * u

    p31(u, _) = f3(u) * g(u) * h1(u) * u
    p32(u, _) = f3(u) * g2(u) * h(u) * u
    p33(u, _) = f33(u) * g(u) * h(u)
    
    paths = [(0.0, H+im), (H+im, H+1.0), (H+1.0, Inf)]
    M = 0.0
    M1, M2, M3 = 0.0, 0.0, 0.0
    M11, M12, M13 = 0.0, 0.0, 0.0
    M21, M22, M23 = 0.0, 0.0, 0.0
    M31, M32, M33 = 0.0, 0.0, 0.0
    for path in paths
        IP = IntegralProblem(p, path)
        
        IP1 = IntegralProblem(p1, path)
        IP2 = IntegralProblem(p2, path)
        IP3 = IntegralProblem(p3, path)

        IP11 = IntegralProblem(p11, path)
        IP12 = IntegralProblem(p12, path)
        IP13 = IntegralProblem(p13, path)

        IP21 = IntegralProblem(p21, path)
        IP22 = IntegralProblem(p22, path)
        IP23 = IntegralProblem(p23, path)

        IP31 = IntegralProblem(p31, path)
        IP32 = IntegralProblem(p32, path)
        IP33 = IntegralProblem(p33, path)
        
        M += solve(IP, QuadGKJL(); reltol=rtol, abstol=atol, maxiters=imax).u
        
        M1 += solve(IP1, QuadGKJL(); reltol=rtol, abstol=atol, maxiters=imax).u
        M2 += solve(IP2, QuadGKJL(); reltol=rtol, abstol=atol, maxiters=imax).u
        M3 += solve(IP3, QuadGKJL(); reltol=rtol, abstol=atol, maxiters=imax).u

        M11 += solve(IP11, QuadGKJL(); reltol=rtol, abstol=atol, maxiters=imax).u
        M12 += solve(IP12, QuadGKJL(); reltol=rtol, abstol=atol, maxiters=imax).u
        M13 += solve(IP13, QuadGKJL(); reltol=rtol, abstol=atol, maxiters=imax).u

        M21 += solve(IP21, QuadGKJL(); reltol=rtol, abstol=atol, maxiters=imax).u
        M22 += solve(IP22, QuadGKJL(); reltol=rtol, abstol=atol, maxiters=imax).u
        M23 += solve(IP23, QuadGKJL(); reltol=rtol, abstol=atol, maxiters=imax).u

        M31 += solve(IP31, QuadGKJL(); reltol=rtol, abstol=atol, maxiters=imax).u
        M32 += solve(IP32, QuadGKJL(); reltol=rtol, abstol=atol, maxiters=imax).u
        M33 += solve(IP33, QuadGKJL(); reltol=rtol, abstol=atol, maxiters=imax).u
    end

    M = real(M)
    GM = [real(M1), real(M2), real(M3)*hd]
    HM = [
        real(M11)    real(M12)    real(M13)*hd;
        real(M12)    real(M22)    real(M23)*hd;
        real(M13)*hd real(M23)*hd (real(M33)*hd + real(M3)*log(10))*hd;
    ]
    
    return M, GM, HM
end