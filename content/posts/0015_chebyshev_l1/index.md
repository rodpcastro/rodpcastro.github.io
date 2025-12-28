---
date: '2025-12-27'
draft: true
title: 'Approximating functions with Chebyshev polynomials'
author: 'Rodrigo Castro'
summary: 'An improper integral with three parameters is approximated to 14 digits of precision with Chebyshev polynomials.'
tags: ['Chebyshev', 'Approximation Theory', Julia']
---

## Introduction
This article is a development on the two-dimensional free-surface Green function. The infinite-depth case was already discussed in an [earlier post][2d_fsgf], while this text is concerned with the finite-depth one, focusing on the efficient evaluation of a complicated component of its real part by means of Chebyshev polynomials.

## Methods
The function we aim to approximate is the following improper integral evaluated in the principal value sense

$$\eq{
\begin{split}
& L_1(A, B, H) = \\
& \mathrm{p.v.} \int_{0}^{\infty} \left[ 
f(u, H) \left( e^{-u(2+B)} + e^{-u(2-B)} \right) \cos{(u A)} +  e^{-u} \right] \, \frac{du}{u},
\end{split}
}$$

where $f$ is defined by

$$\eq{
f(u, H) = \frac{u + H}{(u - H) - (u + H)e^{-2u}}.
}$$

$A \in [0, 0.5]$, $B \in [0, 1]$ and $H \in [10^{-2}, 10^{2}]$ are dimensionless quantities. $A$ and $B$ are, respectively, the dimensionless horizontal and vertical distances between the source and field points, while $H$ is the dimensionless depth.

The integral in $(1)$ must be evaluated by contour integration to avoid numerical issues close to the poles of $f$. The path of integration suggested by *Mackay* (2021) is the line connecting the points $\{0, H+i, H+1, \infty\}$.

The numerical computation of this integral can be performed in [Julia] with the [Integrals.jl] library, with absolute and relative tolerances of $10^{-16}$ to achieve machine double precision, as displayed below:

```julia
using Integrals

function L₁(x::Vector{Float64})
    A, B, logH = x
    H = 10^logH

    f(u) = (u + H) / (u - H - (u + H)*exp(-2u))
    g(u) = (exp(-u*(2+B)) + exp(-u*(2-B))) * cos(u*A)
    h(u, _) = (f(u) * g(u) + exp(-u)) / u
    
    paths = [(0.0, H+im), (H+im, H+1.0), (H+1.0, Inf)]
    M = 0.0
    for path in paths
        IP = IntegralProblem(h, path)
        M += solve(IP, QuadGKJL(); reltol=1e-16, abstol=1e-16).u
    end
    return real(M)
end
```

To prepare for the Chebyshev approximation, the arguments of $L_1$ are placed in a vector `x` containing three values, corresponding to $A$, $B$ and $\log_{10} H$, respectively. The $\log_{10} H$ is choosen as input instead of $H$ because it makes the function variation smaller. Notice that $\log_{10} H \in [-2, 2]$.

An **important observation** to make at this point is that the numerical evaluation of $L_1$, as provided by the code above, does not always guarantee the desired precision. For unknown reasons, the results obtained with the `QuadGKJL` algorithm present deviations from the smooth behavior for a small number of test points, making the actual residue as bad as $10^{-10}$. Despite this inconvenience, the major part of the test points behave smoothly, therefore, the residual analysis employed to validate the polynomial approximation will be based on statistical quantities.

### Chebyshev polynomial fitting
In one variable $x \in [-1, 1]$ (not the same $x$ defined for $L_1$), the Chebyshev polynomials of the first kind are defined by the following recurrence relation

$$\eq{
T_0(x) = 1, \quad T_1(x) = x, \quad T_{n} = 2x T_{n-1}(x) - T_{n-2}(x),
}$$

which can be implemented in Julia as follows.

```julia
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
```




```julia
function L₁_cheby(coefs::Array{Float64, 3}, A::Float64, B::Float64, logH::Float64)
    x = (2A - A2 - A1) / (A2 - A1)
    y = (2B - B2 - B1) / (B2 - B1)
    z = (2logH - logH2 - logH1) / (logH2 - logH1)

    Tx = reshape([T(p, x) for p in 0:size(coefs, 1)-1], :, 1, 1)
    Ty = reshape([T(q, y) for q in 0:size(coefs, 2)-1], 1, :, 1)
    Tz = reshape([T(r, z) for r in 0:size(coefs, 3)-1], 1, 1, :)
    
    sum(coefs .* Tx .* Ty .* Tz)
end
```

To approximate $L_1$ with a base of [Chebyshev polynomials][cheby] of the first kind, the Julia library [FastChebInterp.jl] was used. This library provides two options to obtain an approximation function: interpolation with `chebinterp` and fitting with `chebregression`. The interpolation option was tried first, but crashed the computer before getting close to the desired precision, probably because of the inconsistent results obtained with `QuadGKJL` as explained earlier. Fitting presented itself as the best option, because it smooths out these inconsistencies, and never crashed my computer.

## Results


## Conclusion


## References
1. Ed Mackay. 2021. The Green function for diffraction and radiation of regular waves by two-dimensional structures. European Journal of Mechanics - B/Fluids 87 (May 2021), 151–160. https://doi.org/10.1016/j.euromechflu.2021.01.012

## Appendices
* {{< post_files_view >}}
* {{< post_files_download >}}

<!--Links-->
[Julia]: https://julialang.org/
[2d_fsgf]: ../0009_2d_inf_depth_fsurface_gfunction/
[cheby]: https://en.wikipedia.org/wiki/Chebyshev_polynomials
[Integrals.jl]: https://docs.sciml.ai/Integrals/stable/
[FastChebInterp.jl]: https://github.com/JuliaMath/FastChebInterp.jl
