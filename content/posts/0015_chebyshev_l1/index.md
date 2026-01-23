---
date: '2025-12-28'
draft: false
title: 'Approximating functions with Chebyshev polynomials'
author: 'Rodrigo Castro'
summary: 'An improper integral with three parameters is approximated to 14 digits of decimal precision with Chebyshev polynomials.'
tags: ['Chebyshev', 'Approximation Theory', 'Julia']
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

The numerical computation of this integral can be performed in [Julia] with the [Integrals.jl] library, with absolute and relative tolerances of $10^{-16}$ to achieve double precision results, as displayed below:

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

Using this base of polynomials, $L_1$ can be approximated by the following series

$$\eq{
L_1(x, y, z) = \sum_{p=0}^{\infty} \sum_{q=0}^{\infty} \sum_{r=0}^{\infty} a_{pqr} T_p(x) T_q(y) T_r(z),
}$$

where $(x, y, z) \in [-1, 1]^3$ are the transformed version of the variables $(A, B, \log_{10} H)$ to satisfy the standard domain of Chebyshev polynomials.

Since $T_n(x) \le 1$ for $x \in [-1, 1]$, the series $(4)$ can be truncated for absolute values of the coefficients $a_{pqr}$ lower than the intended precision. To calculate these coefficients, the Julia library [FastChebInterp.jl] was used.

`FastChebInterp` provides two options to obtain the coefficients $a_{pqr}$: interpolation with `chebinterp` and fitting with `chebregression`. The interpolation option was tried first, but crashed the computer before getting close to the desired precision, probably due to the inconsistent results obtained with `QuadGKJL`, as explained earlier. Fitting presented itself as the best option, because it smooths out these deviations, and never crashed my computer.

To reduce the total number of coefficients and, consequently, the computational cost of expression $(4)$, the range of $\log_{10} H$ was subdvided in three parts: $[-2, 0.15]$, $[0.15, 1]$ and $[1, 2]$. The following code computes the coefficients for the first part of the domain, and the same process was followed for the other two parts.

```julia
using FastChebInterp  # version 1.2.0

A1, A2 = 0.0, 0.5
B1, B2 = 0.0, 1.0
logH1, logH2 = -2.0, 0.15

# Number of points for regression.
np, nq, nr = 15, 19, 41

lb, ub = [A1, B1, logH1], [A2, B2, logH2]
x = vec(chebpoints((np, nq, nr), lb, ub))
c = chebregression(x, L₁.(x), lb, ub, (np, nq, nr))

# The coefficients are stored in the coefs attribute.
# coefs.size returns (16, 20, 42)
coefs = c.coefs;

# I stopped increasing np, nq and nr when the regression computed
# coefficients lower than machine epsilon for double precision.
# maximum(coefs[end, :, :]) returns 3.820990154979531e-17
# maximum(coefs[:, end, :]) returns 2.1863957579554176e-16
# maximum(coefs[:, :, end]) returns 2.2126464867121333e-16
```

The size of the coefficients arrays in each domain are $(16, 20, 42)$, $(17, 21, 42)$ and $(17, 21, 40)$, respectively. The coefficients can then be stored in a file and expression $(4)$ can be evaluated in Julia in the following way for the first domain, and similarly for the other two.

```julia
function L₁_cheby(A::Float64, B::Float64, logH::Float64)
    # The coefficients are defined in the global scope.

    A1, A2 = 0.0, 0.5
    B1, B2 = 0.0, 1.0
    logH1, logH2 = -2.0, 0.15  # Domain 1

    x = (2A - A2 - A1) / (A2 - A1)
    y = (2B - B2 - B1) / (B2 - B1)
    z = (2logH - logH2 - logH1) / (logH2 - logH1)

    Tx = reshape([T(p, x) for p in 0:size(coefs, 1)-1], :, 1, 1)
    Ty = reshape([T(q, y) for q in 0:size(coefs, 2)-1], 1, :, 1)
    Tz = reshape([T(r, z) for r in 0:size(coefs, 3)-1], 1, 1, :)

    return sum(coefs .* Tx .* Ty .* Tz)
end
```

## Results
The three figures below present the maximum absolute value of the coefficients along each axis of the coefficients array, which indicates the total number of coefficients necessary in the variables $A$, $B$ and $\log_{10} H$ to approximate $L_1$ with high level of precision.

{{< figure src="images/max_abs_coefs1.svg" alt="Coefs Domain 1" align="center" >}}

{{< figure src="images/max_abs_coefs2.svg" alt="Coefs Domain 2" align="center" >}}

{{< figure src="images/max_abs_coefs3.svg" alt="Coefs Domain 3" align="center" >}}

The Chebyshev approximation was tested against the original expression $(1)$ for 30 points along $A$, 40 points along $B$ and 80 points along $\log_{10} H$, providing a total of 96000 test points. The mean absolute residue was of $1.07 \times 10^{-15}$, with a standard deviation of $8.78 \times 10^{-16}$. Also, 99.99% of the residual data is smaller than $10^{-14}$, making the approximation function precise to 14 decimal digits. In addition to being very precise, the approximated version is, on average, 13 times faster than the original one.

## Conclusion
Approximation with Chebyshev polynomials proved to be very precise and efficient. However, the algorithm presented in this post can be further improved in speed. The coefficients arrays still have several elements smaller than machine precision, making their use unnecessary. For example, 43% of the coefficients in domain 1 are smaller than $10^{-17}$, and similar numbers are found in the other two domains.

## References
1. Ed Mackay. 2021. The Green function for diffraction and radiation of regular waves by two-dimensional structures. European Journal of Mechanics - B/Fluids 87 (May 2021), 151–160. https://doi.org/10.1016/j.euromechflu.2021.01.012
2. Steven Johnson. 2024. FastChebInterp.jl. https://github.com/JuliaMath/FastChebInterp.jl.

## Appendices
* {{< post_files_view >}}
* {{< post_files_download >}}

<!--Links-->
[Julia]: https://julialang.org/
[2d_fsgf]: ../0009_2d_inf_depth_fsurface_gfunction/
[cheby]: https://en.wikipedia.org/wiki/Chebyshev_polynomials
[Integrals.jl]: https://docs.sciml.ai/Integrals/stable/
[FastChebInterp.jl]: https://github.com/JuliaMath/FastChebInterp.jl
