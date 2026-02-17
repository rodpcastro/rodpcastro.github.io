---
date: '2026-02-16'
draft: false
title: 'Evaluating Chebyshev series and its derivatives with the Clenshaw algorithm'
author: 'Rodrigo Castro'
summary: 'A three-dimensional Chebyshev series, its gradient and Hessian are evaluated with the Clenshaw algorithm.'
tags: ['Chebyshev', 'Approximation Theory', 'Julia']
---

## Introduction
This post is an improvement on the [last article][chebyl1] about the use of [Chebyshev Polynomials][chebypoly] to approximate functions. This time, the Chebyshev series is evaluated with the Clenshaw algorithm, and a extension of it is used to compute the first and second order derivatives.

## Methods
First we define the function, followed by its Chebyshev series representation and the algorithm used to evaluate the series and its derivatives.

### Function

Once again, the function we aim to approximate is the following improper integral

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

The domain of interest is given by the points $A \in [0, 0.5]$, $B \in [0, 1]$ and $\log_{10} H \in [-2, 2]$. The integral $(1)$ is computed by contour integration along the path that connects the points $\{0, H+i, H+1, \infty\}$.

The following [Julia] code implements the computation of $L_1$, its gradient and Hessian using the [Integrals.jl] library.

{{< dropdown_file title="$L_1$ and derivatives computed with Integrals.jl" src="l1.jl" fmt="julia" >}}

### Chebyshev series
Using [Chebyshev polynomials][chebypoly] of the first kind, $L_1$ can be approximated by the following series

$$\eq{
L_1(x, y, z) \simeq \sum_{p=0}^{n_p} \sum_{q=0}^{n_q} \sum_{r=0}^{n_r} a_{pqr} T_p(x) T_q(y) T_r(z),
}$$

where $(x, y, z) \in [-1, 1]^3$ is the transformed version of the variables $(A, B, \log_{10} H)$ to satisfy the standard domain of Chebyshev polynomials. The coefficients $a_{pqr}$ are computed with the [FastChebInterp.jl] library, as described in a [previous post][chebyl1].

The range of $\log_{10} H$ is subdvided in three parts: $[-2, 0.15]$, $[0.15, 1]$ and $[1, 2]$. The three subdomains have the same number of coefficients $(n_p, n_q, n_r) = (16, 20, 41)$ and are organized in the following composite type and saved to disk and loaded from it with the [JLD2.jl] library. The three subdomains are named `L₁1`, `L₁2` and `L₁3`, respectively.

```julia
using StaticArrays
using JLD2

struct ChebyshevSeries{T, N}
    coefs::Array{T, N}
    lb::SVector{N, T}  # domain lower bound
    ub::SVector{N, T}  # domain upper bound
end

@load "l1domains.jld2" L₁1 L₁2 L₁3
```

### Clenshaw algorithm
The [Clenshaw algorithm][clenshaw], as described by *Clenshaw* (1955), is a recursive method to evaluate a linear combination of Chebyshev polynomials, and is defined below.

Consider the truncated Chebyshev series

$$\eq{
s_n = \sum_{k=0}^{n} a_k T_k(x).
}$$

Now, we define the recurrence

$$
\eq{
\begin{split}
& b_n = a_n, \quad b_{n-1} = a_{n-1} + 2 x b_n, \\
& b_k = a_k + 2 x b_{k+1} - b_{k+2}, \quad k = n-2, \ldots, 1
\end{split}
}$$

Then, $s_n$ is given by

$$\eq{
s_n = a_0 + x b_1 - b_2.
}$$

An extension of this algorithm was obtained by *Skrzipek* (1998) to get the derivatives of any order. For the first and second order derivatives, we define the recurrence relations

$$
\eq{
\begin{split}
& c_{n-1} = 2 b_n, \quad c_{n-2} = 2 b_{n-1} + 2 x c_{n-1}, \\
& c_k = 2 b_{k+1} + 2 x c_{k+1} - c_{k+2}, \quad k = n-3, \ldots, 1,
\end{split}
}$$

$$
\eq{
\begin{split}
& d_{n-2} = 2 c_{n-1}, \quad d_{n-3} = 2 c_{n-2} + 2 x d_{n-2}, \\
& d_k = 2 c_{k+1} + 2 x d_{k+1} - d_{k+2}, \quad k = n-4, \ldots, 1,
\end{split}
}$$

Then, the first and second order derivatives are

$$\eq{
\frac{d s_n}{d x} = b_1 + x c_1 - c_2,
}$$

$$\eq{
\frac{d^2 s_n}{d x^2} = 2(c_1 + x d_1 - d_2).
}$$

The algorithm is implemented in the following Julia function. Notice that, in Julia, the first index is 1, not 0 like in the series defined in $(4)$.

```julia
function hessian_clenshaw(f::ChebyshevSeries{T, N}, x::T) where {T, N}
    a = f.coefs
    n = size(a, N)
    dx = 2x

    aₖ, aₙ₋₁, aₙ = (selectdim(a, N, i) for i in n-2:n)
    bₖ, bₖ₊₁ = (Array{T, N-1}(undef, a.size[1:N-1]) for _ in 1:2)
    cₖ, cₖ₊₁ = (Array{T, N-1}(undef, a.size[1:N-1]) for _ in 1:2)
    dₖ, dₖ₊₁ = (Array{T, N-1}(undef, a.size[1:N-1]) for _ in 1:2)

    # bₖ used on the right-hand side actually represents bₖ₊₂.
    # bₖ₊₂ is ommited to reduce allocations. Idem for cₖ₊₂ and dₖ₊₂.
    
    # k = n - 2
    @. bₖ = aₙ  # Here, bₖ is bₖ₊₂
    @. bₖ₊₁ = aₙ₋₁ + dx*bₖ
    @. bₖ = aₖ + dx*bₖ₊₁ - bₖ
    bₖ, bₖ₊₁ = bₖ₊₁, bₖ
    
    # k = n - 3
    @. cₖ = 2aₙ  # Here, cₖ is cₖ₊₂
    @. cₖ₊₁ = 2bₖ + dx*cₖ
    
    aₖ = selectdim(a, N, n-3)
    @. bₖ = aₖ + dx*bₖ₊₁ - bₖ
    @. cₖ = 2bₖ₊₁ + dx*cₖ₊₁ - cₖ
    bₖ, bₖ₊₁ = bₖ₊₁, bₖ
    cₖ, cₖ₊₁ = cₖ₊₁, cₖ
    
    # k = n-4 to 2
    @. dₖ = 4aₙ  # Here, dₖ is dₖ₊₂
    @. dₖ₊₁ = 2cₖ + dx*dₖ
    
    for k in n-4:-1:2
        aₖ = selectdim(a, N, k)
        @. bₖ = aₖ + dx*bₖ₊₁ - bₖ
        @. cₖ = 2bₖ₊₁ + dx*cₖ₊₁ - cₖ
        @. dₖ = 2cₖ₊₁ + dx*dₖ₊₁ - dₖ
        bₖ, bₖ₊₁ = bₖ₊₁, bₖ
        cₖ, cₖ₊₁ = cₖ₊₁, cₖ
        dₖ, dₖ₊₁ = dₖ₊₁, dₖ
    end

    # k = 1
    aₖ = selectdim(a, N, 1)
    @. bₖ = aₖ + x*bₖ₊₁ - bₖ
    @. cₖ = bₖ₊₁ + x*cₖ₊₁ - cₖ
    @. dₖ = 2.0*(cₖ₊₁ + x*dₖ₊₁ - dₖ)

    lbc = SVector(ntuple(i -> f.lb[i], Val(N-1)))
    ubc = SVector(ntuple(i -> f.ub[i], Val(N-1)))
    
    fc = ChebyshevSeries(bₖ, lbc, ubc)
    gc = ChebyshevSeries(cₖ, lbc, ubc)
    hc = ChebyshevSeries(dₖ, lbc, ubc)
    
    return fc, gc, hc
end
```

The function above handles Chebyshev series of any dimension and returns the series and its derivatives with respect to the last dimension evaluated at a value `x` of its last dimension. So, for the three-dimensional case of $L_1$ defined in $(3)$, evaluating this function recursively for $z$, $y$ and $x$ provides the values of $L_1$, its gradient and Hessian matrix at a point $(x, y, z) \in [-1, 1]^3$.

See the [source code](#appendices) for more details on how this is done.

## Results
In the following `L₁c` represents $L_1$ as computed by the Chebyshev series $(3)$, and `∇L₁c` and `HL₁c` are its gradient and Hessian matrix, respectively. The function $L_1$ was computed for $10^5$ points randomly distributed in the domain of interest, first with the [function that uses numerical integration](#function), then with the [Clenshaw algorithm](#clenshaw-algorithm) described previously. The output below shows the mean absolute error between the two.

```console
L₁c mean absolute error = 1e-15

∇L₁c mean absolute error = [3e-14  1e-14  2e-14]

HL₁c mean absolute error = [3e-12  7e-13  9e-13]
                           [7e-13  1e-12  7e-13]
                           [9e-13  7e-13  3e-12]

```

Statistical quantities are used here because, as explained in the [previous post][chebyl1], `Integrals.jl` does not always converge to the desired accuracy. Fortunately, it does converge for most part of the points in the domain of interest, as demonstrated by the next results.

For `L₁c`, `∇L₁c` and `HL₁c`, the error value for each point is compared with a corresponding reference value. The reference value is the mean absolute error plus three times the standard deviation. The following result shows that for 98-99% of the $10^5$ points, the error between the integral expression and the Chebyshev series is lower than the reference value.

```console
Percentage of L₁c error values lower than the reference = 99.2%

Percentage of ∇L₁c error values lower than the reference = [98.6%  98.6%  99.5%]

Percentage of HL₁c error values lower than the reference = [98.8%  98.9%  99.5%]
                                                           [98.9%  99.0%  99.1%]
                                                           [99.5%  99.1%  99.7%]
```

## Conclusion
The Clenshaw algorithm was successfully implemented and a nice framework for working with Chebyshev series was established. The next step is using this base code to make developments in other projects, like the free-surface Green function in two and three dimensions.

## References
1. C. W. Clenshaw. 1955. A note on the summation of Chebyshev series. Math. Comp. 9 (July 1955), 118–120. https://doi.org/10.1090/S0025-5718-1955-0071856-0
2. M. R. Skrzipek. 1998. Polynomial evaluation and associated polynomials. Numer. Math. 79, 4 (June 1998), 601–613. https://doi.org/10.1007/s002110050354

## Appendices
* {{< post_files_view >}}
* {{< post_files_download >}}

<!--Links-->
[chebyl1]: ../0015_chebyshev_l1/
[Julia]: https://julialang.org/
[chebypoly]: https://en.wikipedia.org/wiki/Chebyshev_polynomials
[Integrals.jl]: https://docs.sciml.ai/Integrals/stable/
[FastChebInterp.jl]: https://github.com/JuliaMath/FastChebInterp.jl
[JLD2.jl]: https://juliaio.github.io/JLD2.jl/
[StaticArrays.jl]: https://github.com/JuliaArrays/StaticArrays.jl
[clenshaw]: https://en.wikipedia.org/wiki/Clenshaw_algorithm
