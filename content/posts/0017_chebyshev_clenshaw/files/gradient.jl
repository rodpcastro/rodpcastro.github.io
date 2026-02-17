"""
    gradient_clenshaw(f::ChebyshevSeries{T, N}, x::T) where {T, N} -> ChebyshevSeries{T, N-1}

Implements the Clenshaw algorithm to evaluate the series `f` and its partial derivative
with respect to the `N`-th dimension at a value `x` of its `N`-th dimension.

# Arguments
- `f::ChebyshevSeries{T, N}`: `N`-dimensional series to be evaluated
- `x::T`: Value of the `N`-th coordinate in the domain [-1, 1]

# Returns
- `ChebyshevSeries{T, N-1}`: `f` evaluated at `x`
- `ChebyshevSeries{T, N-1}`: partial derivative of `f` with
  respect to the `N`-th dimension evaluated at `x`
"""
function gradient_clenshaw(f::ChebyshevSeries{T, N}, x::T) where {T, N}
    a = f.coefs
    n = size(a, N)
    dx = 2x

    aₖ, aₙ₋₁, aₙ = (selectdim(a, N, i) for i in n-2:n)
    bₖ, bₖ₊₁ = (Array{T, N-1}(undef, a.size[1:N-1]) for _ in 1:2)
    cₖ, cₖ₊₁ = (Array{T, N-1}(undef, a.size[1:N-1]) for _ in 1:2)

    # bₖ used on the right-hand side actually represents bₖ₊₂.
    # bₖ₊₂ is ommited to reduce allocations. Idem for cₖ₊₂.
    
    # k = n - 2
    @. bₖ = aₙ  # Here, bₖ is bₖ₊₂
    @. bₖ₊₁ = aₙ₋₁ + dx*bₖ
    @. bₖ = aₖ + dx*bₖ₊₁ - bₖ
    bₖ, bₖ₊₁ = bₖ₊₁, bₖ
    
    # k = n-3 to 2
    @. cₖ = 2aₙ  # Here, cₖ is cₖ₊₂
    @. cₖ₊₁ = 2bₖ + dx*cₖ
    
    for k in n-3:-1:2
        aₖ = selectdim(a, N, k)
        @. bₖ = aₖ + dx*bₖ₊₁ - bₖ
        @. cₖ = 2bₖ₊₁ + dx*cₖ₊₁ - cₖ
        bₖ, bₖ₊₁ = bₖ₊₁, bₖ
        cₖ, cₖ₊₁ = cₖ₊₁, cₖ
    end

    # k = 1
    aₖ = selectdim(a, N, 1)
    @. bₖ = aₖ + x*bₖ₊₁ - bₖ
    @. cₖ = bₖ₊₁ + x*cₖ₊₁ - cₖ

    lbc = SVector(ntuple(i -> f.lb[i], Val(N-1)))
    ubc = SVector(ntuple(i -> f.ub[i], Val(N-1)))
    
    fc = ChebyshevSeries(bₖ, lbc, ubc)
    gc = ChebyshevSeries(cₖ, lbc, ubc)
    
    return fc, gc
end


function gradient_clenshaw(f::ChebyshevSeries{T, N}, x::SVector{N, T}) where {T, N}
    fc, gc = gradient_clenshaw(f, x[N])
    xc = SVector(ntuple(i -> x[i], Val(N-1)))
    return gradient_clenshaw(fc, xc)..., clenshaw(gc, xc)
end


function gradient_clenshaw(f::ChebyshevSeries{T, 1}, x::SVector{1, T}) where T
    return gradient_clenshaw(f, x[])
end


"""
    gradient(f::ChebyshevSeries{T, N}, x::SVector{N, T}) where {T, N} -> T, SVector{N, T}

Evaluates the series `f` and its gradient at a point `x`.

# Arguments
- `x::SVector{N, T}`: evaluation point

# Returns
- `T`: `f` evaluated at `x`
- `SVector{N, T}`: gradient of `f` evaluated at `x`
"""
function gradient(f::ChebyshevSeries{T, N}, x::SVector{N, T}) where {T, N}
    x̄ = normalize(f, x)
    dx̄_dx = @. 2.0 / (f.ub - f.lb)
    
    res = gradient_clenshaw(f, x̄)
    
    y = res[1]
    ∇y = SVector{N, T}(ntuple(i -> res[i+1], Val(N))) .* dx̄_dx
    
    return y, ∇y
end


function gradient(f::ChebyshevSeries{T, N}, x::AbstractVector{T}) where {T, N}
    return gradient(f, SVector{N, T}(x))
end


function gradient(f::ChebyshevSeries{T, 1}, x::T) where T
    res = gradient(f, SVector{1, T}(x))
    return res[1][], res[2][]
end


function gradient(g::ChebyshevCluster{T, N, M}, x::Union{AbstractVector{T}, T}) where {T, N, M}
    x_in_g, i = contains(g, x)
    
    if x_in_g
        return gradient(g.series[i], x)
    else
        throw(DomainError(x))
    end
end