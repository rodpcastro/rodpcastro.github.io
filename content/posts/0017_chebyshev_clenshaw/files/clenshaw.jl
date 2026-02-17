"""
    clenshaw(f::ChebyshevSeries{T, N}, x::T) where {T, N} -> ChebyshevSeries{T, N-1}

Implements the Clenshaw algorithm to evaluate the
series `f` at a value `x` of its `N`-th dimension.

# Arguments
- `f::ChebyshevSeries{T, N}`: `N`-dimensional series to be evaluated
- `x::T`: Value of the `N`-th coordinate in the domain [-1, 1]

# Returns
- `ChebyshevSeries{T, N-1}`: `f` evaluated at `x`
"""
function clenshaw(f::ChebyshevSeries{T, N}, x::T) where {T, N}
    a = f.coefs
    n = size(a, N)
    dx = 2x
    
    aₙ₋₁, aₙ = (selectdim(a, N, i) for i in n-1:n)
    bₖ, bₖ₊₁ = (Array{T, N-1}(undef, a.size[1:N-1]) for _ in 1:2)

    # bₖ used on the right-hand side actually represents bₖ₊₂.
    # bₖ₊₂ is ommited to reduce allocations.
    
    # k = n-2 to 2. 
    @. bₖ = aₙ  # Here, bₖ is bₖ₊₂
    @. bₖ₊₁ = aₙ₋₁ + dx*bₖ
    
    for k in n-2:-1:2
        aₖ = selectdim(a, N, k)
        @. bₖ = aₖ + dx*bₖ₊₁ - bₖ
        bₖ, bₖ₊₁ = bₖ₊₁, bₖ
    end

    # k = 1
    aₖ = selectdim(a, N, 1)
    @. bₖ = aₖ + x*bₖ₊₁ - bₖ

    lbc = SVector(ntuple(i -> f.lb[i], Val(N-1)))
    ubc = SVector(ntuple(i -> f.ub[i], Val(N-1)))
    
    fc = ChebyshevSeries(bₖ, lbc, ubc)
    
    return fc
end


function clenshaw(f::ChebyshevSeries{T, N}, x::SVector{N, T}) where {T, N}
    fc = clenshaw(f, x[N])
    xc = SVector(ntuple(i -> x[i], Val(N-1)))
    return clenshaw(fc, xc)
end


function clenshaw(f::ChebyshevSeries{T, 1}, x::SVector{1, T}) where T
    return clenshaw(f, x[])
end


"""
    (f::ChebyshevSeries{T, N})(x::SVector{N, T}) where {T, N} -> T

Evaluates the series `f` at a point `x`.

# Arguments
- `x::SVector{N, T}`: evaluation point

# Returns
- `T`: `f` evaluated at `x`
"""
function (f::ChebyshevSeries{T, N})(x::SVector{N, T}) where {T, N}
    x̄ = normalize(f, x)
    y = clenshaw(f, x̄)
    return y
end


function (f::ChebyshevSeries{T, N})(x::AbstractVector{T}) where {T, N}
    return f(SVector{N, T}(x))
end


function (f::ChebyshevSeries{T, 1})(x::T) where T
    return f(SVector{1, T}(x))[]
end


function (g::ChebyshevCluster{T, N, M})(x::Union{AbstractVector{T}, T}) where {T, N, M}
    x_in_g, i = contains(g, x)
    
    if x_in_g
        return g.series[i](x)
    else
        throw(DomainError(x))
    end
end