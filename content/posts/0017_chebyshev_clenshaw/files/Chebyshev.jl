module Chebyshev

export ChebyshevSeries, ChebyshevCluster, gradient, hessian

using StaticArrays


"""
    ChebyshevSeries{T, N}

The Chebyshev Series approximation of a `N`-dimensional
function defined in a bounded domain.

# Fields
- `coefs::Array{T, N}`: coefficients
- `lb::SVector{N, T}`: domain lower bound
- `ub::SVector{N, T}`: domain upper bound
"""
struct ChebyshevSeries{T, N}
    coefs::Array{T, N}
    lb::SVector{N, T}
    ub::SVector{N, T}
end


function ChebyshevSeries(coefs::Array{T, 0}, lb::SVector{0}, ub::SVector{0}) where T
    return coefs[]
end


"""A collection of `ChebyshevSeries` objects."""
struct ChebyshevCluster{T, N, M}
    series::NTuple{M, ChebyshevSeries{T, N}}
end


function ChebyshevCluster(series::ChebyshevSeries{T, N}...) where {T, N}
    M = length(series)
    return ChebyshevCluster{T, N, M}(series)
end


"""Converts a point `x` to its normalized coordinates in ``[-1, 1]^N``."""
function normalize(f::ChebyshevSeries{T, N}, x::SVector{N, T}) where {T, N}
    @. (2x - f.lb - f.ub) / (f.ub - f.lb)
end


function normalize(f::ChebyshevSeries{T, 1}, x::T) where T
    return normalize(f, SVector{1, T}(x))[]
end


"""Checks if the point `x` is in the domain of `f`."""
function contains(f::ChebyshevSeries{T, N}, x::SVector{N, T}) where {T, N}
    return all(f.lb .≤ x .≤ f.ub)
end


function contains(f::ChebyshevSeries{T, N}, x::AbstractVector{T}) where {T, N}
    return contains(f, SVector{N, T}(x))
end


function contains(f::ChebyshevSeries{T, 1}, x::T) where T
    return contains(f, SVector{1, T}(x))
end


function contains(g::ChebyshevCluster{T, N, M}, x::Union{AbstractVector{T}, T}) where {T, N, M}
    for i in 1:M
        if contains(g.series[i], x)
            return true, i
        end
    end
    return false, nothing
end

include("clenshaw.jl")
include("gradient.jl")
include("hessian.jl")

end