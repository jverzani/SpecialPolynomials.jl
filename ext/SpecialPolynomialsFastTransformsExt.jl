module SpecialPolynomialsFastTransformsExt

using SpecialPolynomials, FastTransforms

# code for FastTransforms.jl until that package is compatible with Polynomials v2.0.0

## Use FastTransform for conversion, as  possible
## FastTransforms.kind2string.(0:15)
# Legendre<--Chebyshev
function _convert(
    ::Type{Q},
    p::P,
) where {
    T<:AbstractFloat,
    Q<:Union{Legendre,OrthonormalLegendre},
    P<:Union{Chebyshev{T},OrthonormalChebyshev{T}},
}
    ps = coeffs(p)
    Q(cheb2leg(ps, normcheb=isorthonormal(P), normleg=isorthonormal(Q)))
end

# Chebyshev<--Legendre
function _convert(
    ::Type{Q},
    p::P,
) where {
    T<:AbstractFloat,
    Q<:Union{Chebyshev,OrthonormalChebyshev},
    P<:Union{Legendre{T},OrthonormalLegendre{T}},
}
    ps = coeffs(p)
    Q(leg2cheb(ps, normleg=isorthonormal(P), normcheb=isorthonormal(Q)))
end

# ultraspherical<--ultraspherical
function _convert(
    ::Type{Q},
    p::P,
) where {
    β,
    α,
    T<:AbstractFloat,
    Q<:Union{Gegenbauer{β},OrthonormalGegenbauer{β}},
    P<:Union{Gegenbauer{α,T},OrthonormalGegenbauer{α,T}},
}
    ps = coeffs(p)
    Q(ultra2ultra(ps, α, β, norm1=isorthonormal(P), norm2=isorthonormal(Q)))
end

# Jacobi<--Jacobi
function _convert(
    ::Type{Q},
    p::P,
) where {
    γ,
    δ,
    α,
    β,
    T<:AbstractFloat,
    Q<:Union{Jacobi{γ,δ},OrthonormalJacobi{γ,δ}},
    P<:Union{Jacobi{α,β,T},OrthonormalJacobi{α,β,T}},
}
    ps = coeffs(p)
    Q(jac2jac(ps, α, β, γ, δ, norm1=isorthonormal(P), norm2=isorthonormal(Q)))
end

# Laguerre<--Laguerre
function _convert(
    ::Type{Q},
    p::P,
) where {
    α,
    β,
    T<:AbstractFloat,
    Q<:Union{Laguerre{β},OrthonormalLaguerre{β}},
    P<:Union{Laguerre{α,T},OrthonormalLaguerre{α,T}},
}
    ps = coeffs(p)
    Q(lag2lag(ps, α, β, norm1=isorthonormal(P), norm2=isorthonormal(Q)))
end

# Jacobi<--ultraspherical
function _convert(
    ::Type{Q},
    p::P,
) where {
    γ,
    δ,
    α,
    T<:AbstractFloat,
    Q<:Union{Jacobi{γ,δ},OrthonormalJacobi{γ,δ}},
    P<:Union{Gegenbauer{α,T},OrthonormalGegenbauer{α,T}},
}
    ps = coeffs(p)
    Q(ultra2jac(ps, α, γ, δ, normultra=isorthonormal(P), normjac=isorthonormal(Q)))
end

# ultraspherical<--Jacobi
function _convert(
    ::Type{Q},
    p::P,
) where {
    α,
    γ,
    δ,
    T<:AbstractFloat,
    Q<:Union{Gegenbauer{α,T},OrthonormalGegenbauer{α,T}},
    P<:Union{Jacobi{γ,δ},OrthonormalJacobi{γ,δ}},
}
    ps = coeffs(p)
    Q(jac2ultra(ps, γ, δ, α, normjac=isorthonormal(P), normultra=isorthonormal(Q)))
end

# Jacobi<--Chebyshev
function _convert(
    ::Type{Q},
    p::P,
) where {
    γ,
    δ,
    T<:AbstractFloat,
    Q<:Union{Jacobi{γ,δ},OrthonormalJacobi{γ,δ}},
    P<:Union{Chebyshev{T},OrthonormalChebyshev{T}},
}
    ps = coeffs(p)
    Q(cheb2jac(ps, γ, δ, normcheb=isorthonormal(P), normjac=isorthonormal(Q)))
end

# Chebyshev<--Jacobi
function _convert(
    ::Type{Q},
    p::P,
) where {
    γ,
    δ,
    T<:AbstractFloat,
    Q<:Union{Chebyshev,OrthonormalChebyshev},
    P<:Union{Jacobi{γ,δ,T},OrthonormalJacobi{γ,δ,T}},
}
    ps = coeffs(p)
    Q(jac2cheb(ps, γ, δ, normjac=isorthonormal(P), normcheb=isorthonormal(Q)))
end

# ultraspherical<--Chebyshev
function _convert(
    ::Type{Q},
    p::P,
) where {
    α,
    T<:AbstractFloat,
    Q<:Union{Gegenbauer{α},OrthonormalGegenbauer{α}},
    P<:Union{Chebyshev{T},OrthonormalChebyshev{T}},
}
    ps = coeffs(p)
    Q(cheb2ultra(ps, α, normcheb=isorthonormal(P), normultra=isorthonormal(Q)))
end

# Chebyshev<--ultraspherical
function _convert(
    ::Type{Q},
    p::P,
) where {
    α,
    T<:AbstractFloat,
    Q<:Union{Chebyshev,OrthonormalChebyshev},
    P<:Union{Gegenbauer{α,T},OrthonormalGegenbauer{α,T}},
}
    ps = coeffs(p)
    Q(ultra2cheb(ps, α, normultra=isorthonormal(P), normcheb=isorthonormal(Q)))
end

end
