
## Laguerre

struct Laguerre{T <: Number} <: OrthogonalPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function Laguerre{T}(coeffs::AbstractVector{T}, var::Symbol) where {T <: Number}
        length(coeffs) == 0 && return new{T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{T}(coeffs[1:last], var)
    end
end

export Laguerre
Polynomials.@register Laguerre

basis_symbol(::Type{<:Laguerre}) = "L"

# Pn = (An*x + Bn) * P_{n-1} + Cn P_{n-2}
# (k+1) L_{k+1} = (2k+1 - x) L_k  - k L_{k-1}
# n Ln = (2n-1 -  x)L_{n-1} - (n-2)L_{n-2}
# An  = -1.n
# Bn = (2n-1)/n
# Cn  = (n-2)/n
An(::Type{Laguerre{T}}, n) where {T <: Integer} = iszero(n) ? -1//1 : -1//n
An(::Type{<:Laguerre}, n) = iszero(n) ? -1/1 : -1/n
Bn(::Type{Laguerre{T}}, n) where {T <: Integer} = iszero(n) ? 1//1 : (2n-1)//n
Bn(::Type{<:Laguerre}, n) = iszero(n) ? 1/1 : (2n-1)/n
Cn(::Type{Laguerre{T}}, n) where {T <: Integer} = -(n-1)//n
Cn(::Type{<:Laguerre}, n) = -(n-1)/n

Polynomials.domain(::Type{<:Laguerre}) = Polynomials.Interval(0, Inf)
weight_function(::Type{Laguerre{T}}) where {T} = x -> exp(-x)
generating_function(::Type{<:Laguerre}) = (t, x)  -> exp(-t*x/(1-t)) / (1-t)

Polynomials.variable(::Type{P}, var::Polynomials.SymbolLike=:x) where {P <: Laguerre} = P([1, -1], var)
norm2(::Type{Laguerre{T}}, n) where {T} =  one(T)
(ch::Laguerre{T})(x::S) where {T,S} = orthogonal_polyval(ch, x)


## https://en.wikipedia.org/wiki/Laguerre_polynomials#Derivatives_of_generalized_Laguerre_polynomials
function Base.convert(P::Type{<:Laguerre}, p::Polynomial)

    d = degree(p)
    R = eltype(one(eltype(p))/1)
    qs = zeros(R, d+1)
    for i in 0:d
        qs[i+1] =  (iseven(i) ? 1 : -1) * sum(p[j] *   gamma(j+1)  * binomial(j, i) for j in i:d)
    end
    Laguerre(qs, p.var)
end
