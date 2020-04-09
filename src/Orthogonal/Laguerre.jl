"""
     Laguerre{T}

Implements the [Laguerre](https://en.wikipedia.org/wiki/Laguerre_polynomials) polynomials. These have weight function `w(x) = exp(-x)` over the domain `[0,oo)`.

```jldoctest
julia> p = Laguerre([1,2,3])
Laguerre(1⋅L_0(x) + 2⋅L_1(x) + 3⋅L_2(x))

julia> convert(Polynomial, p)
Polynomial(6//1 - 8//1*x + 3//2*x^2)
```
"""
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

Polynomials.domain(::Type{<:Laguerre}) = Polynomials.Interval(0, Inf)
Polynomials.variable(::Type{P}, var::Polynomials.SymbolLike=:x) where {P <: Laguerre} = P([1, -1], var)

weight_function(::Type{<: Laguerre}) = x -> exp(-x)
generating_function(::Type{<:Laguerre}) = (t, x)  -> exp(-t*x/(1-t)) / (1-t)

# (k+1) L_{k+1} = (2k+1 - x) L_k  - k L_{k-1}
# An  = -1/(n+1)
# Bn = (2n+1)/(n+1)
# Cn  =k/k+1
An(::Type{Laguerre{T}}, n) where {T <: Integer} = -1//(n+1)
An(::Type{<:Laguerre}, n) = -1/(n+1)
Bn(::Type{Laguerre{T}}, n) where {T <: Integer} =  (2n+1)//(n+1)
Bn(::Type{<:Laguerre}, n) =  (2n+1)/(n+1)
Cn(::Type{Laguerre{T}}, n) where {T <: Integer} = -n//(n+1)
Cn(::Type{<:Laguerre}, n) = -n/(n+1)

norm2(::Type{Laguerre{T}}, n) where {T} =  one(T)
norm2(::Type{Laguerre},n) =  1

# for gauss nodes
pqr(p::Laguerre) = (x,n) -> (p=x^2, q=x, r=-(1/4*x^2-(n+1/2)*x), dp=2x, dq=1, dr=-x/2+(n+1/2))
pqr_scale(p::Laguerre) = (x,n) ->  (one(x),  iszero(n) ? exp(-x/2) :  one(x), zero(x), iszero(n)  ? -1/2*exp(-x/2) : zero(x))
pqr_start(p::Laguerre, n) = 2/(4n+2)
pqr_symmetry(p::Laguerre) = false
pqr_weight(p::Laguerre, n, x, dπx) = 1/(x*dπx^2)
gauss_nodes_weights(p::P, n) where {P <: Laguerre} = glaser_liu_rokhlin_gauss_nodes(Polynomials.basis(P,n))

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
