# https://arxiv.org/pdf/1901.01648.pdf
## We have the Chebyshev-Hermite, or Probabilist's Hermite polynomials

"""
    Hermite{T}

Implements the probabilists form of the [Hermite](https://en.wikipedia.org/wiki/Hermite_polynomials) polynomials. These are orthogonal with respect to weight function `w(x) = exp(-x^2/2)` over the real line.

```jldoctest
julia> p = Hermite([1,2,3])
Hermite(1⋅H_0(x) + 2⋅H_1(x) + 3⋅H_2(x))

julia> convert(Polynomial, p)
Polynomial(-2 + 2*x + 3*x^2)
```

"""
struct Hermite{T <: Number} <: OrthogonalPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function Hermite{T}(coeffs::AbstractVector{T}, var::Symbol) where {T <: Number}
        length(coeffs) == 0 && return new{T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{T}(coeffs[1:last], var)
    end
end

export Hermite

Polynomials.@register Hermite

basis_symbol(::Type{<:Hermite}) = "H"

Polynomials.domain(::Type{<:Hermite}) = Polynomials.Interval(-Inf, Inf)

weight_function(::Type{Hermite{T}}) where {T} = x -> exp(-x^2/2)
generating_function(::Type{<:Hermite}) = (t, x)  -> exp(t*x - t^2/2)


# https://en.wikipedia.org/wiki/Hermite_polynomials
# we use Probabalists version
An(::Type{<:Hermite}, n)  = 1
Bn(::Type{<:Hermite}, n) = 0
Cn(::Type{<:Hermite}, n) = -n

norm2(::Type{<:Hermite}, n) = sqrt(2pi) * gamma(n+1)

(ch::Hermite{T})(x::S) where {T,S} = orthogonal_polyval(ch, x)



## https://arxiv.org/pdf/1901.01648.pdf
##  x^n  = n! sum(H_{n-2j}/ (2^j(n-2j)!j!) j = 0:floor(n/2))
function Base.convert(P::Type{<:Hermite}, p::Polynomial)
    d = degree(p)
    R = eltype(one(eltype(p))/1)
    qs = zeros(R, d+1)
    for i in 0:d
        qs[i+1] = sum(p[jj] * _hermite_lambda(jj, j-1) for (j, jj) in enumerate(i:2:d))
    end
    Hermite(qs, p.var)
end

function _hermite_lambda(n,j)
    tot = 1/1
    for i in 1:j
        tot  /=  2
        tot *= n
        n -= 1
    end
    for i in 1:j
        tot /= i
        tot *= n
        n -= 1
    end
    tot
end



function Polynomials.derivative(p::Hermite{T}, order::Integer = 1) where {T}
    order < 0 && throw(ArgumentError("Order of derivative must be non-negative"))
    order == 0 && return convert(Hermite{T}, p)
    hasnan(p) && return Hermite(T[NaN], p.var)
    order > length(p) && return zero(Hermite{T})

    d = degree(p)
    qs = zeros(T, d+1-order)
    for i in order:d # Hn' =  n Hn-1
        qs[i+1-order] = prod((1.0 + i - j) for j in 1:order)  * p[i]
    end

    q = Hermite(qs, p.var)


end



function Polynomials.integrate(p::Hermite{T}, C::Number=0) where {T}
    # int H_n = 1/(n+1) H_{n+1}
    R = eltype(one(T)/1)
    d = degree(p)
    qs = zeros(R, d+2)
    q = Hermite(qs, p.var)

    for i in 0:d
        q[i+1] = p[i]/(i+1)
    end

    q = q - q(0) + R(C)

    return q
end
