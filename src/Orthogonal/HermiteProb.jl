# https://arxiv.org/pdf/1901.01648.pdf
## We have the Chebyshev-Hermite, or Probabilist's Hermite polynomials

"""
    ChebyshevHermite{T}

Implements the probabilists form of the probabilists
[Hermite](https://en.wikipedia.org/wiki/Hermite_polynomials)
polynomials, also known as the Chebyshev-Hermite polynomials.. These
are orthogonal with respect to weight function `w(x) = exp(-x^2/2)`
over the real line.

```jldoctest
julia> p = ChebyshevHermite([1,2,3])
ChebyshevHermite(1⋅e_0(x) + 2⋅He_1(x) + 3⋅He_2(x))

julia> convert(Polynomial, p)
Polynomial(-2 + 2*x + 3*x^2)
```

"""
struct ChebyshevHermite{T <: Number} <: OrthogonalPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function ChebyshevHermite{T}(coeffs::AbstractVector{T}, var::Symbol) where {T <: Number}
        length(coeffs) == 0 && return new{T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{T}(coeffs[1:last], var)
    end
end

export ChebyshevHermite

Polynomials.@register ChebyshevHermite

basis_symbol(::Type{<:ChebyshevHermite}) = "He"

Polynomials.domain(::Type{<:ChebyshevHermite}) = Polynomials.Interval(-Inf, Inf)

weight_function(::Type{ChebyshevHermite{T}}) where {T} = x -> exp(-x^2/2)
generating_function(::Type{<:ChebyshevHermite}) = (t, x)  -> exp(t*x - t^2/2)


# https://en.wikipedia.org/wiki/Hermite_polynomials
# we use Probabalists version
An(::Type{<:ChebyshevHermite}, n)  = 1
Bn(::Type{<:ChebyshevHermite}, n) = 0
Cn(::Type{<:ChebyshevHermite}, n) = -n

norm2(::Type{<:ChebyshevHermite}, n) = sqrt(2pi) * gamma(n+1)

(ch::ChebyshevHermite{T})(x::S) where {T,S} = orthogonal_polyval(ch, x)



## https://arxiv.org/pdf/1901.01648.pdf
##  x^n  = n! sum(H_{n-2j}/ (2^j(n-2j)!j!) j = 0:floor(n/2))
function Base.convert(P::Type{<:ChebyshevHermite}, p::Polynomial)
    d = degree(p)
    R = eltype(one(eltype(p))/1)
    qs = zeros(R, d+1)
    for i in 0:d
        qs[i+1] = sum(p[jj] * _hermite_lambda_e(jj, j-1) for (j, jj) in enumerate(i:2:d))
    end
    ChebyshevHermite(qs, p.var)
end

function _hermite_lambda_e(n,j)
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



function Polynomials.derivative(p::ChebyshevHermite{T}, order::Integer = 1) where {T}
    order < 0 && throw(ArgumentError("Order of derivative must be non-negative"))
    order == 0 && return convert(ChebyshevHermite{T}, p)
    hasnan(p) && return ChebyshevHermite(T[NaN], p.var)
    order > length(p) && return zero(ChebyshevHermite{T})

    d = degree(p)
    qs = zeros(T, d+1-order)
    for i in order:d # Hn' =  n Hn-1
        qs[i+1-order] = prod((1.0 + i - j) for j in 1:order)  * p[i]
    end

    q = ChebyshevHermite(qs, p.var)


end



function Polynomials.integrate(p::ChebyshevHermite{T}, C::Number=0) where {T}
    # int H_n = 1/(n+1) H_{n+1}
    R = eltype(one(T)/1)
    d = degree(p)
    qs = zeros(R, d+2)
    q = ChebyshevHermite(qs, p.var)

    for i in 0:d
        q[i+1] = p[i]/(i+1)
    end

    q = q - q(0) + R(C)

    return q
end
