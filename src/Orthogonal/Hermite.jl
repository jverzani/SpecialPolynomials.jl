# https://arxiv.org/pdf/1901.01648.pdf
## We have the Chebyshev-Hermite, or Probabilist's Hermite polynomials

"""
    Hermite{T}

Implements the (physicists) [Hermite](https://en.wikipedia.org/wiki/Hermite_polynomials) polynomials. These are orthogonal with respect to weight function `w(x) = exp(-x^2)` over the real line.

The relation between the two is `H_n(x) = 2^(n/2)*H_{e_n}(sqrt(2)*x)`, with `H_{e_n}` being the probabilists Hermite polynomials.

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

weight_function(::Type{<: Hermite})  = x -> exp(-x^2)
generating_function(::Type{<:Hermite}) = (t, x)  -> exp(2*t*x - t^2)

pqr(p::Hermite) = (x,n) -> (p=1, q=0, r=(2n+1-x^2), dp=0, dq=0, dr=-2x)
pqr_scale(p::Hermite) = (x,n) -> (sqrt(1/2(n+1)),
                                   n == 0 ? exp(-x^2/2)/(pi^(1/4)*sqrt(2)) : sqrt(1/(n*(n+1)))/2,
                                   0,
                                  n == 0 ? -x*exp(-x^2/2)/(pi^(1/4)*sqrt(2)) : 0)
pqr_start(p::Hermite) = 0
pqr_symmetry(p::Hermite) = true
pqr_weight(p::Hermite, n, x, dπx) = sqrt(pi)*2^(n+1)*gamma(n+1)/dπx^2
gauss_nodes_weights(p::P, n) where {P <: Hermite} = glaser_liu_rokhlin_gauss_nodes(Polynomials.basis(P,n))

# https://en.wikipedia.org/wiki/Hermite_polynomials
# Probabalist's version
#An(::Type{<:Hermite}, n)  = 1
#Bn(::Type{<:Hermite}, n) = 0
#Cn(::Type{<:Hermite}, n) = -n
#norm2(::Type{<:Hermite}, n) = sqrt(2pi) * gamma(n+1)

# Physicist's version
An(::Type{<:Hermite}, n)  = 2
Bn(::Type{<:Hermite}, n) = 0
Cn(::Type{<:Hermite}, n) = -2n
norm2(::Type{<:Hermite}, n) = sqrt(pi) * 2^n * gamma(n+1)


(ch::Hermite{T})(x::S) where {T,S} = orthogonal_polyval(ch, x)



## https://arxiv.org/pdf/1901.01648.pdf
##  (2x)^n  = n! sum(H_{n-2j}/ ((n-2j)!j!) j = 0:floor(n/2))
function Base.convert(P::Type{<:Hermite}, p::Polynomial)
    d = degree(p)
    R = eltype(one(eltype(p))/1)
    qs = zeros(R, d+1)
    for i in 0:d
        qs[i+1] = sum(p[jj] * _hermite_lambda(jj, j-1) for (j, jj) in enumerate(i:2:d))
    end
    Hermite(qs, p.var)
end

# compute
# n!/(2^n ⋅ (n-2j)! ⋅ j!)
function _hermite_lambda(n,j)
    tot = 1/1
    nn = 1
    for jj in 1:j
        tot *= nn/2/jj
        nn += 1
    end
    for jj in 1:(n-2j)
        tot *= nn/2/jj
        nn += 1
    end
    for nn in n-j+1:n
        tot *= nn/2
    end
    return tot
end



function Polynomials.derivative(p::Hermite{T}, order::Integer = 1) where {T}
    order < 0 && throw(ArgumentError("Order of derivative must be non-negative"))
    order == 0 && return convert(Hermite{T}, p)
    hasnan(p) && return Hermite(T[NaN], p.var)
    order > length(p) && return zero(Hermite{T})

    d = degree(p)
    qs = zeros(T, d+1-order)
    for i in order:d # Hn' =  2n Hn-1
        qs[i+1-order] = prod(2(1.0 + i - j) for j in 1:order)  * p[i]
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
        q[i+1] = p[i]/(2(i+1))
    end

    q = q - q(0) + R(C)

    return q
end
