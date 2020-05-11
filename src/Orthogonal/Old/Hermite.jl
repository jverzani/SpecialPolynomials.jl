# https://arxiv.org/pdf/1901.01648.pdf
## We have the Chebyshev-Hermite, or Probabilist's Hermite polynomials

abstract type AbstractHermite{T} <: OrthogonalPolynomial{T} end

"""
    Hermite{T}

Implements the (physicists) [Hermite](https://en.wikipedia.org/wiki/Hermite_polynomials) polynomials. These are orthogonal with respect to weight function `w(x) = exp(-x^2)` over the real line.

The relation between the two is `H_n(x) = 2^(n/2)*H_{e_n}(sqrt(2)*x)`, with `H_{e_n}` being the probabilists Hermite polynomials.

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> p = Hermite([1,2,3])
Hermite(1⋅H_0(x) + 2⋅H_1(x) + 3⋅H_2(x))

julia> convert(Polynomial, p)
Polynomials.Polynomial(-5 + 4*x + 12*x^2)

julia> p4, p5 = basis.(Hermite, (4,5)) # verify orthogonality of H4, H5
(Hermite(1⋅H_4(x)), Hermite(1⋅H_5(x)))

julia> SpecialPolynomials.innerproduct(Hermite, p4, p5)
-8.709352683489158e-14
julia> n = 8
8

julia> Hn, Hn_1  = basis.(Hermite, (n, n-1)) # verify Hn' = 2n H_{n-1}
(Hermite(1⋅H_8(x)), Hermite(1⋅H_7(x)))

julia> derivative(Hn) - 2n*Hn_1
Hermite(0)
```
"""
struct Hermite{T <: Number} <: AbstractHermite{T}
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

Polynomials.domain(::Type{<:AbstractHermite}) = Polynomials.Interval(-Inf, Inf)

weight_function(::Type{<: Hermite})  = x -> exp(-x^2)
generating_function(::Type{<:Hermite}) = (t, x)  -> exp(2*t*x - t^2)
leading_coefficient(::Type{<:Hermite}) = 2^n

gauss_nodes_weights(P::Type{<:Hermite}, n)  = glaser_liu_rokhlin_gauss_nodes(basis(ScaledHermite,n))
has_fast_gauss_nodes_weights(::Type{<: Hermite}) = true

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
classical_σ(::Type{<:Hermite}) = x -> 1
classical_τ(::Type{<:Hermite}) = x -> -2x
function classical_hypergeometric(::Type{<:Hermite}, n, x)
    as = iseven(n) ? (-n ÷ 2, -(n-1)/2) : (-n/2, -(n-1)÷2)
    bs = ()
    (2x)^n * pFq(as, bs, -1/x^2)
end

(ch::Hermite{T})(x::S) where {T,S} = orthogonal_polyval(ch, x)

## Linearization
## leverage that of ChebyshevHermite
# Function for multiplication using linearization defined below
function Base.iterate(o::Linearization{P}, state =  nothing) where {P <: AbstractHermite}

    l, m, n = o.l, o.m, o.n
    l  > m + n && return nothing
    
    if state == nothing

        #k =  0
        j = 0
        p = min(l, n)
        q = l - p
        val = linearization_α(P, j, l, p, q) 

    else
        # we work with l + 2j as there is a parity consideration
        j, p, q, val  = state
        s = l + j # l +2j + p + q = l +2j + l = 2l + 2j, so s=l+j
        if p == 0 ||  q == m || q  > s
            j += 1
            l + 2j > n+m && return nothing #  l + 2j  is too  big now
            p = min(l + j, n)  # p <= s
            q = l + 2j - p
            q > s + 1 && return nothing
            val = linearization_α(P, j, l, p, q)
        else
            p -= 1
            q += 1

            λ = q/(p+1)*(s-(q-1))/(s-p)            
            val *= λ
        end

    end
                
    return (p,q,val), (j, p, q, val)
end


#  https://arxiv.org/pdf/1901.01648.pdf equation (74)
function linearization_α(P::Type{<:AbstractHermite}, j, l, n, m)
    s = l + j # pass in j to avoid s = divrem(l + m + n,  2) call
    val = gamma(1+m)*gamma(1+n)/gamma(1+s-m)/gamma(1+s-n)/gamma(1+s-l)
    val *= linearization_λ(P, l, m,n)
end
linearization_λ(::Type{<:Hermite}, l, m, n) = 2.0^((m+n-l)/2)


# A connection α(n,k) returns (k, α(k,k)), (k+1, α(k+1, k)), ..., (n, α(n,k))
#
# https://arxiv.org/pdf/1901.01648.pdf  Exercise 4
# shows formulas (13) and (14) with x^n = ∑ α(n, j) H_{n-2j}
# where α(n,j) depends on Hermite or ChebyshevHermite
# Here we push that difference into `connection_α(P,Q, i,j)`
#
# x^n = ∑_{j=0^n/2} n!/2^n 1/(n-2j)! 1/j! H_(n-2j) 
# let i = n - 2j j = (n-i)/2
# x^n = \sum_{n,n-2,... rem(n,2)} n!/2^n 1/i! 1/((n-i)/2)!  H_i
# 
# x^i = ∑ α(i,k) H_k with α(i,k) = i!/2^i * 1/(i! (i-k)/2)
# n = k  + 2j 
function Base.iterate(o::Connection{P, Q}, state=nothing) where
    {P <: AbstractHermite,
     Q <: Polynomials.StandardBasisPolynomial}

    k, n = o.k, o.n

    if state == nothing
        j  = 0
    else
        j = state + 1
    end

    i = k + 2j
    i > n && return nothing

    α = connection_α(P, Q, i, j)
    return(i, α), j
    
end

# compute
# n!/(2^n ⋅ (n-2j)! ⋅ j!)
function connection_α(::Type{<:Hermite},
                      ::Type{<:Polynomials.StandardBasisPolynomial},
                      n,j)
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

##
## ----
##

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

##
## --------------------------------------------------
##

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
julia> using Polynomials, SpecialPolynomials

julia> p = ChebyshevHermite([1,2,3])
ChebyshevHermite(1⋅He_0(x) + 2⋅He_1(x) + 3⋅He_2(x))

julia> convert(Polynomial, p)
Polynomials.Polynomial(-2 + 2*x + 3*x^2)
```

"""
struct ChebyshevHermite{T <: Number} <: AbstractHermite{T}
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

weight_function(::Type{ChebyshevHermite{T}}) where {T} = x -> exp(-x^2/2)
generating_function(::Type{<:ChebyshevHermite}) = (t, x)  -> exp(t*x - t^2/2)


# https://en.wikipedia.org/wiki/Hermite_polynomials
# we use Probabalists version
An(::Type{<:ChebyshevHermite}, n)  = 1
Bn(::Type{<:ChebyshevHermite}, n) = 0
Cn(::Type{<:ChebyshevHermite}, n) = -n

norm2(::Type{<:ChebyshevHermite}, n) = sqrt(2pi) * gamma(n+1)

(ch::ChebyshevHermite{T})(x::S) where {T,S} = orthogonal_polyval(ch, x)

linearization_λ(::Type{<:ChebyshevHermite}, l, m, n) = 1    


# compute α(n,j) for the ChebyshevHermite inversion
# See https://arxiv.org/pdf/1901.01648.pdf equations (13) and (14)
# connection defined in Hermite.jl
function connection_α(::Type{<:ChebyshevHermite},
                      ::Type{<:Polynomials.StandardBasisPolynomial},
                      n,j)
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


##
## --------------------------------------------------
##


##  ScaledHermite
## not exported, used for quadrature
## not even a polynomial, but satisfies 3-point recursion
"""
    ScaledHermitee

`H̃n(x) = exp(-x^2/2)/π^(1/4) * 1 / sqrt(2^n n!) Hn(x)`
"""
struct ScaledHermite{T <: Number} <: AbstractHermite{T}
    coeffs::Vector{T}
    var::Symbol
    function ScaledHermite{T}(coeffs::AbstractVector{T}, var::Symbol) where {T <: Number}
        length(coeffs) == 0 && return new{T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{T}(coeffs[1:last], var)
    end
end

Polynomials.@register ScaledHermite
basis_symbol(::Type{<:ScaledHermite}) = "H̃"


An(::Type{<:ScaledHermite}, n) = sqrt(2/(n+1))
Bn(::Type{<:ScaledHermite}, n) = 0
Cn(::Type{<:ScaledHermite}, n) = -sqrt(n/(n+1))
P0(::Type{<:ScaledHermite}, x) = one(x) * exp(-x^2/2)/pi^(1/4)
P1(P::Type{<:ScaledHermite}, x) = (An(P,0)*x .+ Bn(P,0)) * P0(P,x)
dP0(P::Type{<:ScaledHermite}, x) = -x*P0(P,x)
(ch::ScaledHermite{T})(x::S) where {T,S} = orthogonal_polyval(ch, x)

pqr(p::ScaledHermite) = (x,n) -> (p=1, q=0, r=(2n+1-x^2), dp=0, dq=0, dr=-2x)
pqr_start(p::ScaledHermite) = 0
pqr_symmetry(p::ScaledHermite) = true
pqr_weight(p::ScaledHermite, n, x, dπx) = 2*exp(-x^2)/(dπx*dπx)

linearization_λ(::Type{<:ScaledHermite}, l, m, n) = throw(MethodError())
connection_α(::Type{<:ScaledHermite}, n,j) = throw(MethodError())
