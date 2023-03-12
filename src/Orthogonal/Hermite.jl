## Hermite
@register0 Hermite AbstractCCOP0
export Hermite

"""
    Hermite

The [Hermite](https://en.wikipedia.org/wiki/Hermite_polynomials) polynomials have two versions the physicists (`Hermite`  or  `H`) and the probablalists (`ChebyshevHermite` or  `Hₑ`). They are  related through  `Hᵢ(x) =  2^(i/2) Hₑᵢ(√2 x)`.

The Hermite   polynomials have weight  function `w(x)=exp(-x^2/2)` and domain the real line.

## Examples
```jldoctest
julia> using Polynomials,  SpecialPolynomials

julia> x = variable(Polynomial{Rational{Int}})
Polynomials.Polynomial(x)

julia> [basis(Hermite, i)(x) for i in 0:5]
6-element Vector{Polynomial{Float64, :x}}:
 Polynomials.Polynomial(1.0)
 Polynomials.Polynomial(2.0*x)
 Polynomials.Polynomial(-2.0 + 4.0*x^2)
 Polynomials.Polynomial(-12.0*x + 8.0*x^3)
 Polynomials.Polynomial(12.0 - 48.0*x^2 + 16.0*x^4)
 Polynomials.Polynomial(120.0*x - 160.0*x^3 + 32.0*x^5)

julia> [basis(ChebyshevHermite, i)(x) for i in 0:5]
6-element Vector{Polynomial{Float64, :x}}:
 Polynomials.Polynomial(1.0)
 Polynomials.Polynomial(1.0*x)
 Polynomials.Polynomial(-1.0 + 1.0*x^2)
 Polynomials.Polynomial(-3.0*x + 1.0*x^3)
 Polynomials.Polynomial(3.0 - 6.0*x^2 + 1.0*x^4)
 Polynomials.Polynomial(15.0*x - 10.0*x^3 + 1.0*x^5)
```

!!! note
    The Hermite family needs help, as the computed values for `Bn`,and,`Cn` are  both 0.
"""
Hermite

basis_symbol(::Type{<:Hermite}) = "H"
Polynomials.domain(::Type{<:Hermite}) = Polynomials.Interval(-Inf, Inf)

abcde(::Type{<:Hermite}) = NamedTuple{(:a, :b, :c, :d, :e)}((1, 0, 0, -2, 0))

k0(P::Type{<:Hermite}) = one(eltype(P))
function k1k0(P::Type{<:Hermite}, k::Int)
    val = 2 * one(eltype(P))
    val
end

norm2(::Type{<:Hermite}, n) = sqrt(pi) * 2^n * gamma(n + 1)
weight_function(::Type{<:Hermite}) = x -> exp(-x^2)
generating_function(::Type{<:Hermite}) = (t, x) -> exp(2 * t * x - t^2)
function classical_hypergeometric(::Type{<:Hermite}, n, x)
    as = iseven(n) ? (-n ÷ 2, -(n - 1) / 2) : (-n / 2, -(n - 1) ÷ 2)
    bs = ()
    (2x)^n * pFq(as, bs, -inv(x)^2)
end

## Overrides
# Use override here, as we get  An,Bn,Cn =  1, 0,0, otherwise
An(P::Type{<:Hermite}, n::Int) = 2 * one(eltype(P))
Bn(P::Type{<:Hermite}, n::Int) = zero(eltype(P))
Cn(P::Type{<:Hermite}, n::Int) = 2 * n * one(eltype(P))

β̃n(P::Type{<:Hermite}, n::Int) = zero(eltype(P))
γ̃n(P::Type{<:Hermite}, n::Int) = zero(eltype(P))
b̂̃n(P::Type{<:Hermite}, ::Val{N}) where {N} = zero(eltype(P)) #where {M} = error("Don't call me")#zero(S)
ĉ̃n(P::Type{<:Hermite}, ::Val{N}) where {N} = zero(eltype(P)) #where {M} = error("Don't call me")#zero(S)

## https://arxiv.org/pdf/1901.01648.pdf. Connection formula (14)
##  x^n  = n! sum(H_{n-2j}/ (2^j(n-2j)!j!) j = 0:floor(n/2))

Base.convert(P::Type{<:Hermite}, q::Polynomial) = connection(P, q)
function Base.iterate(
    o::Connection{P,Q},
    state=nothing,
) where {P<:Hermite,Q<:Polynomials.StandardBasisPolynomial}
    k, n = o.k, o.n

    if state == nothing
        i = k
        j = 0
        i > n && return nothing
        val = __hermite_lambda(i, k)
    elseif state[1] + 2 > n # terminate
        return nothing
    else
        j, val1 = state
        #val1 *= (i+1)*(i+2)/4/(j+1/4)
        j += 1
        i = k + 2j
        val = __hermite_lambda(i, k)
    end

    return (i, val), (j, val)
end

function __hermite_lambda(n, k)
    val = gamma(1 + n) / 2^n
    val /= gamma(1 + k)
    val /= gamma(1 + (n - k) / 2)
    val
end

# 2x more performant
function Polynomials.integrate(p::P, C::Number=0) where {T,X,P<:Hermite{T,X}}
    # int H_n = 1/(n+1) H_{n+1}
    R = eltype(one(T) / 1)
    d = degree(p)
    qs = zeros(R, d + 2)

    for i in 0:d
        qs[i + 2] = p[i] / (2(i + 1))
    end

    q = ⟒(P){R,X}(qs)
    q = q - q(0) + R(C)

    return q
end

##
## --------------------------------------------------
##

@register0 MonicHermite AbstractCCOP0
export MonicHermite
ϟ(::Type{<:MonicHermite}) = Hermite
ϟ(::Type{<:MonicHermite{T}}) where {T} = Hermite{T}
@register_monic(MonicHermite)

pqr_start(::Type{<:MonicHermite}) = 0
pqr_symmetry(::Type{<:MonicHermite}) = true
pqr_weight(P::Type{<:MonicHermite}, n, x, dπx) = (one(P) * 2) * exp(-x^2) / (dπx * dπx)

##
## --------------------------------------------------
##
@register0 ChebyshevHermite AbstractCCOP0
export ChebyshevHermite
"""
    ChebyshevHermite

Type for the Probabalist's  [Hermite](https://en.wikipedia.org/wiki/Hermite_polynomials) polynomials.

"""
ChebyshevHermite

basis_symbol(::Type{<:ChebyshevHermite}) = "Hₑ"
Polynomials.domain(::Type{<:ChebyshevHermite}) = Polynomials.Interval(-Inf, Inf)

# https://arxiv.org/pdf/1901.01648.pdf eqn 17
abcde(::Type{<:ChebyshevHermite}) = NamedTuple{(:a, :b, :c, :d, :e)}((0, 0, 1, -1, 0))

kn(P::Type{<:ChebyshevHermite}, n::Int)    = one(eltype(P))
k1k0(P::Type{<:ChebyshevHermite}, k::Int)  = one(eltype(P))
k1k_1(P::Type{<:ChebyshevHermite}, k::Int) = one(eltype(P))
ismonic(::Type{<:ChebyshevHermite})        = true

weight_function(P::Type{ChebyshevHermite}) = x -> exp(-(one(eltype(P)) * x^2) / 2)
generating_function(::Type{<:ChebyshevHermite}) = (t, x) -> exp(t * x - t^2 / 2)
function classical_hypergeometric(::Type{<:ChebyshevHermite}, n, x)
    2^(-n / 2) * classical_hypergeometric(Hermite, n, x / sqrt(2))
end

function gauss_nodes_weights(P::Type{<:ChebyshevHermite}, n)
    xs, ws = glaser_liu_rokhlin_gauss_nodes(basis(ChebyshevHermite, n))
    xs, ws
end

pqr_start(::Type{<:ChebyshevHermite}) = 0
pqr_symmetry(::Type{<:ChebyshevHermite}) = true
pqr_weight(P::Type{<:ChebyshevHermite}, n, x, dπx) =
    sqrt(2) * gamma(1 + n) * sqrt(pi) / (dπx)^2

## Overrides
#An(P::Type{<:ChebyshevHermite}, n::Int) = one(eltype(P))
#Bn(P::Type{<:ChebyshevHermite}, n::Int) = zero(eltype(P))
#Cn(P::Type{<:ChebyshevHermite}, n::Int) = -one(eltype(P))*n

##
## --------------------------------------------------
##

##  Multiplication
## TODO: work  on  types
AbstractHermite{T,X} = Union{Hermite{T,X},ChebyshevHermite{T,X}}

# Default is slow. Instead
# directly use a linearization formula from
# https://arxiv.org/pdf/1901.01648.pdf
# and  for Hermite, Hn(x) =2^(n/2)Hₑ_n(sqrt(2) x)
# so Hm⋅Hn = ∑ C_{m,n,j} 2^j H_{m+n-2j}
function ⊗(
    ::Type{<:AbstractHermite},
    p::P,
    q::Q,
) where {T,X,S,Y,P<:AbstractHermite{T,X},Q<:AbstractHermite{S,Y}}
    R = eltype(one(promote_type(T, S)) / 1)
    N, M = degree(p), degree(q)
    N == -1 && return zero(⟒(P){R,Y})
    M == -1 && return zero(⟒(P){R,X})
    N == 0 && return q * p[0]
    M == 0 && return p * q[0]
    X != Y && throw(ArgumentError("Variable names must match"))

    as = zeros(R, 1 + N + M)
    for i in eachindex(p)
        pᵢ = p[i]
        iszero(pᵢ) && continue
        for j in eachindex(q)
            qⱼ = q[j]
            pᵢqⱼ = pᵢ * qⱼ
            iszero(qⱼ) && continue
            λᵤ = one(R)
            for k in 0:min(i, j)
                as[1 + i + j - 2k] += pᵢqⱼ * λᵤ
                λᵤ *= (i - k) * (j - k) / (k + 1) * hermite_α(P)
            end
        end
    end
    ⟒(P){R,X}(as)
end

hermite_α(P::Type{<:ChebyshevHermite}) = 1
hermite_α(P::Type{<:Hermite}) = 2

#⊗(p::Hermite, q::Hermite) = linearization_product(p, q)
#⊗(p::ChebyshevHermite, q::ChebyshevHermite) = linearization_product(p, q)

# function Base.iterate(o::Linearization{P,R}, state =  nothing) where {P <: AbstractHermite,R}

#     l, m, n = o.l, o.m, o.n
#     l  > m + n && return nothing

#     if state == nothing

#         #k =  0
#         j = 0
#         p = min(l, n)
#         q = l - p
#         val = linearization_α(P, R, j, l, p, q)

#     else
#         # we work with l + 2j as there is a parity consideration
#         j, p, q, val  = state
#         s = l + j # l +2j + p + q = l +2j + l = 2l + 2j, so s=l+j
#         if p == 0 ||  q == m || q  > s
#             j += 1
#             l + 2j > n+m && return nothing #  l + 2j  is too  big now
#             p = min(l + j, n)  # p <= s
#             q = l + 2j - p
#             q > s + 1 && return nothing
#             val = linearization_α(P, R, j, l, p, q)
#         else
#             p -= 1
#             q += 1

#             λ = q/(p+1)*(s-(q-1))/(s-p)
#             val *= λ
#         end

#     end

#     return (p,q,val), (j, p, q, val)
# end

# #  https://arxiv.org/pdf/1901.01648.pdf equation (74)
# function linearization_α(P::Type{<:AbstractHermite}, R, j, l, n, m)
#     s = l + j # pass in j to avoid s = divrem(l + m + n,  2) call
#     val = one(R)
#     val *= gamma(1+m)*gamma(1+n)/gamma(1+s-m)/gamma(1+s-n)/gamma(1+s-l)
#     val *= linearization_λ(P, l, m,n)
# end
# linearization_λ(::Type{<:Hermite}, l, m, n) = 2.0^((m+n-l)/2)
# linearization_λ(::Type{<:ChebyshevHermite}, l, m, n) = 1
