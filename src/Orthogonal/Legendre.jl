## Legendre Polynomials = Gegenbauer{1/2}

struct LegendreBasis <: AbstractCCOPBasis end

"""
    Legendre{T}

Implements the [Legendre](https://en.wikipedia.org/wiki/Legendre_polynomials) polynomials. These have weight function `w(x) = 1` over the domain `[-1,1]`.

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> p = Legendre([1,2,3])
Legendre(1⋅P₀(x) + 2⋅P₁(x) + 3⋅P₂(x))

julia> convert(Polynomial, p)
Polynomial(-0.5 + 2.0*x + 4.5*x^2)

julia> p2m, p2m1 = basis.(Legendre, (8,9)) # evaluation P_{2m+k}(-1) =  (-1)^k
(Legendre(1.0⋅P₈(x)), Legendre(1.0⋅P₉(x)))

julia> p2m(-1) == 1
false

julia> p2m1(-1) == -1
false

julia> n = 5  # verify  Rodrigues' formula
5

julia> x = Polynomial(:x)
Polynomial(1.0*x)

julia> derivative((x^2-1)^n, n) - 2^n *  factorial(n) * basis(Legendre, n)
LaurentPolynomial(0.0)

julia> p4, p5  =  basis.(Legendre, (4,5)) # verify  orthogonality  of  P₄,P₅
(Legendre(1.0⋅P₄(x)), Legendre(1.0⋅P₅(x)))

julia> SpecialPolynomials.innerproduct(Legendre, p4,  p5)
-1.111881270890998e-16
```
"""
Legendre = MutableDensePolynomial{LegendreBasis}
export Legendre
Polynomials._typealias(::Type{P}) where {P<:Legendre} = "Legendre"
Polynomials.basis_symbol(::Type{<:AbstractUnivariatePolynomial{LegendreBasis}})  = "P"

Polynomials.domain(::Type{<:LegendreBasis}) = Polynomials.Interval(-1, 1)

abcde(::Type{<:LegendreBasis}) = (a=-1, b=0, c=1, d=-2, e=0)
k0(::Type{<:LegendreBasis}) = 1
k1k0(::Type{<:LegendreBasis}, n::Int) =  (2n + 1) // (n + 1) #k1k0(Gegenbauer{1/2, eltype(P)}, n)

norm2(::Type{<:LegendreBasis}, n) = 2 / (2n + 1)
ω₁₀(::Type{<:LegendreBasis}, n) = sqrt(2n + 1) / (2n + 3)
weight_function(::Type{<:LegendreBasis}) = x -> one(x)
generating_function(::Type{<:LegendreBasis}) = (t, x) -> 1 / sqrt(1 - 2x * t + t^2)
function classical_hypergeometric(::Type{<:LegendreBasis}, n, x)
    as = (-n, n + 1)
    bs = (1,)
    pFq(as, bs, (1 - x) / 2)
end

# overrides

B̃n(::Type{<:LegendreBasis}, n::Int) = 0 #1
C̃n(::Type{<:LegendreBasis}, n::Int) = n^2 // (4n^2 - 1)

B̃n(::Type{<:LegendreBasis}, ::Val{0}) = B̃n(GegenbauerBasis{1 / 2}, 0)
C̃n(::Type{<:LegendreBasis}, ::Val{0}) = 0
b̂̃n(::Type{<:LegendreBasis}, ::Val{0}) = b̂̃n(GegenbauerBasis{1 / 2}, 0)
ĉ̃n(::Type{<:LegendreBasis}, ::Val{0}) = ĉ̃n(GegenbauerBasis{1 / 2}, 0)

function Polynomials.derivative(p::P) where {B<:LegendreBasis, T,X, P<:AbstractUnivariatePolynomial{B,T,X}}
    hasnan(p) && return ⟒(P){T,X}(T[NaN])

    d = degree(p)
    qs = zeros(T, d)
    for i in 0:(d - 1)
        gamma = 2 * i + 1
        qs[i + 1] = gamma * sum(p[j] for j in (i + 1):2:d)
    end

    dp = ⟒(P){T,X}(qs)
    return dp
end

##
## --------------------------------------------------
##

# Monic
struct MonicLegendreBasis <: AbstractCCOPBasis end
ϟ(::Type{MonicLegendreBasis}) = LegendreBasis
@register_monic(MonicLegendreBasis)

MonicLegendre = MutableDensePolynomial{MonicLegendreBasis}
Polynomials._typealias(::Type{P}) where {P<:MonicLegendre} = "MonicLegendre"
export MonicLegendre

## fast  gauss nodes
pqr_symmetry(::Type{<:MonicLegendreBasis}) = true
pqr_weight(P::Type{<:MonicLegendreBasis}, n, x, dπx) = (2) / (1 - x^2) / dπx^2


# Orthonormal
struct OrthonormalLegendreBasis <: AbstractCCOPBasis end
ϟ(::Type{OrthonormalLegendreBasis}) = LegendreBasis
@register_orthonormal(OrthonormalLegendreBasis)

OrthonormalLegendre = MutableDensePolynomial{OrthonormalLegendreBasis}
Polynomials._typealias(::Type{P}) where {P<:OrthonormalLegendre} = "OrthonormalLegendre"
export OrthonormalLegendre

# Shifted
struct ShiftedLegendreBasis <: AbstractCCOPBasis end
ϟ(::Type{ShiftedLegendreBasis}) = LegendreBasis
@register_shifted(ShiftedLegendreBasis, 2, -1)

ShiftedLegendre = MutableDensePolynomial{ShiftedLegendreBasis}
Polynomials._typealias(::Type{P}) where {P<:ShiftedLegendre} = "ShiftedLegendre"
export ShiftedLegendre

# MonicShifted
struct MonicShiftedLegendreBasis <: AbstractCCOPBasis end
ϟ(::Type{MonicShiftedLegendreBasis}) = ShiftedLegendreBasis
@register_monic(MonicShiftedLegendreBasis)

MonicShiftedLegendre = MutableDensePolynomial{MonicShiftedLegendreBasis}
Polynomials._typealias(::Type{P}) where {P<:MonicShiftedLegendre} = "MonicShiftedLegendre"
export MonicShiftedLegendre
