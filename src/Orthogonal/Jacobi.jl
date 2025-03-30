## Jacobi Polynomials

struct JacobiBasis{α, β} <: AbstractCCOPBasis end


"""
    Jacobi{α,  β, T}

Implements the [Jacobi](https://en.wikipedia.org/wiki/Jacobi_polynomials) polynomials. These have weight function `w(x) = (1-x)^α ⋅ (1+x)^β` over the domain `[-1,1]`. Many orthogonal polynomial types are special cases. The parameters are  specified to the constructors:

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> p = Jacobi{-1/2, -1/2}([0,0,1])
Jacobi{-0.5,-0.5}(1⋅Jᵅᵝ₂(x))

julia> convert(Polynomial, p)
Polynomial(-0.375 + 0.75*x^2)

julia> monic(p) = (q=convert(Polynomial,p); q/q[end])
monic (generic function with 1 method)

julia> monic(p) ≈  monic(basis(Chebyshev, 2))
true

```

Special cases include `Jacobi{α-1/2,α-1/2} = Gegenbauer{α}`, `Jacobi{-1/2,-1/2} = Chebyshev`, `Jacobi{1/2,1/2} = ChebyshevU`, and `Jacobi{0,0} = Legendre`.

"""
Jacobi = MutableDensePolynomial{JacobiBasis{α,β}} where {α,β}
export Jacobi
Polynomials._typealias(::Type{P}) where {α, β, P<:Jacobi{α, β}} = "Jacobi{$α,$β}"
Polynomials.basis_symbol(::Type{<:AbstractUnivariatePolynomial{JacobiBasis{α,β}}}) where {α,β} =  "Jᵅᵝ"

Polynomials.domain(::Type{<:JacobiBasis{α,β}}) where {α,β} =
    Polynomials.Interval{β >= 0 ? Open : Closed,α >= 0 ? Open : Closed}(-1, 1)

abcde(::Type{<:JacobiBasis{α,β}}) where {α,β} = (a=-1, b=0, c=1, d=-(α + β + 2), e=β - α)

k0(::Type{<:JacobiBasis}) = 1
function k1k0(::Type{<:JacobiBasis{α,β}}, n::Int) where {α,β}
    n == -1 && return 1
    γ = 2n + α + β
    val = 1
    if n == 0
        val *= (α + β) / 2 + 1
    else
        num = (γ + 1) * (γ + 2)
        den = 2 * (n + 1) * (γ - n + 1)
        iszero(den) && return k1k0(B, Val(n))
        val *= num / den
    end
    val
end

#kn =  1/2^n ⋅ choose(2n+α+β, n) = (2n+α+β)_n / (2^b  n!)
function leading_term(B::Type{<:JacobiBasis{α,β}}, n::Int) where {α,β}
    a = α + β + n + 1
    nn = n
    tot = 1
    for i in 0:(n - 1)
        tot *= a / (2 * nn)
        a += 1
        nn -= 1
    end
    tot
end

weight_function(::Type{<:JacobiBasis{α,β}}) where {α,β} = x -> (1 - x)^α * (1 + x)^β
generating_function(::Type{<:JacobiBasis{α,β}}) where {α,β} =
    (t, x) -> begin
        R = sqrt(1 - 2x * t + t^2)
        2^(α + β) * 1 / R * (1 - t + R)^(-α) * (1 + t + R)^(-β)
    end

"""
    jacobi_eval(α, β, n, z)

Evaluate the nth basis element of `Pᵅᵝ` at `z` using

`Pₙᵅᵝ(z) = 2⁻ⁿ∑ (α+n, k) × (β+n, n-k) ⋅ (z+1)ᵏ ⋅ (z-1)ⁿ⁻ᵏ`

This alternate evaluation is  useful when `α` or `β ≤ -1` (and consequently the polynomials are not orthogonal)
"""
function jacobi_eval(α, β, n, z)
    T = eltype(z / 2)
    tot = zero(T)
    for k in 0:n
        tot +=
            generalized_binomial′(α + n, k) *
            generalized_binomial′(β + n, n - k) *
            ((z + 1) / 2)^k *
            ((z - 1) / 2)^(n - k)
    end
    tot
end

function classical_hypergeometric(::Type{<:JacobiBasis{α,β}}, n, x) where {α,β}
    (α ≤ -1 || β ≤ -1) && throw(ArgumentError("α and β must be > -1"))

    as = (-n, n + α + β + 1)
    bs = (α + 1,)

    Pochhammer_factorial(α + 1, n) * pFq(as, bs, (1 - x) / 2)
end

function norm2(::Type{<:JacobiBasis{α,β}}, n) where {α,β}
    α > -1 && β > -1 || throw(ArgumentError("α, β > -1 is necessary"))
    2^(α + β + 1) / (2n + α + β + 1) * (Γ(n + α + 1) * Γ(n + β + 1)) /
    (Γ(n + α + β + 1) * Γ(n + 1))
end
function ω₁₀(::Type{<:JacobiBasis{α,β}}, n) where {α,β}
    val = Γ(n + 1 + α + 1) / Γ(n + α + 1)
    val *= Γ(n + 1 + β + 1) / Γ(n + β + 1)
    val *= Γ(n + α + β + 1) / Γ(n + 1 + α + β + 1)
    val *= Γ(n + 1) / Γ(n + 1 + 1)
    val *= (2n + α + β + 1) / (2(n + 1) + α + β + 1)
    sqrt(val)
end

# overrides
B̃n(::Type{<:JacobiBasis{α,β}}, ::Val{0}) where {α,β} =
    iszero(α + β) ? (α - β)  / 2 : (α - β) / (2(α + β + 2))
C̃n(::Type{<:JacobiBasis{α,β}}, ::Val{1}) where {α,β} =
    -((α - β)^2 - (α + β + 2)^2) / ((α + β + 2)^2 * (α + β + 3))

b̂̃n(::Type{<:JacobiBasis{α,β}}, ::Val{0}) where {α,β} = NaN
ĉ̃n(::Type{<:Jacobi}, ::Val{0}) = 0
ĉ̃n(B::Type{<:JacobiBasis{α,β}}, ::Val{1}) where {α,β} =
    ((α - β)^2 - (α + β + 2)^2) /
    ((α + β + 1) * (α + β + 2)^2 * (α + β + 3))

##
## --------------------------------------------------
##
struct ShiftedJacobiBasis{α,β} <: AbstractCCOPBasis end
ϟ(::Type{ShiftedJacobiBasis{α,β}}) where {α,β}= JacobiBasis{α,β}
@register_shifted(ShiftedJacobiBasis, 2, -1)
"""
    ShiftedJacobi

Shifted Jacobi polynomial constructor. `P̃(x) = P(2x-1)`.
Shifted Jacobi are orthogonal on ``[0,1]``.
"""
ShiftedJacobi = MutableDensePolynomial{ShiftedJacobiBasis{α,β}} where {α,β}
Polynomials._typealias(::Type{P}) where {P<:ShiftedJacobi} = "ShiftedJacobi"
export ShiftedJacobi
Polynomials.domain(::Type{<:ShiftedJacobiBasis{α,β}}) where {α,β} =
    Polynomials.Interval{β >= 0 ? Open : Closed,α >= 0 ? Open : Closed}(0, 1)
weight_function(::Type{<:ShiftedJacobiBasis{α,β}}) where {α,β} = x -> (1 - x)^α * x^β

# helpers for n = 0
# use ζ₀(n)Rₙ + ζ₁(n)Rₙ₊₁ + ζ₂(n)Rₙ₊₂ = 0
function ζ₀(P::Type{<:ShiftedJacobiBasis{α,β}}, n::Int) where {α,β}
    σ = α + β + 1
    -2(n+α+1)*(n+β+1)*(2n+σ+3)
end
function ζ₁ₓ(P::Type{<:ShiftedJacobiBasis{α,β}}, n::Int) where {α,β}
    σ = α + β + 1
    (2n+σ+2) * (2n+σ+1) * (2n+σ+3) * 2
end
function ζ₁c(P::Type{<:ShiftedJacobiBasis{α,β}}, n::Int) where {α,β}
    σ = α + β + 1
    (2n+σ+2) * ((2n+σ+1) * (2n+σ+3) * (-1) + α^2 - β^2)
end
function ζ₂(P::Type{<:ShiftedJacobiBasis{α,β}}, n::Int) where {α,β}
    σ = α + β + 1
    -2 * (n+2) * (n+σ+1) * (2n+σ+1)
end

# Pₙ₊₁ = (Aₙ*x + Bₙ)Pₙ - CₙPₙ₋₁
function An(P::Type{<:ShiftedJacobiBasis{α,β}}, n::Int) where {α,β}
    n == 0 && return α + β + 2
    -ζ₁ₓ(P,n-1)/ζ₂(P,n-1)
end
function Bn(P::Type{<:ShiftedJacobiBasis{α,β}}, n::Int) where {α,β}
    if iszero(n)
        num,den = (α^2 - β^2 - (α + β)*(α + β + 2)), (2*(α + β))
        iszero(num) && return -1 - β
        return num/den
    end
    -ζ₁c(P,n-1)/ζ₂(P,n-1)
end
function Cn(P::Type{<:ShiftedJacobiBasis{α,β}}, n::Int) where {α,β}
    if n == 0
        num, den = α*β*(α + β + 2), ((α + β)*(α + β + 1))
        if iszero(num)
            if iszero(α)
                d = (β+2)^2
                return 2/d2 + 5β/(2d) + β^2/(2d)
            elseif iszero(β)
                d = (α + 2)^2
                return 2/d2 + 5α/(2d) + α^2/(2d)
            end
        end
        return num/den
    end
    ζ₀(P,n-1)/ζ₂(P,n-1)
end

# --- monic
struct MonicJacobiBasis{α, β} <: AbstractCCOPBasis end
ϟ(::Type{MonicJacobiBasis{α, β}}) where {α, β} = JacobiBasis{α, β}
@register_monic(MonicJacobiBasis)

MonicJacobi = MutableDensePolynomial{MonicJacobiBasis{α, β}} where {α, β}
Polynomials._typealias(::Type{P}) where {α,β,P<:MonicJacobi{α,β}} = "MonicJacobi{$α, $β}"
export MonicJacobi

#--- orthogonal
struct OrthonormalJacobiBasis{α,β} <: AbstractCCOPBasis end
ϟ(::Type{OrthonormalJacobiBasis{α,β}}) where {α, β} = JacobiBasis{α, β}
@register_orthonormal(OrthonormalJacobiBasis)

OrthonormalJacobi = MutableDensePolynomial{OrthonormalJacobiBasis{α, β}} where {α, β}
Polynomials._typealias(::Type{P}) where {α, β, P<:OrthonormalJacobi{α, β}} = "OrthonormalJacobi{$α, $β}"
export OrthonormalJacobi
