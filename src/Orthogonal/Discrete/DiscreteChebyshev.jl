#@registerN DiscreteChebyshev AbstractCDOP2 α β

struct DiscreteChebyshevBasis{α,β} <: AbstractCDOPBasis end
export DiscreteChebyshev

#=
[DiscreteWeightFunction is on pause right now]
See  the example  in [`DiscreteWeightFunction`](@ref) for implementing  the [DiscreteChebyshev](https://en.wikipedia.org/wiki/Discrete_Chebyshev_polynomials)  polynomials  from Wikipedia.
=#

"""
    DiscreteChebyshev

This uses p22 of [Koepf and Schmersau](https://arxiv.org/pdf/math/9703217.pdf)  to define a two-parameter  family of *non* orthogonal polynomials.

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> import SpecialPolynomials: Δₓ, ∇ₓ

julia> α,β = 1/2, 1
(0.5, 1)

julia> P  = DiscreteChebyshev{α,β}
Polynomials.MutableDensePolynomial{SpecialPolynomials.DiscreteChebyshevBasis{0.5, 1}}

julia> i = 5
5

julia> yᵢ = basis(P, i)
DiscreteChebyshev{0.5,1}(1.0⋅K⁽ᵅᵝ⁾₅(x))

julia> x = variable(P)
DiscreteChebyshev{0.5,1}(- 2.0⋅K⁽ᵅᵝ⁾₀(x) + 2.0⋅K⁽ᵅᵝ⁾₁(x))

julia> a,b,c,d,e = SpecialPolynomials.abcde(P)
(a = 0, b = 0, c = 1, d = 0.5, e = 1)

julia> λᵢ  = -(a*i*(i-1)  + d*i)
-2.5

julia> Δₓ(∇ₓ(yᵢ)) +  (α*x + β) * Δₓ(yᵢ) ≈ -λᵢ*yᵢ # p22: "are not orthogonal, but satisfy the difference equation..."
true
```

"""
const DiscreteChebyshev = MutableDensePolynomial{DiscreteChebyshevBasis{α, β}} where {α, β}
Polynomials._typealias(::Type{P}) where {α,β,P<:DiscreteChebyshev{α,β}} = "DiscreteChebyshev{$α,$β}"
Polynomials.basis_symbol(::Type{<:AbstractUnivariatePolynomial{DiscreteChebyshevBasis{α,β}}}) where {α,β} = "K⁽ᵅᵝ⁾"


Polynomials.domain(::Type{<:DiscreteChebyshevBasis{α,β}}) where {α,β} =
    Polynomials.Interval(-Inf, Inf)
abcde(::Type{<:DiscreteChebyshevBasis{α,β}}) where {α,β} =
    NamedTuple{(:a, :b, :c, :d, :e)}((0, 0, 1, α, β))

function k1k0(::Type{<:DiscreteChebyshevBasis{α,β}}, n::Int) where {α,β}
    α
end
