@registerN DiscreteChebyshev AbstractCDOP2 α β
export DiscreteChebyshev

"""
    DiscreteChebyshev

This uses p22 of [](https://arxiv.org/pdf/math/9703217.pdf)  to define a two-parameter  family of *non* orthogonal polynomials.
See  the example  in [`DiscreteWeightFunction`](@ref) for implementing  the [DiscreteChebyshev](https://en.wikipedia.org/wiki/Discrete_Chebyshev_polynomials)  polynomials  from Wikipedia.

```jldoctest
using SpecialPolynomials
import SpecialPolynomials: Δₓ, ∇ₓ
α,β = 1/2, 1
P  = DiscreteChebyshev{α,β}
i = 5
yᵢ = basis(P, i)
x = variable(P)
a,b,c,d,e = SpecialPolynomials.abcde(P)
λᵢ  = -(a*i*(i-1)  + d*i)
Δₓ(∇ₓ(yᵢ)) +  (α*x + β) * Δₓ(yᵢ) ≈ -λᵢ*yᵢ # p22: "are not orthogonal, but satisfy the difference equation..."
```

"""
DiscreteChebyshev

basis_symbol(::Type{<:DiscreteChebyshev{α,β}}) where {α,β} = "K⁽ᵅᵝ⁾"
Polynomials.domain(::Type{<:DiscreteChebyshev{α, β}}) where {α, β} = Polynomials.Interval(-Inf, Inf)
abcde(::Type{<:DiscreteChebyshev{α,β}})  where {α,β} = NamedTuple{(:a,:b,:c,:d,:e)}((0,0,1,α,β))

function k1k0(P::Type{<:DiscreteChebyshev{α,β}}, n::Int) where {α,β}
    α
end


