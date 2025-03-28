#@registerN Hahn AbstractCDOP3 α β 𝐍

export Hahn
struct HahnBasis{α,β,𝐍} <: AbstractCDOPBasis end
"""
    Hahn{α,β,𝐍}

References: [Koekoek and Swarttouw §1.5](https://arxiv.org/pdf/math/9602214.pdf)

!!! note
    In  [Koekoek and Swarttouw](https://arxiv.org/pdf/math/9602214.pdf) sections 1.3, 1.4, and 1.6  are other  Hahn-type polynomials,  not  implemented here.

"""
Hahn{α,β,𝐍} = MutableDensePolynomial{HahnBasis{α,β,𝐍}} where {α,β,𝐍}

abcde(::Type{<:HahnBasis{α,β,𝐍}}) where {α,β,𝐍} =
    NamedTuple{(:a, :b, :c, :d, :e)}((1, -(β + 𝐍 + 1), 0, α + β + 2, -𝐍 * (α + 1)))

Polynomials.basis_symbol(::Type{<:Hahn{α,β,𝐍}}) where {α,β,𝐍} = "Q⁽ᵅᵝ⁾"
Polynomials.domain(::Type{<:HahnBasis{α,β,𝐍}}) where {α,β,𝐍} = Polynomials.Interval(0, 𝐍)
weight_function(::Type{<:HahnBasis{α,β,𝐍}}) where {α,β,𝐍} =
    x -> generalized_binomial(α + x, x) * generalized_binomial(𝐍 + β - x, 𝐍 - x)

# use   recurrence relation  1.5.3 of Koekoek and Swarttouw
# -xQᵢ = AQᵢ₊₁ -  (A+C)Qᵢ + CQᵢ₋₁
#  or Qᵢ₊₁ = (1/A ⋅ x - (1 + C/A)) Qᵢ + C/A ⋅ Qᵢ₋₁
# Using kᵢ = A⁻¹ᵢ₋₁ ⋅ A⁻¹ᵢ₋₂ ⋯ A⁻¹₀; so kᵢ₊₁/kᵢ = A⁻¹ᵢ, kᵢ₊₁/kᵢ₋ᵢ = A⁻¹ᵢA⁻¹ᵢ₋₁
#
# function Aᵢ⁻¹(n,α,β,𝐍)
#     num = (2n+α+β + 1) * (2n+α+β+2)
#     den = (n+α + β + 1)*(n + α + 1) *(𝐍-n)
#     num/den
# end

# define kn, k1k_1 through defaults
function k1k0(B::Type{<:HahnBasis{α,β,𝐍}}, n::Int) where {α,β,𝐍}
    num = (2 * n + α + β + 1) * (2 * n + α + β + 2)
    den = (-n + 𝐍) * (n + α + 1) * (n + α + β + 1)

    - num / den
end

# end

function classical_hypergeometric(::Type{<:HahnBasis{α,β,𝐍}}, n::Int, x) where {α,β,𝐍}
    pFq((-n, -x, n + 1 + α + β), (α + 1, -𝐍), 1)
end
