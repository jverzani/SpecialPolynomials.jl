#@registerN Charlier AbstractCDOP1 μ

struct CharlierBasis{μ} <: AbstractCDOPBasis end
export Charlier

"""
    Charlier{μ}

References: [Koekoek and Swarttouw §1.12](https://arxiv.org/pdf/math/9602214.pdf)
"""
const Charlier = MutableDensePolynomial{CharlierBasis{μ}} where {μ}
Polynomials._typealias(::Type{P}) where {P<:Charlier} = "Charlier"
Polynomials.basis_symbol(::Type{<:AbstractUnivariatePolynomial{CharlierBasis{μ}}}) where {μ} = "Cᵐᵘ"

abcde(::Type{<:CharlierBasis{μ}}) where {μ} = NamedTuple{(:a, :b, :c, :d, :e)}((0, 1, 0, -1, μ))

basis_symbol(::Type{<:CharlierBasis{μ}}) where {μ} = "cᵘ"
Polynomials.domain(::Type{<:CharlierBasis{μ}}) where {μ} = Polynomials.Interval(0, Inf)
weight_function(::Type{<:CharlierBasis{μ}}) where {μ} = x -> μ^x / gamma(1 + x)

function k1k0(::Type{<:CharlierBasis{μ}}, n::Int) where {μ}
    -1 * inv(μ)
end

function classical_hypergeometric(::Type{<:CharlierBasis{μ}}, n::Int, x) where {μ}
    pFq((-n, -x), (), -1 / μ)
end
