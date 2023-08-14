#@registerN Meixner AbstractCDOP2 γ μ

struct MeixnerBasis{γ, μ} <: AbstractCDOPBasis end
export Meixner

"""
    Meixner{γ,μ}

References: [Koekoek and Swarttouw §1.9](https://arxiv.org/pdf/math/9602214.pdf)
"""
Meixner{γ, μ} = MutableDensePolynomial{MeixnerBasis{γ, μ}} where {γ,μ}
Polynomials._typealias(::Type{P}) where {P<:Meixner{γ, μ}} where {γ, μ} = "Meixner"
Polynomials.basis_symbol(::Type{<:Meixner{γ, μ}}) where {γ, μ} = "m⁽ᵞᵝ⁾"

Polynomials.domain(::Type{<:MeixnerBasis{γ,μ}}) where {γ,μ} = Polynomials.Interval(0, Inf)
abcde(::Type{<:MeixnerBasis{γ,μ}}) where {γ,μ} =
    NamedTuple{(:a, :b, :c, :d, :e)}((0, 1, 0, μ - 1, γ * μ))

function kn(::Type{<:MeixnerBasis{γ,μ}}, n::Int) where {γ,μ}
    (1 - 1 / μ)^n
end

function k1k0(::Type{<:MeixnerBasis{γ,μ}}, n::Int) where {γ,μ}
    (1 - 1 / μ)
end
function k1k_1(::Type{<:MeixnerBasis{γ,μ}}, n::Int) where {γ,μ}
    (1 - 1 / μ)^(n > 0 ? 2 : 1)
end

function classical_hypergeometric(::Type{<:MeixnerBasis{γ,μ}}, n::Int, x) where {γ,μ}
    Pochhammer(γ, n) * pFq((-n, -x), γ, 1 - 1 / μ)
end

# Overrides
