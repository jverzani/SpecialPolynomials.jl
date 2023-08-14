## Gegenbauer Polynomials
struct GegenbauerBasis{α} <: AbstractCCOPBasis end
export Gegenbauer

"""
   Gegenbauer{α, T <: Number}

The Gegenbauer polynomials have weight function
`(1-x^2)^(α-1/2)` over the domain `[-1,1]`. The parameter `α` is
specified in the constructor. These are also called the ultra-spherical polynomials.
The Legendre polynomials are the specialization  `Gegenbauer{1/2}`.

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> p =  Gegenbauer{1/2}([1,2,3])
Gegenbauer(1⋅Cᵅ₀(x) + 2⋅Cᵅ₁(x) + 3⋅Cᵅ₂(x))

julia> convert(Polynomial, p)
Polynomial(-0.5 + 2.0*x + 4.5*x^2)
```

"""
Gegenbauer = MutableDensePolynomial{GegenbauerBasis{α}} where {α}
Polynomials._typealias(::Type{P}) where {α,P<:Gegenbauer{α}} = "Gegenbauer{$(α)}"
Polynomials.basis_symbol(::Type{<:AbstractUnivariatePolynomial{GegenbauerBasis{α}}}) where {α} = "Cᵅ"

Polynomials.domain(::Type{<:GegenbauerBasis{α}}) where {α} = Polynomials.Interval(-1, 1)

abcde(::Type{<:GegenbauerBasis{α}}) where {α} =(a=-1, b=0, c=1, d=-(2α + 1), e=0)

k0(::Type{<:Gegenbauer}) = 1
k1k0(::Type{<:GegenbauerBasis{α}}, k::Int) where {α} = (2 * (α + k)) / (k + 1)

function norm2(::Type{<:GegenbauerBasis{α}}, n) where {α}
    pi * 2^(1 - 2α) * gamma(n + 2α) / (gamma(n + 1) * (n + α) * gamma(α)^2)
end
function ω₁₀(::Type{<:GegenbauerBasis{α}}, n) where {α}
    val = gamma(n + 1 + 2α) / gamma(n + 2α)
    val /= (n + 1)
    val *= (n + α) / (n + 1 + α)
    sqrt(val)
end

weight_function(::Type{<:GegenbauerBasis{α}}) where {α} = x -> (1 - x^2)^(α - 1 / 2)
generating_function(::Type{<:GegenbauerBasis{α}}) where {α} =
    (t, x) -> begin
        1 / (1 - 2x * t + t^2)^α
    end
function classical_hypergeometric(::Type{<:GegenbauerBasis{α}}, n, x) where {α}
    (α > -1 / 2 && !iszero(α)) || throw(ArgumentError("α > -1/2 and α ≠ 0 is necessary"))

    as = (-n, 2α + n)
    bs = (α + 1 / 2,)
    return Pochhammer_factorial(2α, n) * pFq(as, bs, (1 - x) / 2)
end

# overrides; big speedup if we avoid generating  from  a,b,c,d,e
B̃n(::Type{<:GegenbauerBasis{α}}, n::Int) where {α} = 0
C̃n(::Type{<:GegenbauerBasis{α}}, n::Int) where {α} =
    ( n * (n + 2 * α - 1)) / (4 * (n^2 + 2 * n * α - n + α^2 - α))
b̂̃n(::Type{<:GegenbauerBasis{α}}, n::Int) where {α} = 0 # one(S) *  4/α/(α+2)
b̂̃n(::Type{<:GegenbauerBasis{α}}, ::Val{N}) where {α,N} = 0 # one(S) *  4/α/(α+2)
ĉ̃n(::Type{<:GegenbauerBasis{α}}, ::Val{0}) where {α} = 0 #one(S) *  4/α/(α+2)

# Orthonormal version
struct OrthonormalGegenbauerBasis{α} <: AbstractCCOPBasis end
ϟ(::Type{OrthonormalGegenbauerBasis{α}}) where {α} = GegenbauerBasis{α}
@register_orthonormal(OrthonormalGegenbauerBasis)

OrthonormalGegenbauer = MutableDensePolynomial{OrthonormalGegenbauerBasis{α}} where {α}
Polynomials._typealias(::Type{P}) where {α,P<:OrthonormalGegenbauer{α}} = "OrthonormalGegenbauer{$α}"
export OrthonormalGegenbauer
