struct LaguerreBasis{α} <: AbstractCCOPBasis end

"""
   Laguerre{α, T <: Number}

The  Laguerre polynomials have weight function `x^α * exp(-x)` over the domain `[0, oo)`. The parameter `α` is specified through the constructor.

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> p = Laguerre{1/2}([1,2,3])
Laguerre(1⋅Lᵅ₀(x) + 2⋅Lᵅ₁(x) + 3⋅Lᵅ₂(x))

julia> convert(Polynomial, p)
Polynomial(9.625 - 9.5*x + 1.5*x^2)
```

The Laguerre polynomials are the case `α=0`.

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> p = Laguerre{0}([1,2,3])
Laguerre(1⋅L₀(x) + 2⋅L₁(x) + 3⋅L₂(x))

julia> convert(Polynomial, p)
Polynomial(6.0 - 8.0*x + 1.5*x^2)

julia> phi(u, i) = derivative(u) -  u # verify Rodrigues' formula for small n; n! L_n = (d/dx-1)^n x^n
phi (generic function with 1 method)

julia> x = Polynomial(:x)
Polynomial(1.0*x)

julia> n = 7
7

julia> factorial(n) * basis(Laguerre{0}, n) - foldl(phi, 1:n, init=x^n)
Polynomial(-5.4569682106375694e-12 + 1.4551915228366852e-11*x - 7.275957614183426e-12*x^2)
```


"""
Laguerre = MutableDensePolynomial{LaguerreBasis{α}} where {α}
export Laguerre
Polynomials._typealias(::Type{P}) where {α, P<:Laguerre{α}} = "Laguerre{$α}"
Polynomials.basis_symbol(::Type{<:AbstractUnivariatePolynomial{LaguerreBasis{α}}}) where {α} = "Lᵅ"
Polynomials.basis_symbol(::Type{<:AbstractUnivariatePolynomial{LaguerreBasis{0}}}) = "L"


Polynomials.domain(::Type{<:LaguerreBasis}) = Polynomials.Interval(0, Inf)

abcde(::Type{<:LaguerreBasis{α}}) where {α} = (a=0, b=1, c=0, d=-1, e=α + 1)

k0(::Type{<:LaguerreBasis}) = 1
k1k0(::Type{<:LaguerreBasis{α}}, k::Int) where {α} = -1 / (k + 1) # k >=0

function norm2(::Type{<:LaguerreBasis{α}}, n) where {α}
    iszero(α) && return one(α)
    Γ(n + α + 1) / Γ(n + 1)
end

function ω₁₀(::Type{<:LaguerreBasis{α}}, n) where {α}
    iszero(α) && return 1.0
    val = Γ(n + 1 + α + 1) / Γ(n + α + 1)
    val /= n + 1 # Γ(n+1)/Γ(n+2)
    sqrt(val)
end

weight_function(::Type{<:LaguerreBasis{α}}) where {α} = x -> x^α * exp(-x)
generating_function(::Type{<:LaguerreBasis{α}}) where {α} =
    (t, x) -> begin
        exp(-t * x / (1 - t)) / (1 - t)^(α + 1)
    end
function classical_hypergeometric(::Type{<:LaguerreBasis{α}}, n, x) where {α}
    α > -1 || throw(ArgumentError("α > -1 is required"))
    as = -n
    bs = α + 1
    Pochhammer_factorial(α + 1, n) * pFq(as, bs, x)
end

## Overrides

# default  connection between Laguerre is popping out  0s
function Base.iterate(
    o::Connection{P,Q},
    state=nothing,
) where {β,P<:Laguerre{β},α,Q<:Laguerre{α}}
    k, n = o.k, o.n

    if state == nothing
        i = k
        i > n && return nothing
    elseif state + 1 > n
        return nothing
    else
        i = state + 1
    end

    # α(i,k)
    K = i - k
    N = α - β - 1 + K

    return (i, generalized_binomial(N, K)), i
end

##
## --------------------------------------------------
##


struct MonicLaguerreBasis{α} <: AbstractCCOPBasis end
ϟ(::Type{MonicLaguerreBasis{α}}) where {α} = LaguerreBasis{α}
@register_monic(MonicLaguerreBasis)

MonicLaguerre = MutableDensePolynomial{MonicLaguerreBasis{α}} where {α}
Polynomials._typealias(::Type{P}) where {α,P<:MonicLaguerre{α}} = "MonicLaguerre{$α}"
export MonicLaguerre

#
struct OrthonormalLaguerreBasis{α} <: AbstractCCOPBasis end
ϟ(::Type{OrthonormalLaguerreBasis{α}}) where {α} = LaguerreBasis{α}
@register_orthonormal(OrthonormalLaguerreBasis)

OrthonormalLaguerre = MutableDensePolynomial{OrthonormalLaguerreBasis{α}} where {α}
Polynomials._typealias(::Type{P}) where {α, P<:OrthonormalLaguerre{α}} = "OrthonormalLaguerre{α}"
export OrthonormalLaguerre
