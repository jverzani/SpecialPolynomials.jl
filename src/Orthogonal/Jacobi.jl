## Jacobi Polynomials
@registerN Jacobi AbstractCCOP2 α β
export Jacobi

"""
    Jacobi{α,  β, T}

Implements the [Jacobi](https://en.wikipedia.org/wiki/Jacobi_polynomials) polynomials. These have weight function `w(x) = (1-x)^α ⋅ (1+x)^β` over the domain `[-1,1]`. Many orthogonal polynomial types are special cases. The parameters are  specified to the constructors:

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> p = Jacobi{-1/2, -1/2}([0,0,1])
typename(Jacobi){-0.5,-0.5}(1⋅Jᵅᵝ₂(x))

julia> convert(Polynomial, p)
Polynomials.Polynomial(-0.375 + 0.75*x^2)

julia> monic(p) = (q=convert(Polynomial,p); q/q[end])
monic (generic function with 1 method)

julia> monic(p) ≈  monic(basis(Chebyshev, 2))
true

```

"""
Jacobi

basis_symbol(::Type{<:Jacobi{α,β}}) where {α,β} = "Jᵅᵝ"
Polynomials.domain(::Type{<:Jacobi{α, β}}) where {α, β} =
    Polynomials.Interval{β >= 0 ? Open : Closed, α >= 0 ? Open : Closed}(-1, 1)
abcde(::Type{<:Jacobi{α,β}})  where {α,β} = NamedTuple{(:a,:b,:c,:d,:e)}((-1,0,1,-(α+β+2),β-α))

k0(P::Type{<:Jacobi}) = one(eltype(P))
function k1k0(P::Type{<:Jacobi{α,β}}, n::Int) where {α,β}
    n == -1  &&  return one(eltype(P))/1
    γ = 2n + α + β
    val = one(eltype(P))
    if n == 0
        val *= (α+β)/2 + 1
    else
        num = (γ + 1) * (γ + 2)
        den = 2 * (n+1) * (γ - n + 1)
        iszero(den) && return k1k0(P, Val(n))
        val *= num / den
    end
    val
end

#kn =  1/2^n ⋅ choose(2n+α+β, n) = (2n+α+β)_n / (2^b  n!)
function leading_term(P::Type{<:Jacobi{α,β}}, n::Int) where {α,β}
    a = α+β+n+1
    nn = n
    tot = one(eltype(P))/1
    for i in 0:n-1
        tot *= a/(2*nn)
        a += 1
        nn -= 1
    end
    tot

end


weight_function(::Type{<:Jacobi{α, β}}) where {α, β} = x -> (1-x)^α *  (1+x)^β
generating_function(::Type{<:Jacobi{α, β}}) where {α, β} = (t,x) -> begin
    R = sqrt(1 - 2x*t+t^2)
    2^(α + β) * 1/R  * (1 - t + R)^(-α) * (1 + t + R)^(-β)
end

"""
    jacobi_eval(α, β, n, z)

Evaluate the nth basis element of `Pᵅᵝ` at `z` using

`Pₙᵅᵝ(z) = 2⁻ⁿ∑ (α+n, k) × (β+n, n-k) ⋅ (z+1)ᵏ ⋅ (z-1)ⁿ⁻ᵏ`

This alternate evaluation is  useful when `α` or `β ≤ -1` (and consequently the polynomials are not orthogonal)
"""
function jacobi_eval(α, β, n, z)
    T = eltype(z/2)
    tot = zero(T)
    for k in 0:n
        tot += generalized_binomial′(α+n,k) * generalized_binomial′(β+n, n-k) * ((z+1)/2)^k * ((z-1)/2)^(n-k)
    end
    tot
end


# gauss_nodes_weights(p::Type{P}, n) where {α, β, P <: Jacobi{α, β}} =
#     FastGaussQuadrature.gaussjacobi(n, α, β)

function classical_hypergeometric(::Type{<:Jacobi{α, β}}, n, x) where {α,β}

    (α ≤ -1 || β ≤ -1) && throw(ArgumentError("α and β must be > -1"))

    as = (-n, n+α+β+1)
    bs = (α+1,)

    Pochhammer_factorial(α+1,n) * pFq(as, bs, (1 - x)/2)
end



function norm2(::Type{<:Jacobi{α, β}}, n) where{α, β}
    α > -1 && β > -1 || throw(ArgumentError("α, β > -1 is necessary"))
    2^(α+β+1)/(2n+α+β+1) * (Γ(n+α+1) *  Γ(n + β +1))/(Γ(n+α+β+1)*Γ(n+1))
end
function ω₁₀(::Type{<:Jacobi{α, β}}, n) where{α, β}
    val = Γ(n+1+α+1)/Γ(n+α+1)
    val *= Γ(n + 1+ β +1) / Γ(n + β +1)
    val *= Γ(n+α+β+1) / Γ(n+1+α+β+1)
    val *= Γ(n+1) / Γ(n+1+1)
    val *= (2n+α+β+1) / (2(n+1)+α+β+1)
    sqrt(val)
end

# overrides
B̃n(P::Type{<:Jacobi{α,β}}, ::Val{0}) where {α, β} = iszero(α+β) ? (α-β)*one(eltype(P))/2 : (α-β)*one(eltype(P))/(2(α+β+2))
C̃n(P::Type{<:Jacobi{α,β}}, ::Val{1}) where {α,β} = -one(eltype(P)) * ((α - β)^2 - (α + β + 2)^2)/((α + β + 2)^2*(α + β + 3))

b̂̃n(P::Type{<:Jacobi{α,β}}, ::Val{0}) where {α,β} = one(eltype(P)) * NaN
ĉ̃n(P::Type{<:Jacobi}, ::Val{0})  = zero(eltype(P))
ĉ̃n(P::Type{<:Jacobi{α,β}}, ::Val{1}) where {α,β} = one(eltype(P)) * ((α - β)^2 - (α + β + 2)^2)/((α + β + 1)*(α + β + 2)^2*(α + β + 3))




##
## --------------------------------------------------
##

@registerN MonicJacobi AbstractCCOP2 α β
export MonicJacobi
ϟ(::Type{<:MonicJacobi{α,β}}) where {α,β} = Jacobi{α,β}
ϟ(::Type{<:MonicJacobi{α,β,T}}) where {α,β,T} = Jacobi{α,β,T}

@register_monic(MonicJacobi)


@registerN OrthonormalJacobi AbstractCCOP2 α β
export OrthonormalJacobi
ϟ(::Type{<:OrthonormalJacobi{α,β}}) where {α,β} = Jacobi{α,β}
ϟ(::Type{<:OrthonormalJacobi{α,β,T}}) where {α,β,T} = Jacobi{α,β,T}

@register_orthonormal(OrthonormalJacobi)
