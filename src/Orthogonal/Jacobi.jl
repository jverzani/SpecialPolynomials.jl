## Jacobi Polynomials
@register2 Jacobi AbstractCCOP2
export Jacobi

"""
    Jacobi{α,  β, T}

Implements the [Jacobi](https://en.wikipedia.org/wiki/Jacobi_polynomials) polynomials. These have weight function `w(x) = (1-x)^α ⋅ (1+x)^β` over the domain `[-1,1]`. Many orthogonal polynomial types are special cases. The parameters are  specified to the constructors:

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> p = Jacobi{-1/2, -1/2}([0,0,1])
Jacobi(1⋅J^(α, β)_2(x))

julia> convert(Polynomial, p)
Polynomials.Polynomial(-0.375 + 0.75*x^2)

julia> monic(p::Polynomial) = p/p[end];

julia> (monic ∘ convert).(Polynomial, (p, basis(Chebyshev, 2))) # scaled version of each other
(Polynomials.Polynomial(-0.5 + 1.0*x^2), Polynomials.Polynomial(-0.5 + 1.0*x^2))
```

"""
Jacobi

basis_symbol(::Type{<:Jacobi{α,β}}) where {α,β} = "Jᵅᵝ"
Polynomials.domain(::Type{<:Jacobi{α, β}}) where {α, β} = Polynomials.Interval(-1, 1, β >= 0, α >= 0)
abcde(::Type{<:Jacobi{α,β}})  where {α,β} = NamedTuple{(:a,:b,:c,:d,:e)}((-1,0,1,-(α+β+2),β-α))

#kn =  1/2^n ⋅ choose(2n+α+β, n) = (2n+α+β)_n / (2^b  n!)
function kn(P::Type{<:Jacobi{α,β}}, n::Int) where {α,β}
    a = α+β+n+1
    nn = n
    tot = one(eltype(P))/1
    for i in 0:n-1
        tot *= a/(2*nn)
        a += 1
        nn -= 1
    end
    tot
    
    #(one(α)*one(β) * generalized_binomial(2n+α+β, n))/2^n
end

function k1k0(P::Type{<:Jacobi{α,β}}, n::Int) where {α,β}
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
function k1k_1(P::Type{<:Jacobi{α,β}}, n) where {α, β}
    @assert n > 0

    γ = 2n + α + β
    val = one(eltype(P))

    if n == 1
        num = (α+β+3)*(α+β+4)
        den = 8
    elseif n == 2
        num = (α + β + 4) * (α + β + 5) * (α + β + 6)
        den = 24 * (α + β + 2)
    else
        num = (γ - 1) * (γ - 0) * (γ + 1) * (γ + 2)
        den = 4 * n * (n+1) * (n + α + β) * (n + α + β + 1)
    end
    
    iszero(den) && return k1k_1(P, Val(n), S)

    val *= num /den
    return val
end

weight_function(::Type{<:Jacobi{α, β}}) where {α, β} = x -> (1-x)^α *  (1+x)^β
generating_function(::Type{<:Jacobi{α, β}}) where {α, β} = (t,x) -> begin
    R = sqrt(1 - 2x*t+t^2)
    2^(α + β) * 1/R  * (1 - t + R)^(-α) * (1 + t + R)^(-β)
end
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

# overrides
Bn(P::Type{<:Jacobi{α,β}}, ::Val{0}) where {α, β} = iszero(α+β) ? (α-β)*one(eltype(P))/2 : (α-β)*one(eltype(P))/(2(α+β+2))
Cn(P::Type{<:Jacobi{α,β}}, ::Val{1}) where {α,β} = -one(eltype(P)) * ((α - β)^2 - (α + β + 2)^2)/((α + β + 2)^2*(α + β + 3))

b̂n(P::Type{<:Jacobi{α,β}}, ::Val{0}) where {α,β} = one(eltype(P)) * NaN
ĉn(P::Type{<:Jacobi}, ::Val{0})  = zero(eltype(P))
ĉn(P::Type{<:Jacobi{α,β}}, ::Val{1}) where {α,β} = one(eltype(P)) * ((α - β)^2 - (α + β + 2)^2)/((α + β + 1)*(α + β + 2)^2*(α + β + 3))

