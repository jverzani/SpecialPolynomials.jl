## Jacobi Polynomials
@register2 Jacobi
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
weight_function(::Type{<:Jacobi{α, β}}) where {α, β} = x -> (1-x)^α *  (1+x)^β
generating_function(::Type{<:Jacobi{α, β}}) where {α, β} = (t,x) -> begin
    R = sqrt(1 - 2x*t+t^2)
    2^(α + β) * 1/R  * (1 - t + R)^(-α) * (1 + t + R)^(-β)
end

abcde(::Type{<:Jacobi{α,β}})  where {α,β} = NamedTuple{(:a,:b,:c,:d,:e)}((-1,0,1,-(α+β+2),β-α))

#kn =  1/2^n ⋅ choose(2n+α+β, n) = (2n+α+β)_n / (2^b  n!)
function kn(::Type{<:Jacobi{α,β}}, n::Int, ::Type{S}=Float64) where {S,α,β}
    a = α+β+n+1
    nn = n
    tot = one(n)/1
    for i in 0:n-1
        tot *= a/(2*nn)
        a += 1
        nn -= 1
    end
    tot
    
    #(one(α)*one(β) * generalized_binomial(2n+α+β, n))/2^n
end

function k1k0(P::Type{<:Jacobi{α,β}}, n::Int, ::Type{S}=Float64) where {S,α,β}
    γ = 2n + α + β
    val = one(S)
    if n == 0
        val *= (α+β)/2 + 1
    else
        num = (γ + 1) * (γ + 2)
        den = 2 * (n+1) * (γ - n + 1)
        iszero(den) && return k1k0(P, Val(n), S)
        val *= num / den
    end
    val
end
function k1k_1(P::Type{<:Jacobi{α,β}}, n, ::Type{S}=Float64) where {S,α, β}
    @assert n > 0

    γ = 2n + α + β
    val = one(S)

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

# overrides
Bn(::Type{<:Jacobi{α,β}}, ::Val{0}, ::Type{S}) where {α, β, S} = iszero(α+β) ? (α-β)/2 : (α-β)/(2(α+β+2))
Cn(P::Type{<:Jacobi{α,β}}, ::Val{1}, ::Type{S}) where {α,β,S} = -((α - β)^2 - (α + β + 2)^2)*(α + β + 4)/(8*(α + β + 2)^2)
ĉn(P::Type{<:Jacobi}, ::Val{0}, ::Type{S}) where {S} = zero(S)
ĉn(P::Type{<:Jacobi{α,β}}, ::Val{1}, ::Type{S}) where {α,β,S} = one(S) * ((α - β)^2 - (α + β + 2)^2)/((α + β + 1)*(α + β + 2)^2*(α + β + 3))
