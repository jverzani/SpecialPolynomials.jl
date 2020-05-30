## Jacobi Polynomials
@registerN Jacobi AbstractCCOP2 α β
export Jacobi

"""
    Jacobi{α,  β, T}

Implements the [Jacobi](https://en.wikipedia.org/wiki/Jacobi_polynomials) polynomials. These have weight function `w(x) = (1-x)^α ⋅ (1+x)^β` over the domain `[-1,1]`. Many orthogonal polynomial types are special cases. The parameters are  specified to the constructors:

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> p = Jacobi{-1/2, -1/2}([0,0,1])
Jacobi{-0.5,-0.5}(1⋅Jᵅᵝ₂(x))

julia> convert(Polynomial, p)
Polynomial(-0.375 + 0.75*x^2)

julia> monic(p) = (q=convert(Polynomial,p); q/q[end])
monic (generic function with 1 method)

julia> monic(p) ≈  monic(basis(Chebyshev, 2))
true

```

"""
Jacobi

basis_symbol(::Type{<:Jacobi{α,β}}) where {α,β} = "Jᵅᵝ"
Polynomials.domain(::Type{<:Jacobi{α, β}}) where {α, β} = Polynomials.Interval(-1, 1, β >= 0, α >= 0)
abcde(::Type{<:Jacobi{α,β}})  where {α,β} = NamedTuple{(:a,:b,:c,:d,:e)}((-1,0,1,-(α+β+2),β-α))


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

# function k1k_1(P::Type{<:Jacobi{α,β}}, n::Int) where {α, β}
#     @assert n > 0

#     γ = 2n + α + β
#     val = one(eltype(P))

#     if n == 1
#         num = (α+β+3)*(α+β+4)
#         den = 8
#     elseif n == 2
#         num = (α + β + 4) * (α + β + 5) * (α + β + 6)
#         den = 24 * (α + β + 2)
#     else
#         num = (γ - 1) * (γ - 0) * (γ + 1) * (γ + 2)
#         den = 4 * n * (n+1) * (n + α + β) * (n + α + β + 1)
#     end
    
#     iszero(den) && return k1k_1(P, Val(n), S)

#     val *= num /den
#     return val
# end

weight_function(::Type{<:Jacobi{α, β}}) where {α, β} = x -> (1-x)^α *  (1+x)^β
generating_function(::Type{<:Jacobi{α, β}}) where {α, β} = (t,x) -> begin
    R = sqrt(1 - 2x*t+t^2)
    2^(α + β) * 1/R  * (1 - t + R)^(-α) * (1 + t + R)^(-β)
end

function  gauss_nodes_weights(P::Type{<:Jacobi{α,β}}, n)  where {α,β}
    # we don't  have a  good starting point  unless α=β
    if α == β
        xs, ws = glaser_liu_rokhlin_gauss_nodes(basis(MonicJacobi{α,β},n))
        λ = kn(P,n)^2
        xs, ws/λ
    else
        J = jacobi_matrix(P, n)
        eig = eigen(J, extrema(P)...)
        wts =  π̃βn(P,0) * (eig.vectors[1,:]).^2
        eig.values,  wts
    end

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

# # overrides
# function B̃n(P::Type{<:Jacobi{α,β}}, n::Int) where {α, β}
#     val=one(eltype(P))
#     if  iszero(n)
#         return val * (α-β)/(α+β+2)
#     end
#     val *= (α - β)*(α + β)
#     val /= (2*n + α + β)*(2*n + α + β + 2)
#     val
# end
# function C̃n(P::Type{<:Jacobi{α,β}}, n::Int) where {α,β}
#     val=one(eltype(P))
#     val *= -n*(n + α + β)*(-4*n*(n + α + β) + 4*α + 4*β + (α - β)^2 - (α + β + 2)^2 + 4)
#     val /= (2*n + α + β)^2*(2*n + α + β - 1)*(2*n + α + β + 1)
#     val
# end

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

pqr_start(P::Type{MonicJacobi{α,α}}, n) where {α} =  0.0 # zero(eltype(P))
pqr_symmetry(P::Type{<:MonicJacobi{α,α}}) where {α} = true
function pqr_weight(P::Type{<:MonicJacobi{α,α}}, n, x, dπx) where {α}
    # https://www.chebfun.org/publications/HaleTownsend2013a.pdf
    # Eqn  (1.4)
    # Cnαβ should be  computed  using asymptotic formula  for larger n (§3.2.3)
    # XXX TThere i ss ome constant  that  makes this not work....
    β = α
    Cnαβ = 2^(α+β+1)  
    Cnαβ *= gamma(n + α +  1) * gamma(n + β + 1) / gamma(n + α + β + 1)
    Cnαβ /= gamma(1 + n)
    val = Cnαβ / (1-x^2) / dπx^2
    val
end
                      
