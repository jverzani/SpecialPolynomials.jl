@registerN Laguerre AbstractCCOP1 α
export  Laguerre

"""
   Laguerre{α, T <: Number}

The  Laguerre polynomials have weight function `x^α * exp(-x)` over the domain `[0, oo)`. The parameter `α` is specified through the constructor.

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> p = Laguerre{1/2}([1,2,3])
Laguerre(1⋅L^(0.5)_0(x) + 2⋅L^(0.5)_1(x) + 3⋅L^(0.5)_2(x))

julia> convert(Polynomial, p)
Polynomials.Polynomial(9.625 - 9.5*x + 1.5*x^2)
```

The Laguerre polynomials are the case `α=0`.

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> p = Laguerre{0}([1,2,3])
Laguerre(1⋅L_0(x) + 2⋅L_1(x) + 3⋅L_2(x))

julia> convert(Polynomial, p)
Polynomials.Polynomial(6.0 - 8.0*x + 1.5*x^2)

julia> phi(u, i) = derivative(u) -  u # verify Rodrigues' formula for small n; n! L_n = (d/dx-1)^n x^n
phi (generic function with 1 method)

julia> x = Polynomial(:x)
Polynomials.Polynomial(x)

julia> n = 7
7

julia> factorial(n) * basis(Laguerre{0}, n) - foldl(phi, 1:n, init=x^n)
Polynomials.Polynomial(0.0)
```


"""
Laguerre

basis_symbol(::Type{<:Laguerre{α}}) where {α} = "Lᵅ"
basis_symbol(::Type{<:Laguerre{0}}) = "L"
Polynomials.domain(::Type{<:Laguerre}) = Polynomials.Interval(0, Inf)

abcde(::Type{<:Laguerre{α}})  where {α} = NamedTuple{(:a,:b,:c,:d,:e)}((0,1,0,-1,α+1))

function kn(P::Type{<:Laguerre{α}}, n::Int)  where {α}
    (-one(eltype(P)))^n/factorial(n) #Γ(1+n)
end
k1k0(P::Type{<:Laguerre{α}}, k::Int) where {α} = k <= 0  ? -one(eltype(P)) : -one(eltype(P))/(k+1)
k1k_1(P::Type{<:Laguerre{α}}, k::Int) where {α} =  k <= 0  ? zero(eltype(P)) : one(eltype(P))/(k+1)/k

function norm2(::Type{<:Laguerre{α}}, n) where{α}
    iszero(α) && return one(α)
    gamma(n + α + 1) / gamma(n+1)
end

weight_function(::Type{<:Laguerre{α}}) where {α} = x  -> x^α * exp(-x)
generating_function(::Type{<:Laguerre{α}}) where {α} = (t,x) -> begin
    exp(-t*x/(1-t))/(1-t)^(α-1)
end
function classical_hypergeometric(::Type{<:Laguerre{α}}, n, x) where {α}
    α > -1 || throw(ArgumentError("α > -1 is required"))
    as = -n
    bs = α+1
    Pochhammer_factorial(α+1,n)*pFq(as, bs, x)
end

function gauss_nodes_weights(P::Type{<:Laguerre{α}}, n) where {α}
    α <= -1 && throw(ArgumentError("α  too small"))
    xs,  ws  =  glaser_liu_rokhlin_gauss_nodes(basis(MonicLaguerre{α},n))
    λ = kn(P,n)^2
    xs, ws/λ
end

        
## Overrides

# default  connection between Laguerre is popping out  0s
function Base.iterate(o::Connection{P, Q}, state=nothing) where
    {β, P <: Laguerre{β},
     α, Q <: Laguerre{α}}

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

    return (i, generalized_binomial(N,K)), i
    
end


##
## --------------------------------------------------
##
@registerN MonicLaguerre AbstractCCOP1 α
export MonicLaguerre
ϟ(::Type{<:MonicLaguerre{α}}) where {α} = Laguerre{α}
ϟ(::Type{<:MonicLaguerre{α,T}}) where {α,T} = Laguerre{α,T}
@register_monic(MonicLaguerre)

#
pqr_start(P::Type{MonicLaguerre{α}}, n) where {α} = 2/(4n+2α+2)
pqr_symmetry(P::Type{<:MonicLaguerre{α}})  where {α} = false
function pqr_weight(P::Type{<:MonicLaguerre{α}}, n, x, dπx) where {α}
    val =  one(eltype(P)) 
    val /= (x*dπx^2)
    val = 1/(x*dπx^2)    
    val
end
