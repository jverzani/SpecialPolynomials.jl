abstract type AbstractLaguerre{T} <: OrthogonalPolynomial{T} end


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
struct Laguerre{α, T <: Number} <: AbstractLaguerre{T}
    coeffs::Vector{T}
    var::Symbol
    function Laguerre{α, T}(coeffs::AbstractVector{T}, var::Symbol=:x) where {α, T <: Number}
        α <= -1 && throw(ArgumentError("α > -1 is necessary"))
        length(coeffs) == 0 && return new{α, T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{α, T}(coeffs[1:last], var)
    end
end

export Laguerre

Polynomials.@register1 Laguerre
constructorof(P::Type{<:Laguerre{α}})  where {α} = Laguerre{α}

basis_symbol(::Type{<:Laguerre{α}}) where {α} = iszero(α)  ? "L" : "L^($(round(α,digits=2)))"


Polynomials.domain(::Type{<:AbstractLaguerre}) = Polynomials.Interval(0, Inf)
Polynomials.variable(::Type{P}, var::Polynomials.SymbolLike=:x) where {α, P <: Laguerre{α}} = P([1+α, -1], var)

weight_function(::Type{<:Laguerre{α}}) where {α} = x  -> x^α * exp(-x)
generating_function(::Type{<:Laguerre{α}}) where {α} = (t,x) -> begin
    exp(-t*x/(1-t))/(1-t)^(α-1)
end


An(::Type{<:Laguerre{α}}, n) where {α} = -one(α)/(n+1)
Bn(::Type{<:Laguerre{α}}, n) where {α} = (2n + 1 +  α)/(n+1)
Cn(::Type{<:Laguerre{α}}, n) where {α} = -(n + α)/(n+1)

# compute <Jn, Jn>
function norm2(::Type{<:Laguerre{α}}, n) where{α}
    iszero(α) && return one(α)
    gamma(n + α + 1) / gamma(n+1)
end
classical_σ(::Type{<:Laguerre})  = x -> x
classical_τ(::Type{<:Laguerre{α}}) where {α} = x -> α + 1 - x
function classical_hypergeometric(::Type{<:Laguerre{α}}, n, x) where {α}
    α > -1 || throw(ArgumentError("α > -1 is required"))
    as = -n
    bs = α+1
    Pochhammer_factorial(α+1,n)*pFq(as, bs, x)
end

# for gauss nodes
gauss_nodes_weights(P::Type{<:Laguerre{0}}, n) = glaser_liu_rokhlin_gauss_nodes(basis(ScaledLaguerre,n))
has_fast_gauss_nodes_weights(::Type{<:Laguerre{0}}) = true

(ch::Laguerre)(x::S) where {T,S} = orthogonal_polyval(ch, x)

# Connection α -> β
## We have https://en.wikipedia.org/wiki/Laguerre_polynomials#Generalized_Laguerre_polynomials
## L^α_n = ∑ a(n,k) L^β_k with
## α(n,k) = (α - β - 1 + n - k, n - k)
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

## Inversion formula
# https://arxiv.org/pdf/math/9908148.pdf
# x^n = n! ∑_{k=0}^n (-1)^k gen_binom(n+α, n-k) L_k^α
#
function Base.iterate(o::Connection{P, Q}, state=nothing) where
    {α, P <: Laguerre{α},
     Q <: Polynomials.StandardBasisPolynomial}

    n,k = o.n, o.k

    if state  == nothing
        j = k
        val = connection_α(P,Q,j,k)
    else
        j, val = state
        j += 1  # likely j += 2 if parity condition
        #        val = connection_α(P,Q,j,k)
        λ= j*(j+α)/(j-k)  #  use ratio here
        val *= λ
    end

    j > n && return nothing
    
    (j,val), (j, val)
end

# e.g., https://arxiv.org/pdf/math/9908148.pdf equation (7)
function connection_α(::Type{<:Laguerre{α}},
                      ::Type{<:Polynomials.StandardBasisPolynomial},
                      n,k) where {α}

    tot = gamma(1+n)
    tot *= (-1)^k
    tot *= generalized_binomial(n+α, n-k)
    return tot
end


# p30 attributed to Watson
function linearization_α(::Type{<:P}, k, n, m) where {α, P <: Laguerre{α}}

    if k < m || k < n
        @show :k
        return 0.0 #0.0
        #     return Inf
    end
#    if k < abs(n-m)
#        return 1
   #     return -Inf
#    end

    
    val = (-2)^(n+m-k)
    val *= Γ(1+k)
    val /= Γ(1 + n + m - k) * Γ(1 + k - n) * Γ(1 + k - m)

    as = ((k-m-n)/2, (k - m - n + 1)/2, k + α + 1)
    bs = (k - n + 1, k - m + 1)
    z = 1.0

    val *= pFq(as, bs, z)

    val
end


# ∫L_n= L_n  - L_{n+1}
function Polynomials.integrate(p::P, C::S) where
    {α,T, P<:Laguerre{α,T},
     S <: Number}
    
    R = promote_type(eltype(one(T) / 1), S)
    if hasnan(p) || isnan(C)
        return ⟒(P)([NaN])
    end
    n = degree(p)
    if n == 1
        return ⟒(P)([C, p[0]])
    end
    
    as = zeros(R, n + 2)
    
    @inbounds for d in 0:n
        as[1 + d] += p[d]
        as[1 + d + 1] += -p[d] 
    end
    ∫p = ⟒(P)(as,  p.var)
    ∫p[0] = R(C) - ∫p(0)

    return  ∫p
end

# From ∫L_n =Ln-L_{n+1} we get L_n' = -(L_0+L_1+⋯+L_{n-1})
function Polynomials.derivative(p::P, order::Int=1) where  {α,T, P<:Laguerre{α,T}}
    order < 0 && throw(ArgumentError("Order of derivative must be non-negative"))

    R  = eltype(one(T)/1)

    order == 0 && return convert(⟒(P){R}, p)
    hasnan(p) && return ⟒(P)(R[NaN], p.var)
    order > length(p) && return zero(⟒(P){R})

    d = degree(p)
    as = zeros(R, d)

    for i in  1:d
        # (pᵢL_i)' =  -pᵢ(L_0 + ... + L_{i-1})
        as[1:i] .+= -p[i]
    end
    
    dp = ⟒(P)(as, p.var)
    return order > 1 ?  derivative(dp, order - 1) : dp

end

##
## --------------------------------------------------
##


##  ScaledLaguerre
"""
    ScaledLaguerre

`L̃n(x) = exp(-x/2) ⋅ Ln(x)`

Not a polynomial, but can still use polynomial evaluation machinery
"""
struct ScaledLaguerre{T <: Number} <: AbstractLaguerre{T}
    coeffs::Vector{T}
    var::Symbol
    function ScaledLaguerre{T}(coeffs::AbstractVector{T}, var::Symbol) where {T <: Number}
        length(coeffs) == 0 && return new{T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{T}(coeffs[1:last], var)
    end
end

Polynomials.@register ScaledLaguerre
basis_symbol(::Type{<:ScaledLaguerre}) = "L̃"


An(::Type{<:ScaledLaguerre}, n) = -1/(n+1)
Bn(::Type{<:ScaledLaguerre}, n) =  (2n+1)/(n+1)
Cn(::Type{<:ScaledLaguerre}, n) = -n/(n+1)
P0(::Type{<:ScaledLaguerre}, x) = exp(-x/2)
P1(P::Type{<:ScaledLaguerre}, x) = (An(P,0)*x .+ Bn(P,0)) * P0(P,x)
dP0(P::Type{<:ScaledLaguerre}, x) = -(0.5)*P0(P,x)
(ch::ScaledLaguerre{T})(x::S) where {T,S} = orthogonal_polyval(ch, x)

pqr(p::ScaledLaguerre) = (x,n) -> (p=x^2, q=x, r=-(1/4*x^2-(n+1/2)*x), dp=2x, dq=1, dr=-x/2+(n+1/2))
pqr_start(p::ScaledLaguerre, n) = 2/(4n+2)
pqr_symmetry(p::ScaledLaguerre) = false
pqr_weight(p::ScaledLaguerre, n, x, dπx) = exp(-x)/(x*dπx^2)
gauss_nodes_weights(P::ScaledLaguerre,  n)  = glaser_liu_rokhlin_gauss_nodes(basis(P,n))

