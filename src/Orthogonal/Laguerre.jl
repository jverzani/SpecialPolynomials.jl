@register1 Laguerre
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
Polynomials.domain(::Type{<:Laguerre}) = Polynomials.Interval(0, Inf)
weight_function(::Type{<:Laguerre{α}}) where {α} = x  -> x^α * exp(-x)
generating_function(::Type{<:Laguerre{α}}) where {α} = (t,x) -> begin
    exp(-t*x/(1-t))/(1-t)^(α-1)
end


abcde(::Type{<:Laguerre{α}})  where {α} = NamedTuple{(:a,:b,:c,:d,:e)}((0,1,0,-1,α+1))

function kn(::Type{<:Laguerre{α}}, n::Int, ::Type{S}=Float64)  where {α,S}
    (-one(S))^n/Γ(1+n)
end
k1k0(::Type{<:Laguerre{α}}, k, ::Type{S}=Float64) where {S,α} = -one(S)/(k+1)
k1k_1(::Type{<:Laguerre{α}}, k, ::Type{S}=Float64) where {S,α} =  one(S)/(k+1)/k

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

Base.convert(::Type{Q}, p::P) where  {α,Q<:Laguerre{α},β,T,N,P<:Laguerre{β,T,N}} = connection_n(Q,  p)

