## Gegenbauer Polynomials
@registerN Gegenbauer AbstractCCOP1 α
export  Gegenbauer

"""
   Gegenbauer{α, T <: Number}

The Gegenbauer polynomials have weight function
`(1-x^2)^(α-1/2)` over the domain `[-1,1]`. The parameter `α` is
specified in the constructor. These are also called the ultra-spherical polynomials.
The Legendre polynomials are the specialization  `Gegenbauer{1/2}`.

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> p =  Gegenbauer{1/2}([1,2,3])
Gegenbauer{0.5}(1⋅Cᵅ₀(x) + 2⋅Cᵅ₁(x) + 3⋅Cᵅ₂(x))

julia> convert(Polynomial, p)
Polynomial(-0.5 + 2.0*x + 4.5*x^2)
```

"""
Gegenbauer


basis_symbol(::Type{<:Gegenbauer{α}}) where {α} = "Cᵅ"
Polynomials.domain(::Type{<:Gegenbauer{α}}) where {α} = Polynomials.Interval(-1,1)

abcde(::Type{<:Gegenbauer{α}})  where {α} = NamedTuple{(:a,:b,:c,:d,:e)}((-1,0,1,-(2α+1),0))


k1k0(P::Type{<:Gegenbauer{α}}, k::Int) where {α} =(one(eltype(P))*2*(α+k))/(k+1)

function norm2(::Type{<:Gegenbauer{α}}, n) where{α}
    pi * 2^(1-2α) * gamma(n + 2α) / (gamma(n+1) * (n+α) * gamma(α)^2)
end
weight_function(::Type{<:Gegenbauer{α}}) where {α} = x -> (1-x^2)^(α-1/2)
generating_function(::Type{<:Gegenbauer{α}}) where {α} = (t,x) -> begin
    1/(1-2x*t +t^2)^α
end
function classical_hypergeometric(::Type{<:Gegenbauer{α}}, n, x) where {α}
    (α > -1/2 && !iszero(α)) || throw(ArgumentError("α > -1/2 and α ≠ 0 is necessary"))

    as = (-n, 2α+n)
    bs = (α + 1/2,)
    return Pochhammer_factorial(2α,n) * pFq(as, bs, (1-x)/2)
          
end



# overrides; big speedup if we avoid generating  from  a,b,c,d,e
B̃n(P::Type{<:Gegenbauer{α}}, n::Int) where {α} = zero(eltype(P))
C̃n(P::Type{<:Gegenbauer{α}}, n::Int) where {α} = (one(eltype(P))*n*(n + 2*α - 1))/(4*(n^2 + 2*n*α - n + α^2 - α))
b̂̃n(P::Type{<:Gegenbauer{α}}, n::Int) where {α} = zero(eltype(P))# one(S) *  4/α/(α+2)
b̂̃n(P::Type{<:Gegenbauer{α}}, ::Val{N}) where {α,N} = zero(eltype(P))# one(S) *  4/α/(α+2) 
ĉ̃n(P::Type{<:Gegenbauer{α}}, ::Val{0}) where {α} = zero(eltype(P)) #one(S) *  4/α/(α+2) 

