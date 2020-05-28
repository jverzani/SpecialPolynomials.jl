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
Gegenbaeur

basis_symbol(::Type{<:Gegenbauer{α}}) where {α} = "Cᵅ"
Polynomials.domain(::Type{<:Gegenbauer{α}}) where {α} = Polynomials.Interval(-1,1)

abcde(::Type{<:Gegenbauer{α}})  where {α} = NamedTuple{(:a,:b,:c,:d,:e)}((-1,0,1,-(2α+1),0))


function kn(::Type{<:Gegenbauer{α}}, n::Int) where {α,S}
    SpecialPolynomials.Pochhammer_factorial(α, n) * 2^n
end

function k1k0(P::Type{<:Gegenbauer{α}}, k::Int) where {α}
    S = eltype(P)
    iszero(k) && return 2*one(S)*α
    val = one(S)
    val *= 2(α+k)
    val /= k + 1
    val
end
function k1k_1(P::Type{<:Gegenbauer{α}}, k::Int) where {α}
    #iszero(k) && return k1k0(P,0)
    @assert k > 0
    
    S = eltype(P)
    
    val = one(S)
    if k == 1
        val *= 2α*(α+1)
    else
        val *= 4(α+k)*(α+k-1)
        val /= (k+1)*k
    end
    
    return val
    
end

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
#An(P::Type{<:Gegenbauer{α}}, n::Int) where {α} = 2*one(eltype(P))*(n+α)/(n+1)
#Bn(P::Type{<:Gegenbauer{α}}, n::Int) where {α} = zero(eltype(P))
#Cn(P::Type{<:Gegenbauer{α}}, n::Int) where {α} = -(one(eltype(P)) *  (n - 1 + 2α))/(n+1)


B̃n(P::Type{<:Gegenbauer{α}}, ::Val{0}) where  {α}  =  zero(eltype(P))
C̃n(P::Type{<:Gegenbauer{α}}, ::Val{1}) where  {α}  =  one(eltype(P))/2/(α+2)
b̂̃n(P::Type{<:Gegenbauer{α}}, n::Int) where {α} = zero(eltype(P))# one(S) *  4/α/(α+2)
b̂̃n(P::Type{<:Gegenbauer{α}}, ::Val{N}) where {α,N} = zero(eltype(P))# one(S) *  4/α/(α+2) 
ĉ̃n(P::Type{<:Gegenbauer{α}}, ::Val{0}) where {α} = zero(eltype(P)) #one(S) *  4/α/(α+2) 

