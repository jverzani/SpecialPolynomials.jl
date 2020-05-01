abstract type AbstractBessel{T} <: OrthogonalPolynomial{T} end

"""
    Bessel{α,T}

Implements the [Bessel](https://dlmf.nist.gov/18.34) polynomials, introduced by [Krall and Frink](https://www.ams.org/journals/tran/1949-065-01/S0002-9947-1949-0028473-1/S0002-9947-1949-0028473-1.pdf) (with `b=2`). The  case `a=2` corresponds to the 
[Bessel](https://en.wikipedia.org/wiki/Bessel_polynomials) polynomials of Wikipedia. The Bessel  polynomials are not orthogonal over  a domain of the real  line, rather over an arbitray curve in the complex plane enclosing the  origin.  The weight  function is `ρ(x)=(2πi)^(-1)∑Γ(α)/Γ(α+n-1)(-β/x)^n`,   where `β=2`.

```jldoctest
julia> using Polynomials, SpecialPolynomials
```



"""
struct Bessel{α, T <: Number} <: OrthogonalPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function Bessel{α, T}(coeffs::AbstractVector{T}, var::Symbol=:x) where {α, T <: Number}
        α <= -1 && throw(ArgumentError("α > -1 is necessary"))
        length(coeffs) == 0 && return new{α, T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{α, T}(coeffs[1:last], var)
    end
end

export Bessel

Polynomials.@register1 Bessel
constructorof(P::Type{<:Bessel{α}})  where {α} = Bessel{α}

basis_symbol(::Type{<:Bessel{α}}) where {α} = "B^($(round(α,digits=2)))"


Polynomials.domain(::Type{<:Bessel}) = Polynomials.Interval(0, Inf)
Polynomials.variable(::Type{P}, var::Polynomials.SymbolLike=:x) where {α, P <: Bessel{α}} = P((2/α)*[-1, 1], var)

weight_function(::Type{<:Bessel{α}}) where {α} = z -> (2π*im)^(-1)*z^(α-2) * exp(-2/z)
generating_function(::Type{<:Bessel}) = (t, x)  -> error("XXX")
function leading_coefficient(::Type{<:Bessel{α}}, n) where {α}
    m = n - 1
    prod((α + i)/2 for i in m:2m)
end



## from Krall and  Frink
## taking b=2
function An(::Type{<:Bessel{α}}, n) where {α}
    (iszero(n) && α ∈ (1,2)) && return (2n+1)*(α/2)
    
    den = (n + α - 1) * (2n + α - 2)
    return  (2n+α)*(2n+α-1)*(2n+α-2)/2  / den

end
function Bn(::Type{<:Bessel{α}}, n) where {α}
    (iszero(n) && α ∈ (1,2)) && return 1#zero(α)
    
    den = (n + α - 1) * (2n + α - 2)
    return (α - 2) * (2n + α - 1) / den
end
function Cn(::Type{<:Bessel{α}}, n) where {α}
    (iszero(n) && α ∈ (1,2)) && return 1
    
    den = (n + α - 1) * (2n + α - 2)
    return n * (2n + α) / den
end

norm2(::Type{<:Bessel{α}}, n) where  {α} = -1^(n + α -1) * Γ(1+n) * 2^(α-1) / (Γ(n  + α -2) *  (2n +  α - 1))
classical_σ(::Type{<:Bessel})  = x -> x^2
classical_τ(::Type{<:Bessel{α}}) where {α} = x -> 2 + (α + 2) * x
function classical_hypergeometric(::Type{<:Bessel{α}}, n, x) where {α}
    as = (-n, n + α - 1)
    bs = ()
    pFq(as, bs, -x/2) #/ kn
end



(ch::Bessel{T})(x::S) where {T,S} = orthogonal_polyval(ch, x)

function Bessel{α,T}(q::Bessel{β,S})  where {α, β, T, S}
    β == α && return Bessel{α, T}(T.(coeffs(q)), q.var)
    connection(Bessel{α,T}, q)
end
Bessel{α}(q::Bessel{β,S})  where {α, β, S} = Bessel{α,S}(q)

# Conversion from J{α, β} to J{γ,δ}
# p29 of thesis
function Base.iterate(o::Connection{P, Q}, state=nothing) where
    {β, P <: Bessel{β},
     α, Q <: Bessel{α}}
    
    n,k = o.n, o.k

    if state  == nothing
        j = k
    else
        j, val = state
        j += 1  # likely j += 2 if parity condition
    end
    j > n && return nothing
    val = connection_α(P,Q,j,k)        
    (j,val), (j, val)
end

function connection_α(::Type{<:Bessel{β}},
                      ::Type{<:Bessel{α}},
                      n,m) where {α, β}

    α == β && return (n==m ? one(α) : zero(α))
    
    if m-n+β-α+1 < 0
        @show m,n,β,α
        return 0.0
    end
    val = (-1)^m * (2m+β+1)
    val *= Pochhammer(-n,m) * Pochhammer(n+α+1, m)
    val *= Γ(m + β+1)*Γ(β-α+1)
    val /= Γ(1+m) * Γ(n+m+β+2) * Γ(m-n+β-α+1)
    
    val
    
end

# # Inversion
# # p28 of thesis
# This doesn't work for us
# function Base.iterate(o::Connection{P, Q}, state=nothing) where
#     {α, P <: Bessel{α},
#      Q <: Polynomials.StandardBasisPolynomial}
  
#     n,k = o.n, o.k
#     if state  == nothing
#         j = k
#         val = connection_α(P,Q,j,k)                
#     else
#         j, val = state
#         j += 1
#         #λ = ...
#         #val *= λ
#         val = connection_α(P,Q,j,k)                        
#     end
#     j > n && return nothing


#     (j,val), (j, val)
  
# end

# function connection_α(::Type{<:Bessel{α}},
#                       ::Type{<:Polynomials.StandardBasisPolynomial},
#                       n,m) where {α}
         
#     val = (-2)^n 
#     val *= (2m + α + 1)
#     val *= Pochhammer(-n, m) * Γ(α + m + 1)
#     val /= Γ(1  +  m) * Γ(n + m + α + 2)
  
#     val
  
# end

# Koepf and Schmersau  p7 Corollary 1
# function Polynomials.integrate(p::P, C::S) where
#     {α,T, P<:Bessel{α,T},
#      S <: Number}
    
#     R = promote_type(eltype(one(T) / 1), S)
#     if hasnan(p) || isnan(C)
#         return ⟒(P)([NaN])
#     end
#     d = degree(p)
#     if d == 0
#         return ⟒(P)([C, p[0]])
#     end
    
#     as = zeros(R, d + 2)
    
#     @inbounds for n in 0:d

#         pn = p[n]
#         as[1 + n + 1] += pn * 2(n+1+α) / ((n+1)*(2n+α+1)*(2n+α+2))
#         as[1 + n] += pn * 4 / ((2n+α)*(2n+α+2))
#         if  n > 0
#             as[1 + n - 1] += pn * 2n / ((n+α)*(2n+α)*(2n+α+1))
#         end
#     end
    
#     ∫p = ⟒(P)(as,  p.var)
#     ∫p[0] = R(C) - ∫p(0)

#     return  ∫p
    
# end
