abstract type AbstractBessel{T} <: OrthogonalPolynomial{T} end

"""
    Bessel{T}

Implements the

 [Bessel](https://en.wikipedia.org/wiki/Bessel_polynomials) polynomials. These have weight function `w(x) = 1` over the domain `[-1,1]`.

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

basis_symbol(::Type{<:Bessel{α}}) where {α} = "B^($α)"


Polynomials.domain(::Type{<:Bessel}) = Polynomials.Interval(0, Inf)
Polynomials.variable(::Type{P}, var::Polynomials.SymbolLike=:x) where {P <: Bessel} = P([-1, 1], var)

weight_function(::Type{<:Bessel{α}}) where {α} = x -> x^α*exp(-2/x)
generating_function(::Type{<:Bessel}) = (t, x)  -> XXX
leading_coefficient(::Type{<:Bessel{α}}, n) where {α} = prod((1+α+i)/2 for i=1:n)




function An(::Type{<:Bessel{α}}, n) where {α} 
    iszero(α) && return 2n+1
    (2n+α)*(2n+α-1)/(2*(n+α-1))
end
function Bn(::Type{<:Bessel{α}}, n) where {α}
    iszero(α) && return iszero(n) ? 1 : 0
    (α-2)*(2n+α-1)/((n+α-1)*(2n+α+2))
end
function Cn(::Type{<:Bessel{α}}, n) where {α}
    iszero(α) && return 1
    n*(2n+α)/((n+α-1)*(2n+α-2))
end

norm2(::Type{<:Bessel}, n) = XXX
classical_σ(::Type{<:Bessel})  = x -> x^2
classical_τ(::Type{<:Bessel{α}}) where {α} = x -> 2 + (α + 2) * x
function classical_hypergeometric(::Type{<:Bessel{α}}, n, x) where {α}
    as = (-n, n + α +1)
    bs = ()
    pFq(as, bs, -x/2)
end



(ch::Bessel{T})(x::S) where {T,S} = orthogonal_polyval(ch, x)

# Conversion from J{α, β} to J{γ,δ}
function Base.iterate(o::Connection{P, Q}, state=nothing) where
    {β, P <: Bessel{β},
     α, Q <: Bessel{α}}
    
    n,k = o.n, o.k
    if state  == nothing
        j = k
        val = connection_α(P,Q,j,k)        
    else
        j, val = state
        j += 1  # likely j += 2 if parity condition
        val = connection_α(P,Q,j,k)        # XXX can do ratio
    end
    j > n && return nothing

    (j,val), (j, val)
end

# p29 of thesis
function connection_α(::Type{<:Bessel{β}},
                      ::Type{<:Bessel{α}},
                      n,m) where {α, β}

    val = (-1)*m * (2m+β+1)
    val *= Pochhammer(-n,m) * Pochhammer(n+α+1, m)
    val *= gamma(m + β+1)*gamma(β-α+1)
    val /= gamma(1+m)*gamma(n+m+β+2) * gamma(m-n+β-α+1)
    
    val
    
end
