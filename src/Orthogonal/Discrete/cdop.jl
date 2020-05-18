# generic classical discrete orthogonal polynomial, CDOP
"""
     AbstractCDOP{T,N}

Following [Koepf  and Schmersau](https://arxiv.org/pdf/math/9703217.pdf), a family `y(x)=p_n(x)=k_x⋅x^n +  ...`  
for  `n  ∈  {0, 1,…}, k_n ≠ 0` of polynomials is a family of classic *discrete* orthogonal polynomials if it  is  a
solution of a differential equation

(a⋅x²+b⋅x+c) ⋅ Δ∇y + (d⋅x + e) ⋅ ∇' + λᵢ⋅ y = 0,

where  `Δy(x) = y(x+1) - y(x)` and `∇y(x) = y(x) - y(x-1)`.

A family is characterized by the 5 coefficients: a,b,c,d,e.
Let σ = (a⋅x²+b⋅x+c), τ = (d⋅x + e).

As in the classical coninuous orthogonal polynomial case
[`AbstractCCOP`](@ref), from these 5 values the cofficients in the
there-point recursion, and other structural equations can be
represented. These allow polynomial multiplication, integration,
differentiation, conversion, etc. to be defined generically.

"""
abstract type AbstractCDOP{T,N} <: AbstractCOP{T,N} end

# subtypes  to keep track of number of parameters
# passed to  @registerN macros
abstract type AbstractCDOP0{T,N} <: AbstractCDOP{T,N} end
abstract type AbstractCDOP1{α,T,N} <: AbstractCDOP{T,N} end
abstract type AbstractCDOP2{α, β,T,N} <: AbstractCDOP{T,N}  end
abstract type AbstractCDOP3{α, β,γ,T,N} <: AbstractCDOP{T,N}  end

# compose with FallingFactorial
function Base.convert(::Type{Q}, p::P)  where  {Q <: AbstractCDOP,  P <: AbstractCDOP}
    T = eltype(P)
    _convert_ccop(Q, _convert_ccop(FallingFactorial{T},  p))
end
    
    
⟒(P::Type{<:AbstractCDOP1{α}}) where {α} = constructorof(P){α}
⟒(P::Type{<:AbstractCDOP2{α,β}}) where {α, β} = constructorof(P){α, β}
⟒(P::Type{<:AbstractCDOP3{α,β,γ}}) where {α, β, γ} = constructorof(P){α, β, γ}


Δₓ(p::AbstractCDOP) = p(variable(p)+1)-p(variable(p))
∇ₓ(p::AbstractCDOP) = p(variable(p))-p(variable(p)-1)

## Structural Equations
## https://dlmf.nist.gov/18.22
## we have σ⋅Δ∇p + τ⋅∇p + λp =  0
## Using parameterization   of  above, τ = C, σ = A - C
function _An(P::Type{<:AbstractCDOP}, a,b,c,d,e, n::Int)
    one(eltype(P))
end

function _Bn(P::Type{<:AbstractCDOP}, a,b,c,d,e, n::Int)

    S = eltype(P)
    
    num = n*(d + 2b) * (d + a*n - a) + e * (d - 2a)
    den = (2a*n - 2a   + d )  *  (d + 2a*n)

    iszero(den) && return Bn(P, Val(n))  
    
    val = one(S) * num / den
    
    val
end

function _Cn(P::Type{<:AbstractCDOP}, a,b,c,d,e, n::Int)

    S = eltype(P)
    
    num = one(S)
    num *= (n - 1) * (d + a*n - a)
    num *= a*n*d - d*b - a*d + a^2*n^2 - 2a^2*n + 4c*a + a^2 + 2e*a - b^2
    num += -d*b*e + d^2*c + a*e^2
    num *= (a*n +  d - 2a) * n
    den = (d - a + 2a*n) * (d + 2a*n - 3a) * (2a*n - 2a + d)^2

    iszero(den) && return Cn(P, Val(n))
    
    val = -one(S) * num  / den
    
    val
end


# αn, βn,γn
# σ⋅pn' = [αn, βn,γn] ⋅ [p_{n+1},p_n,p_{n-1}]
function αn(P::Type{<:AbstractCDOP}, n::Int)
    a,b,c,d,e = abcde(P)
    S = eltype(P)

    num = a * n
    val = one(S) *  num
    
    val /= k1k0(P,n) 
    return val
end

function βn(P::Type{<:AbstractCDOP}, n::Int)
    a,b,c,d,e = abcde(P)
    S = eltype(P)

    num = -n *  (d  + a*n -  2) * (2*a*n*d - a*d - d*b + 2e*a - 2a^2*n +  2a^2*n^2)
    den = (2a*n - 2a +  d) * (d + 2a*n)

    iszero(den) &&  return βn(P, Val(n))

    val = one(S) *  num  / den

    return val
end    


function γn(P::Type{<:AbstractCDOP}, n::Int)
    a,b,c,d,e = abcde(P)
    S = eltype(P)
    
    num = (n - 1)  * (d + a*n  -a)
    num *= a*n*d - d*b - a*d + a^2*n^2 - 2a^2*n + 4c*a + a^2 + 2e*a - b^2
    num += -d*b*e + d^2*c + a*e^2
    num *= (d + a*n -  a) *  (a*n + d - 2a)
    den = (d - a + 2n) * (d + 2a*n - 3a) * (2a*n - 2a + d)^2
    iszero(den) &&  return γn(P, Val(n))
    
    val = one(S) * num / den
    val *= k1k0(P,n-1)
    
    return val
end


function abcdeᴵ(P::Type{<:AbstractCDOP}) 
    a,b,c,d,e = abcde(P).a, abcde(P).b, abcde(P).c, abcde(P).d, abcde(P).e
    NamedTuple{(:a,:b,:c,:d,:e)}((a,b,c,d+2a,d + e+ a + b))
end


function ⊗(p::P, q::Q) where {P <: AbstractCDOP, Q <: AbstractCDOP}

    #@assert  ⟒(P) == ⟒(Q)
    #@assert eltype(p) == eltype(q)

    Polynomials.isconstant(p)  && return q * p[0]
    Polynomials.isconstant(q)  && return p * q[0]    
    p.var != q.var && throw(ArgumentError("Variables don't  match"))    
        
    convert(⟒(P), convert(FallingFactorial, p) * convert(FallingFactorial, q))

#    R = eltype(p)
#    convert(⟒(P){R}, convert(Polynomial, p) * convert(Polynomial, q))
end
