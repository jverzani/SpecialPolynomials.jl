module CCOP

using Polynomials
using SpecialFunctions
const Γ = gamma
"""
    AbstractCCOP{T}


A classic continuous orthogonal polynomial (CCOP) is characterized by 5 coefficients: a,b,c,d,e where

(a⋅x²+b⋅x+c)*P₍ᵢ₊₂₎'' + (d⋅x + e) * P₍ᵢ₊₁₎ + λᵢ Pᵢ = 0.

Let σ = (a⋅x²+b⋅x+c), τ  =  (d⋅x + e).

From these several structural  equations are represented.
The three point recursion, expressed as:

P₍ᵢ₊₁) = (Aᵢ⋅x + Bᵢ) * Pᵢ - Cᵢ *  P₍ᵢ₋₁₎

The three point recursion  is  utilized by  Clenshaws method to  evaluate the polynomials.

Following Koepf  and  Schmersau, Representations of Orthogonal Polynomials, https://arxiv.org/pdf/math/9703217.pdf
The other structure equations are

x⋅p_n   = [an, bn, cn]    ⋅ [p_{n+1}, p_n, p_{n-1}]     #  Eqn (7)
σ⋅p'_n  = [αn, βn, γn]    ⋅  [p_{n+1}, p_n, p_{n-1}]    # Eqn (9), n ≥ 1
p_n    = [ân, b̂n, ĉn]    ⋅  [p'_{n+1}, p'_n, p'_{n-1}] # Eqn (19)
x⋅p'_n  = [αᴵn, βᴵn, γᴵn] ⋅  [p'_{n+1}, p'_n, p'_{n-1}] # Eqn  (14) with  α^*, β^*,  γ^* 

Using  (19)  expressions   for  integration are  found  (p7).
Connection coefficients,  C(n,m) satisfying P_n(x) =  ∑  C(n,m)  Q_m(x) (n ≥ 0, 0 ≤  m ≤ n) are  found.

"""
abstract type AbstractCCOP{T,N} <: Polynomials.AbstractPolynomial{T} end


##
## -----
##
## interface for a  given type

basis_symbol(::Type{<:AbstractCCOP}) = "P"

"""
   abcde

A tuple returning  the  constants a,b,c,d,e  for a CCOP type with
(a⋅x²+b⋅x+c)*P₍ᵢ₊₂₎'' + (d⋅x + e) * P₍ᵢ₊₁₎ + λᵢ Pᵢ = 0.
"""
abcde(::Type{<:AbstractCCOP}) = throw(MethodError())

# kn is the leading term (Section 3 table)
kn(::Type{<:AbstractCCOP},  n::Int) =  throw(MethodError())
leading_term(P::Type{<:AbstractCCOP},  n::Int) =  kn(P, n)

# k₍ᵢ₊₁₎/kᵢ
k1k0(::Type{<:AbstractCCOP},  i::Int) =  throw(MethodError())
# k₍ᵢ₊₁₎/k₍ᵢ₋₁₎
k1k_1(::Type{<:AbstractCCOP},  i::Int) =  throw(MethodError())

Polynomials.domain(::Type{<:AbstractCCOP}) = Polynomials.Interval(-Inf,  Inf)



abstract type AbstractCCOP0{T,N} <: AbstractCCOP{T,N} end
abstract type AbstractCCOP1{α,T,N} <: AbstractCCOP{T,N} end
abstract type AbstractCCOP2{α, β,T,N} <: AbstractCCOP{T,N}  end
## constructor of,  but keeps the parameters (strips off T,  N)
## remove once not  a module
@generated function constructorof(::Type{T}) where T
    getfield(parentmodule(T), nameof(T))
end

⟒(P::Type{<:AbstractCCOP}) = constructorof(P)
⟒(P::Type{<:AbstractCCOP1{α}}) where {α} = constructorof(P){α}
⟒(P::Type{<:AbstractCCOP2{α,β}}) where {α, β} = constructorof(P){α,  β}


function Base.setindex!(p::AbstractCCOP, value::Number, idx::Int)
    throw(ArgumentError("CCOPs are immutable"))
end

function Polynomials.showterm(io::IO, ::Type{P}, pj::T, var, j, first::Bool, mimetype) where {N, T, P <: AbstractCCOP}
    iszero(pj) && return false
    !first &&  print(io, " ")
    print(io, Polynomials.hasneg(T)  && Polynomials.isneg(pj) ? "- " :  (!first ? "+ " : ""))
    print(io, "$(abs(pj))⋅$(basis_symbol(P))_$j($(var[]))")
    return true
end



# variables  are  mutable
var(p::AbstractCCOP) = p.var[]
var!(p::AbstractCCOP, x::Symbol) = p.var[] = Symbol(x)

Polynomials.degree(p::AbstractCCOP0{T,N})  where {T,N} = N-1
Polynomials.degree(p::AbstractCCOP1{α,T,N})  where {α,T,N} = N-1
Polynomials.degree(p::AbstractCCOP2{α,β,T,N})  where {α,β,T,N} = N-1
Polynomials.isconstant(p::AbstractCCOP) = degree(p) <=  0

Polynomials.zero(P::Type{<:AbstractCCOP0{T}},  var::Polynomials.SymbolLike=:x)     where {T}     = ⟒(P)(NTuple{0,T}(), Symbol(var))
Polynomials.zero(P::Type{<:AbstractCCOP1{α,T}},  var::Polynomials.SymbolLike=:x)   where {α,T}   = ⟒(P)(NTuple{0,T}(), Symbol(var))
Polynomials.zero(P::Type{<:AbstractCCOP2{α,β,T}},  var::Polynomials.SymbolLike=:x) where {α,β,T} = ⟒(P)(NTuple{0,T}(), Symbol(var))
Polynomials.zero(P::Type{<:AbstractCCOP},  var::Polynomials.SymbolLike=:x) where {α} = zero(⟒(P){Float64}, Symbol(var))

Polynomials.variable(P::Type{<:AbstractCCOP0{T}},  var::Polynomials.SymbolLike=:x)     where {T}     = ⟒(P)((zero(T),one(T)), Symbol(var))
Polynomials.variable(P::Type{<:AbstractCCOP1{α,T}},  var::Polynomials.SymbolLike=:x)   where {α,T}   = ⟒(P)((zero(T),one(T)), Symbol(var))
Polynomials.variable(P::Type{<:AbstractCCOP2{α,β,T}},  var::Polynomials.SymbolLike=:x) where {α,β,T} = ⟒(P)((zero(T),one(T)), Symbol(var))
Polynomials.variable(P::Type{<:AbstractCCOP},  var::Polynomials.SymbolLike=:x) where {α} = zero(⟒(P){Float64}, Symbol(var))

basis(P::Type{<:AbstractCCOP{T}}, n::Int, var::Polynomials.SymbolLike=:x)  where {T} =
    ⟒(P)(NTuple{n+1,T}(i==n+1 ? one(T) : zero(T) for i in 1:n+1), Symbol(var))
basis(P::Type{<:AbstractCCOP}, n::Int, var::Polynomials.SymbolLike=:x) = basis(⟒(P){Float64}, n, Symbol(var))

##
## -----
##

"""
    abcdeᴵ

If P_n is a CCOP, then 

x pᵢ = αᵢ p₍ᵢ₊₁₎ + βᵢ pᵢ + γᵢ p₍ᵢ₋₁₎.

Sp will pᵢ' with this  choice of a,b,d,e derived from  those of  pᵢ. (Eqn  (13)
"""
function abcdeᴵ(P::Type{<:AbstractCCOP}) 
    a,b,c,d,e = abcde(P).a, abcde(P).b, abcde(P).c, abcde(P).d, abcde(P).e
    NamedTuple{(:a,:b,:c,:d,:e)}((a,b,c,d+2a,e+b))
end


"""
    Clenshaw evaluation of an orthogonal polynomial 
"""
function eval_ccop(P::Type{<:AbstractCCOP}, cs::NTuple{N,T}, x::S) where {T,N,S}
    if @generated
        
        N == 0 && return zero(T) * zero(S)
        N == 1 && return cs[1] * one(S)

        Δ0 = :(cs[N-1])
        Δ1 = :(cs[N])
        for i in N-1:-1:2
            a = :(cs[i - 1] - c1 * Cn(P, i-1, eltype(S)))
            b = :(Δ0 + Δ1 * muladd(x, An(P,i-1,eltype(S)),Bn(P,i-1,eltype(S))))
            Δ0 = :(a)
            Δ1 = :(b)
        end
        c0 + c1* muladd(x, An(P,0,eltype(S)), Bn(P,0,eltype(S)))
    else
        _eval_ccop(P,cs,x)
    end
end

function _eval_ccop(P::Type{<:AbstractCCOP}, cs::NTuple{N,T}, x::S) where {N,T,S}

    N == 0 && return zero(T)*zero(S)
    N == 1 && return cs[1] * one(S)
    
    Δ0 = cs[end - 1]
    Δ1 = cs[end]
    @inbounds for i in N-1:-1:2
        Δ0, Δ1 = cs[i - 1] - Δ1 * Cn(P, i-1,eltype(S)), Δ0 + Δ1 * muladd(x, An(P,i-1,eltype(S)),Bn(P,i-1,eltype(S)))
    end

    return Δ0 + Δ1 * muladd(x, An(P,0,eltype(S)),  Bn(P,0,eltype(S)))
end
    

##
## Structural  equations
##
# An, Bn,  Cn
# p_{n+1} = (An*x + Bn)⋅p_n + Cn⋅p_{n-1}
function  An(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    a,b,c,d,e = abcde(P)
    _An(P, a,b,c,d,e ,n, S)
end
            
@generated function _An(P::Type{<:AbstractCCOP}, a,b,c,d,e, n::Int, ::Type{S}=Float64) where {S}
    quote
        k1k0(P, n, S)
    end
end

function Bn(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    a,b,c,d,e = abcde(P)
    _Bn(P, a,b,c,d,e ,n, S)
end

@generated function _Bn(P::Type{<:AbstractCCOP}, a,b,c,d,e, n::Int, ::Type{S}=Float64) where {S}
    
            
    quote
        
        num = (2b*n*(a*n+d-a)-e*(-d+2a))
        den = (d+2a*n) * (d-2a+2a*n)

        iszero(den) && return Bn(P, Val(n), S)  
        
        val = one(S) * num / den
        val *=  k1k0(P,n,S)

        val
    end
end

function Cn(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    a,b,c,d,e = abcde(P)
    _Cn(P, a,b,c,d,e,n,S)
end

@generated function _Cn(P::Type{<:AbstractCCOP}, a,b,c,d,e, n::Int, ::Type{S}=Float64) where {S}
    quote
        num = -((a*n+d-2a)*n*(4c*a-b^2) + 4a^2*c -a*b^2+a*e^2-4*a*c*d+d*b^2-b*e*d+d^2*c)*(a*n+d-2a)*n
        den = (d-2a+2a*n)^2*(2a*n-3a+d)*(2a*n-a+d)
        iszero(den) && return Cn(P, Val(n), S)  
        val = one(S) * num  / den
        val *= k1k_1(P, n, S)

        # oops, this is the discrete case
#        val *= -((n-1)*(d+a*n-a)*(a*n*d-d*b-a*d+a^2*n^2-2a^2*n+4c*a+a^2+2e*a-b^2)-d*b*e+d^2*c+a*e^2)
#        val *= (a*n+d-2a)*n
#        val /=  (d-a+2a*n)*(d+2a*n-3a)*(2a*n-2a+d)^2
#        val *= k1k_1(P, n, S)

        val
    end
end

# an, bn, cn
# 0 = [an,bn,cn] ⋅ [p_{n+1},p_n,p_{n-1}]
@generated function an(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    quote
        1/An(P,n,S)
    end
end

@generated function bn(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    quote
        -Bn(P,n,S)/An(P,n,S)
    end
end

@generated function cn(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    quote
        Cn(P,n,S)/An(P,n,S)
    end
end

# αn, βn,γn
# σ⋅pn' = [αn, βn,γn] ⋅ [p_{n+1},p_n,p_{n-1}]
@generated function αn(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    quote
        a,b,c,d,e = abcde(P)
        num = a * n
        val = one(S) *  num
        
        val /= k1k0(P,n, S) 
        return val
    end
end

@generated function βn(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    quote
        a,b,c,d,e = abcde(P)
        num = -n*(a*n+d-a)*(2e*a-d*b)
        den = (d+2*a*n)(d-2a+2a*n)
        iszero(den) &&  return βn(P, Val(n), S)
        val = one(S) *  num  / den
        return val

    end
end

@generated function γn(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    quote
        a,b,c,d,e = abcde(P)
        num = ((n-1)*(a*n+d-a)*(4c*a-b^2)+a*e^2 + d^2*c - b*e*d)
        num *= (a*n+d-a)*(a*n+d-2a)*n
        den = (d-2a+2a*n)^2*(2a*n-3a+d)*(2a*n-a+d)
        iszero(den) &&  return γn(P, Val(n), S)

        val = one(S) * num / den
        val *= k1k0(P,n-1, S)

        return val
    end
end



# αᴵn, βᴵn, γᴵn  (α^*,...)
# x⋅pn' = [αᴵn, βᴵn, γᴵn] ⋅  [p'_{n+1},p'_n,p'_{n-1}]
function αᴵn(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    a,b,c,d,e = abcdeᴵ(P)
    1/_An(P,a,b,c,d,e,n,S)
end

function βᴵn(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    a,b,c,d,e = abcdeᴵ(P)
    -(_Bn(P,a,b,c,d,e,n,S)/An(P,a,b,c,d,e,n,S))
end

function γᴵn(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    a,b,c,d,e = abcdeᴵ(P)
    _Cn(P,a,b,c,d,e,n,S)/_An(P,a,b,c,d,e,n,S)
end


# for ∫pn = ...
# ân, b̂n, ĉn
# pn = [ân, b̂n, ĉn] ⋅ [p'_{n+1},p'_n,p'_{n-1}]
@generated function ân(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    quote
        a,b,c,d,e = abcde(P)
        val = one(S)
        val /= n+1
        val /= k1k0(P,n, S)
        return val
    end
end

@generated function b̂n(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    quote
        a,b,c,d,e = abcde(P)

        num = (2e*a - d*b)
        den = (d+2a*n)*(d-2a+2a*n)
        iszero(den) && return b̂n(P, Val(n), S)
        val = one(S) * num/ den
        return val

    end
end

@generated function ĉn(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    quote
        a,b,c,d,e = abcde(P)
        num = ((n-1)*(a*n+d-a)*(4c*a-b^2)+a*e^2+d^2*c-b*e*d)*a*n
        den =  (d-2a+2a*n)^2*(2a*n-3a+d)*(2a*n-a+d)
        iszero(den)  &&  return ĉn(P, Val(n),  S)

        val = one(S)  *  num / den        
        val *= k1k0(P,n-1,S)
        return val
    end
end

# conversions convert(P, q::Q) gtrhough connection and iteration


# polynomial operations -p, p+c, p*c, p/c, p+q, p*q
function Base.:+(p::AbstractCCOP, c::S) where {S<:Number}
    P,T,N = typeof(p),eltype(p), degree(p)+1
    R = promote_type(T,S)
    as = NTuple{N,R}(i == 1 ? ai + c : ai for (i,ai) in enumerate(p.coeffs))
    ⟒(P)(as, var(p))
end

function Base.:-(p::AbstractCCOP)
    P,T,N = typeof(p),eltype(p), degree(p)+1
    as = NTuple{N,T}(-ai for ai in p.coeffs)
    ⟒(P)(as, var(p))
end

function Base.:*(p::AbstractCCOP, c::S) where {S<:Number}
    P,T,N = typeof(p),eltype(p), degree(p)+1
    R = promote_type(T,S)
    as = NTuple{N,R}(ai*c for ai in p.coeffs)
    ⟒(P)(as, var(p))
end

function Base.:/(p::AbstractCCOP, c::S) where {S<:Number}
    P,T,N = typeof(p),eltype(p), degree(p)+1
    R = eltype(one(T)/one(S))
    as = NTuple{N,R}(ai/c for ai in p.coeffs)
    ⟒(P)(as, var(p))
end

function Base.:+(p::P, q::Q) where {P <: AbstractCCOP, Q <: AbstractCCOP}

    Polynomials.isconstant(p) && return q + p[0]
    Polynomials.isconstant(q) && return p + q[0]
    var(p) != var(q) && throw(ArgumentError("Variables don't  match"))
    
    if ⟒(P) == ⟒(Q)
        
        T,N,S,M = eltype(p),degree(p)+1, eltype(q), degree(q)+1
        R = promote_type(T,S)
        NM = max(N,M)
        as = NTuple{NM,R}(p[i]+q[i] for i in 0:NM-1)
        return  ⟒(P)(as, var(p))
       
    else
        
        p1,q1 = promote(p1, q1)
        p1 + q1
        
    end
end

function Base.:*(p::P, q::Q) where {P <: AbstractCCOP, Q <: AbstractCCOP}

    Polynomials.isconstant(p) && return q * p[0]
    Polynomials.isconstant(q) && return p * q[0]
    var(p) != var(q) && throw(ArgumentError("Variables don't  match"))
    
    if ⟒(P) == ⟒(Q)
        # XXX
       
    else
        
        p1,q1 = promote(p1, q1)
        p1 * q1
        
    end
end

# Connection coefficients (Thm 2  for  case with σ = σ̄)
# can write as an iterator, if we keep state
function connection(P::Type{<:AbstractCCOP}, Q::Type{<:AbstractCCOP}, m, n)
    m > n && return 0.0
    
    a,b,c,d,e = abcde(P)
    ā,b̄,c̄,d̄,ē = abcde(Q)

    c2 = -(d +2a*n)^2 * (d-a+2a*n) * (d+2a*n+2a) * (d + a + 2a*n) * (-m + n + 2)
    c2 *= (d̄ + a*m  +  a  +a*n)

    c1 = (d - a + 2n)*(d+a+2a*n)*(n+2)*(d+2a*n)
    c1 *= (-2ē*a*d - b*d^2*n -2*m*a^2*e + 2m^2*a^2*e + 2a^2*e*n - 4ē*a^2*n^2 + 2d̄*b*d - 2d̄*e*a
           + m*b*d*a - b*d^2 + 2e*d*a -  4ē*a*d*n + 2*e*d*a*n - a*b*d*n^2 -b*d*n*a + 2a^2*e*n^2 -  m^2*a*b*d
           -4ē*a^2*n - d̄*m*b*d + 2d̄*a*n*b + 2d̄*a*n^2*b +2d̄*d*b*n + 2d̄*m*e*a
           -ē*d^2 + d̄*d*e)

    c0 = (d+2a*n+2a) * (n+2)  * (n+1)  *  (a*n - a*m + d -  d̄)
    c0 *= (a*n +a*m - a + d) *  (b*e*d -a*e^2 -d^2*c -4*a*c*n*d -4a^2*c*n^2 + a*b^2*n^2 + n*b^2 * d)

    (c0, c1, c2)
end
    

# can write as an iterator, if we keep state
function _connection2(P::Type{<:AbstractCCOP}, Q::Type{<:AbstractCCOP}, m, n)
    m > n && return 0.0
    
    a,b,c,d,e = abcde(P)
    ā,b̄,c̄,d̄,ē = abcde(Q)

    c0 = -(m - n)*(a * m + d - a + a * n) * (d̄ + 2a * m) * (d̄ + a + 2a * m) * (d̄ + 3a + 2a * m)
    c0 *= (d̄ + 2a * m + 2a)^2

    c1 = -d * b * n * d̄ + 2d * a * m^2 * b + d * b * d̄ + 2d * a * m * b + 2d * ē * n * a 
    c1 += d * d̄ * e + 2d * d̄ * b * m - m * b * d^2 - e * d̄^2 -  4a^2 * m^2 * e - m^2*a*b*d̄ +b*n*d̄*a - 2e*d̄*a
    c1 -= 4*a^2*m*e - 4*e*d̄*a*m + 2m^2 * a^2 * ē + 2ē * a^2 * n^2 - 2ē*a^2 * n - m * a *  b *  d̄ + 2m * d̄ * ē * a  
    c1 += 2m * ē * a^2 - b*n^2*d̄ *a
    c1 *= (d̄ +2a*m +2a) * (m+1) * (d̄ + a + 2a*m) * (d̄ + 3a + 2a*m)

    c2 = -(d̄ + 2a*m) * (m+1) * (-a*m -2a +a*n - d̄ + d) *( a*m + a*n + a +  d̄)
    c2 *= (a*b^2*m^2 - 4a^2*m^2*c -8a^2*m*c + 2a *m *b^2 - 4a*d̄*m*c + m*b*d̄^2 - 4a*d̄*c - a*ē^2 + a*b^2 -c*d̄^2
           +b * ē * d̄ - 4a^2*c + b^2*d) * (m + 2)

    
    
    (c0, c1, c2)
end

# Return C(n,m): return C(n,n), C(n,n-1), .... C(n, 0)
struct Connection2{P,Q}
    n::Int
end

function Base.iterate(o::Connection2{P,Q}, state=nothing) where {P,Q}
    n = o.n
    if state == nothing
        Cn1, Cn = 0, 1 # eltype
        m = n
        Cm = Cn
        return (m, Cm), (m, Cn1, Cn)
    else
        m, Cn1, Cn = state
        m -= 1
        m < 0 && return nothing
        c0, c1, c2 = _connection2(P,Q,m,n)
        # c0 * Cm + c1 * C_{m+1} * c2*C_{m_2} = 0
        # c0 * Cm = -(c1 * C_{m+1} + c2 * C_{m+2))
        #         = -(c1*Cn, +c2*Cn1)
        Cm = iszero(c0) ? 0 : -(c1*Cn + c2*Cn1)/c0

        return (m, Cm), (m, Cn, Cm)
    end
end
        
    
           
        
                                                                        

function Polynomials.integrate(p::AbstractCCOP, C::Number)

    T,S = eltype(p), typeof(C)
    R = promote_type(eltype(one(T) / 1), S)
    P = ⟒(typeof(p)){R}
    #if hasnan(p) || isnan(C)
    #    error("XXX nan")
    #end
    n = degree(p)
    if n == 1
        return P([C, p[0]], var(p))
    end
    
    as = zeros(R, n + 2)
    
    @inbounds for d in 0:n

        pd = p.coeffs[d+1]
        as[1 + d + 1] += pd * ân(P, d, S)
        as[1 + d]     += pd * b̂n(P, d, S)
        if  d > 0
            as[1 + d - 1] += pd * ĉn(P, d, S)
        end
    end

    # adjust constant
    F = P(as,  var(p))
    as[1] = R(C) - F(0)
    ∫p = P(as, var(p))
    
    return  ∫p
end

function Polynomials.derivative(p::AbstractCCOP, order::Int=1)
    throw(MethodError()) #XXX
end
##
## --------------------------------------------------
##
# Macros to register POLY{T},  POLY{α, T} and POLY{α, β, T}
#
# We use `NTuple{N,T}` with N =  d+1 to store  a degree d  polynomial - no trailing zeros permitted.
macro register0(name)
    poly = esc(name)
    quote
        struct $poly{T,N} <: AbstractCCOP0{T,N}
            coeffs::NTuple{N,T}
            var::Base.RefValue{Symbol}
            
            function $poly{T,N}(coeffs::NTuple{M,S},  var::Symbol=:x) where {T,N,S,M}
                lz = findlast(!iszero, coeffs)
                if lz ==  nothing
                    new{T,0}(NTuple{0,T}(), Ref(var))
                elseif lz != N
                    throw(ArgumentError("coeffs is wrong size"))
                else
                    cs = NTuple{N,T}(i <= lz ? T(coeffs[i]) : zero(T)  for  i in 1:N)
                    new{T,N}(cs,  Ref(var))
                end
            end
            
            function $poly{T}(coeffs::NTuple{M,S},  var::Symbol=:x) where {T,S,M}
                lz = findlast(!iszero, coeffs)
                if lz ==  nothing
                    new{T,0}(NTuple{0,T}(), Ref(var))
                else
                    new{T,lz}(T.(coeffs),  Ref(var))
                end
            end

            function $poly(coeffs::NTuple{N,T},  var::Symbol=:x) where {α, T, N}
                lz = findlast(!iszero,  coeffs)
                if lz == nothing
                    new{T,0}(NTuple{0,T}(), Ref(var))
                else
                    cs = NTuple{N,T}(coeffs[i] for i in  1:lz)
                    new{T,lz}(cs, Ref(var))
                end
            end

            
        end

        (p::$poly)(x::S) where  {S} = eval_ccop(typeof(p), p.coeffs, x)
    
        Base.convert(::Type{P}, p::P) where {P<:$poly} = p
        Base.promote_rule(::Type{$poly{T}}, ::Type{$poly{S}}) where {T,S} =
            $poly{promote_type(T, S)}
        function $poly{T}(x::AbstractVector{S}, var::Polynomials.SymbolLike = :x) where {T,S}
            $poly{T}(NTuple{length(x),T}(T(xi) for xi in x), Symbol(var))
        end
        $poly(coeffs::AbstractVector{T}, var::Polynomials.SymbolLike=:x) where {T} =
            $poly{T}(NTuple{length(coeffs),T}(c for c in coeffs), Symbol(var))
        $poly{T}(n::Number, var::Polynomials.SymbolLike = :x) where {T} = $poly((one(T),), Symbol(var))
        $poly(n::S, var::Polynomials.SymbolLike = :x) where {S<:Number} = $poly{S}(n,var)
        $poly{T}(var::Polynomials.SymbolLike=:x) where {T} = variable($poly{T}, Symbol(var))
        $poly(var::Polynomials.SymbolLike=:x) = variable($poly, Symbol(var))
        
    end
end

macro register1(name)
    poly = esc(name)
    quote
        struct $poly{α,T,N} <: AbstractCCOP1{α,T,N}
            coeffs::NTuple{N,T}
            var::Base.RefValue{Symbol}
            
            function $poly{α,T,N}(coeffs::NTuple{M,S},  var::Symbol=:x) where {α,T,N,S,M}
                lz = findlast(!iszero, coeffs)
                if lz ==  nothing
                    new{α,T,0}(NTuple{0,T}(), Ref(var))
                elseif lz != N
                    throw(ArgumentError("coeffs is wrong size"))
                else
                    cs = NTuple{N,T}(i <= lz ? T(coeffs[i]) : zero(T)  for  i in 1:N)
                    new{α,T,N}(cs,  Ref(var))
                end
            end
            
            function $poly{α,T}(coeffs::NTuple{M,S},  var::Symbol=:x) where {α,T,S,M}
                lz = findlast(!iszero, coeffs)
                if lz ==  nothing
                    new{α,T,0}(NTuple{0,T}(), Ref(var))
                else
                    new{α,T,lz}(T.(coeffs),  Ref(var))
                end
            end
            
            function $poly{α}(coeffs::NTuple{N,T},  var::Symbol=:x) where {α, T, N}
                lz = findlast(!iszero,  coeffs)
                if lz == nothing
                    new{α,T,0}(NTuple{0,T}(), Ref(var))
                else
                    cs = NTuple{N,T}(coeffs[i] for i in  1:lz)
                    new{α,T,lz}(cs, Ref(var))
                end
            end
        end

        (p::$poly)(x::S) where  {S} = eval_ccop(typeof(p), p.coeffs, x)
    
        Base.convert(::Type{P}, p::P) where {P<:$poly} = p
        Base.promote_rule(::Type{$poly{α,T}}, ::Type{$poly{α,S}}) where {α,T,S} =
            $poly{α,promote_type(T, S)}
        Base.promote_rule(::Type{$poly{α,T}}, ::Type{S}) where {α,T,S<:Number} = 
            $poly{α,promote_type(T,S)}
        function $poly{α,T}(x::AbstractVector{S}, var::Polynomials.SymbolLike = :x) where {α,T,S}
            $poly{α,T}(NTuple{length(x),T}(T(xi) for xi in x), Symbol(var))
        end
        $poly{α}(coeffs::AbstractVector{T}, var::Polynomials.SymbolLike=:x) where {α,T} =
            $poly{α,T}(NTuple{length(coeffs),T}(c for c in coeffs), Symbol(var))
        $poly{α,T}(n::Number, var::Polynomials.SymbolLike = :x) where {α,T} = $poly{α}((one(T),), Symbol(var))
        $poly{α}(n::S, var::Polynomials.SymbolLike = :x) where {α,S<:Number} = $poly{α,S}(n,var)
        $poly{α,T}(var::Polynomials.SymbolLike=:x) where {α, T} = variable($poly{α,T}, Symbol(var))
        $poly{α}(var::Polynomials.SymbolLike=:x) where {α} = variable($poly{α}, Symbol(var))
        
    end
end

macro register2(name)
    poly = esc(name)
    quote
        struct $poly{α,β,T,N} <: AbstractCCOP2{α,β,T,N}
            coeffs::NTuple{N,T}
            var::Base.RefValue{Symbol}
            
            function $poly{α,β,T,N}(coeffs::NTuple{M,S},  var::Symbol=:x) where {α,β,T,N,S,M}
                lz = findlast(!iszero, coeffs)
                if lz ==  nothing
                    new{α,β,T,0}(NTuple{0,T}(), Ref(var))
                elseif lz != N
                    throw(ArgumentError("coeffs is wrong size"))
                else
                    cs = NTuple{N,T}(i <= lz ? T(coeffs[i]) : zero(T)  for  i in 1:N)
                    new{α,β,T,N}(cs,  Ref(var))
                end
            end
            
            function $poly{α,β,T}(coeffs::NTuple{M,S},  var::Symbol=:x) where {α,β,T,S,M}
                lz = findlast(!iszero, coeffs)
                if lz ==  nothing
                    new{α,β,T,0}(NTuple{0,T}(), Ref(var))
                else
                    new{α,β,T,lz}(T.(coeffs),  Ref(var))
                end
            end
            
            function $poly{α,β}(coeffs::NTuple{N,T},  var::Symbol=:x) where {α,β, T, N}
                lz = findlast(!iszero,  coeffs)
                if lz == nothing
                    new{α,β,T,0}(NTuple{0,T}(), Ref(var))
                else
                    cs = NTuple{N,T}(coeffs[i] for i in  1:lz)
                    new{α,β,T,lz}(cs, Ref(var))
                end
            end
        end

        (p::$poly)(x::S) where  {S} = eval_ccop(typeof(p), p.coeffs, x)
    
        Base.convert(::Type{P}, p::P) where {P<:$poly} = p
        Base.promote_rule(::Type{$poly{α,β,T}}, ::Type{$poly{α,β,S}}) where {α,β,T,S} =
            $poly{α,β,promote_type(T, S)}
        Base.promote_rule(::Type{$poly{α,β,T}}, ::Type{S}) where {α,β,T,S<:Number} = 
            $poly{α,β,promote_type(T,S)}
        function $poly{α,β,T}(x::AbstractVector{S}, var::Polynomials.SymbolLike = :x) where {α,β,T,S}
            $poly{α,β,T}(NTuple{length(x),T}(T(xi) for xi in x), Symbol(var))
        end
        $poly{α,β}(coeffs::AbstractVector{T}, var::Polynomials.SymbolLike=:x) where {α,β,T} =
            $poly{α,β,T}(NTuple{length(coeffs),T}(c for c in coeffs), Symbol(var))
        $poly{α,β,T}(n::Number, var::Polynomials.SymbolLike = :x) where {α,β,T} = $poly{α,β}((one(T),), Symbol(var))
        $poly{α,β}(n::S, var::Polynomials.SymbolLike = :x) where {α,β,S<:Number} = $poly{α,β,S}(n,var)
        $poly{α,β,T}(var::Polynomials.SymbolLike=:x) where {α,β, T} = variable($poly{α,β,T}, Symbol(var))
        $poly{α,β}(var::Polynomials.SymbolLike=:x) where {α,β} = variable($poly{α,β}, Symbol(var))
        
    end
end

# pull value from  Val
getval(::Val{T}) where {T} = T

##
## -----
##

## abcde for each family
## Families       a   b   c   d       e
## Hermite        1   0   0  -2       0
## Hermiteₑ       0   0   1  -1       0
## Laguerre{α}    0   1   0  -1      α+1
## Bessel{α}      1   0   0  α        2
## Gegenbaeur{α} -1   0   1 -(2α+1)   0
## Jacobi{α,β}   -1   0   1 -(α+β_2) β-α
## standard

## Hermite
@register0 Hermite

basis_symbol(::Type{<:Hermite}) = "H"
abcde(::Type{<:Hermite})  = NamedTuple{(:a,:b,:c,:d,:e)}((1,0,0,-2,0))

# Use override here, as we get  0/0 in  default  defn
An(::Type{Hermite{T,N}}, n::Int, ::Type{S}) where {T,N,S} = 2
Bn(::Type{Hermite{T,N}}, n::Int, ::Type{S}) where {T,N,S} = 0
Cn(::Type{Hermite{T,N}}, n::Int, ::Type{S}) where {T,N,S} = 2n

b̂n(::Type{<:Hermite}, n::Val{M}, ::Type{S}) where {M,S} = zero(S)
ĉn(::Type{<:Hermite}, n::Val{M}, ::Type{S}) where {M,S} = zero(S)

function kn(::Type{<:Hermite}, n::Int)
    2^n
end
function k1k0(::Type{<:Hermite}, k, ::Type{S}=Float64) where {S}
    val = 2*one(S)
    val
end
function k1k_1(P::Type{<:Hermite}, k, ::Type{S}=Float64) where {S}
    @assert k > 0
    k  == 1  && return 2*one(S) #???
    val = 4*one(S)
    return val
end


## ChebyshevHermite
@register0 ChebyshevHermite

basis_symbol(::Type{<:ChebyshevHermite}) = "Hₑ"
# https://arxiv.org/pdf/1901.01648.pdf eqn 17
abcde(::Type{<:ChebyshevHermite})  = NamedTuple{(:a,:b,:c,:d,:e)}((0,0,1,-1,0))

kn(::Type{<:ChebyshevHermite}, n::Int) = 1
k1k0(::Type{<:ChebyshevHermite}, k::Int, ::Type{S}=Float64) where {S} = one(S)
k1k_1(P::Type{<:ChebyshevHermite}, k::Int, ::Type{S}=Float64) where {S} =  one(S)

## Laguerre Polynomials
@register1 Laguerre

#(p::Laguerre)(x::S) where  {S} = eval_ccop(typeof(p), p.coeffs, x)


basis_symbol(::Type{<:Laguerre{α}}) where {α} = "L^($α)"
abcde(::Type{<:Laguerre{α}})  where {α} = NamedTuple{(:a,:b,:c,:d,:e)}((0,1,0,-1,α+1))

function kn(::Type{<:Laguerre{α}}, n::Int)  where {α}
    (-1)^n/gamma(1+n)#Γ(1+n)
end
k1k0(::Type{<:Laguerre{α}}, k, ::Type{S}=Float64) where {S,α} = -one(S)/(k+1)
k1k_1(::Type{<:Laguerre{α}}, k, ::Type{S}=Float64) where {S,α} =  one(S)/(k+1)/k

## Gegenbauer Polynomials
@register1 Gegenbauer

basis_symbol(::Type{<:Gegenbauer{α}}) where {α} = "C^($α)"
abcde(::Type{<:Gegenbauer{α}})  where {α} = NamedTuple{(:a,:b,:c,:d,:e)}((-1,0,1,-(2α+1),0))

# overrides
Bn(::Type{<:Gegenbauer{α}}, ::Val{0}, ::Type{S})  where  {α,S}  =  zero(S)
b̂n(::Type{<:Gegenbauer{α}}, ::Val{0}, ::Type{S}) where {α,S} = one(S) *  4/α/(α+2) 
ĉn(::Type{<:Gegenbauer{α}}, ::Val{0}, ::Type{S}) where {α,S} = one(S) *  4/α/(α+2) 


function kn(::Type{<:Gegenbauer{α}}, n::Int) where {α}
    Pochhammer_factorial(α, n) * 2^n
end

function k1k0(::Type{<:Gegenbauer{α}}, k, ::Type{S}=Float64) where {S,α}
    iszero(k) && return 2*one(S)*α
    val = one(S)
    val *= 2(α+k)
    val /= k + 1
    val
end
function k1k_1(P::Type{<:Gegenbauer{α}}, k, ::Type{S}=Float64) where {S,α}
    @assert k > 0
    val = one(S)
    if k == 1
        val *= 2α*(α+1)
    else
        val *= 4(α+k)*(α+k-1)
        val /= (k+1)*k
    end
    return val
end

const Legendre = Gegenbauer{1/2}

## Bessel
@register1 Bessel

basis_symbol(::Type{<:Bessel{α}}) where {α} = "C^($α)"
abcde(::Type{<:Bessel{α}})  where {α} = NamedTuple{(:a,:b,:c,:d,:e)}((1,0,0,α, 2))

# From https://www.ams.org/journals/tran/1949-065-01/S0002-9947-1949-0028473-1/S0002-9947-1949-0028473-1.pdf we use
# kn = (n + α - 1)_n / 2^n not (n+α+1)_n/2^n
# Koepf suggests kn = (n + α + 1)_n/2^n
function kn(::Type{<:Bessel{α}}, n::Int) where {α}
    one(α)*Pochhammer(n+α-1,n)/2^n
end
function k1k0(::Type{<:Bessel{α}}, k, ::Type{S}=Float64) where {S,α}
    iszero(k) && return α/2
    val = one(S)
    val *=  (2k+α)*(2k+α-1)
    val /= (k+α-1)*2
    val
end
function k1k_1(P::Type{<:Bessel{α}}, k, ::Type{S}=Float64) where {S,α}
    @assert k > 0
    val = one(S)
    if k == 1
        val *= (α + 1)*(α + 2)
        val /= 4
    else
        val *= (2k+α-3) * (2k+α-2) * (2k+α-1) * (2k+α)
        val /= 4 * (k+α-2) * (k+α-1)
    end
    return val
end

## Jacobi Polynomials
@register2 Jacobi

basis_symbol(::Type{<:Jacobi{α,β}}) where {α,β} = "Jᵅᵝ"
abcde(::Type{<:Jacobi{α,β}})  where {α,β} = NamedTuple{(:a,:b,:c,:d,:e)}((-1,0,1,-(α+β+2),β-α))

function kn(::Type{<:Jacobi{α,β}}, n::Int)  where {α,β}
    (one(α)*one(β) * generalized_binomial(2n+α+β, n))/2^n
end
# kn =  1/2^n ⋅ choose(2n+α+β, n) = (2n+α+β)_n / (2^b  n!)
function k1k0(P::Type{<:Jacobi{α,β}}, n, ::Type{S}=Float64) where {S,α,β}
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
        num = (α+β + 4) *(α+β + 6) *(α+β + 6)
        den = 24 * (α + β + 2)
    else
        num = (γ -1) * (γ - 0) * (γ + 1) * (γ + 2)
        den = 4 * n * (n+1) * (n + α + β) * (n + α + β + 1)
    end
    
    iszero(den) && return k1k_1(P, Val(n), S)

    val *= num /den
    return val
end

end


