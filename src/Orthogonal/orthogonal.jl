## Abstract types
##
abstract type AbstractOrthogonalPolynomial{T} <: AbstractSpecialPolynomial{T} end
abstract type AbstractContinuousOrthogonalPolynomial{T} <: AbstractOrthogonalPolynomial{T} end
abstract type AbstractDiscreteOrthogonalPolynomial{T} <: AbstractOrthogonalPolynomial{T} end
abstract type AbstractCOP{T,N} <: AbstractOrthogonalPolynomial{T} end



"""
    AbstractCCOP{T,N}

Following [Koepf  and Schmersau](https://arxiv.org/pdf/math/9703217.pdf), a family `y(x)=p_n(x)=k_x⋅x^n +  ...`  
for  `n  ∈  {0, 1,…}, k_n ≠ 0` of polynomials is a family of classic *continuous* orthogonal polynomials if it  is  a
solution of a differential equation

(a⋅x²+b⋅x+c) ⋅ y'' + (d⋅x + e) ⋅ y' + λᵢ⋅ y = 0.

A family is characterized by the 5 coefficients: a,b,c,d,e.
Let σ = (a⋅x²+b⋅x+c), τ = (d⋅x + e).

From these  5  coefficients several structural  equations are represented. For example
the three-point recusion.

P₍ᵢ₊₁) = (Aᵢ⋅x + Bᵢ) * Pᵢ - Cᵢ *  P₍ᵢ₋₁₎

Here `Aᵢ,Bᵢ,Cᵢ`  can be represented in formulas involving just  `a,b,c,d,e` and `i`.

Rearraging   gives the structural equation

x⋅p_n   = [an, bn, cn]    ⋅ [p_{n+1}, p_n, p_{n-1}]     #  Eqn (7)


The other structural equations are (equation  references are for Koepf  and  Schmerrsaus):

σ⋅p'_n  = [αn, βn, γn]    ⋅  [p_{n+1}, p_n, p_{n-1}]    # Eqn (9), n ≥ 1
p_n    = [ân, b̂n, ĉn]    ⋅  [p'_{n+1}, p'_n, p'_{n-1}] # Eqn (19)
x⋅p'_n  = [αᴵn, βᴵn, γᴵn] ⋅  [p'_{n+1}, p'_n, p'_{n-1}] # Eqn  (14) with  α^*, β^*,  γ^* 

Using (7), Clenshaw polynomial evaluation using the three  point recursion is defined.

Using (19), expressions for derivatives are found.

Using  (19),  expressions   for  integration are  found  (p7).

Using Thms 2,4, and 5, connection coefficients,  C(n,m) satisfying 
P_n(x) =  ∑  C(n,m)  Q_m(x) (n ≥ 0, 0 ≤  m ≤ n) are  found. These 
allow  fallback  definitions for `convert(Polynomial,p)`,  `convert(P, p::Polynomial)`,
`convert(P{α…}, p::P(β…))` and through composition  `p*q`

Subtypes of `AbstractCCOP` are  created through  the  `@registerN` macros, where `N` is the number  of  parameters used to describe the family.

If non-monic versions are desired, then the  leading  term can be  specified through   `kn()`.  

This is  sufficient for many cases, though the general  equations may need specializations when algebraic cancellation is required. 

## Example

For this example, the value of `Bn` at `0` needs help:

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> const SP=SpecialPolynomials
SpecialPolynomials

julia> SP.@register0 MonicLegendre SP.AbstractCCOP0

julia> SP.abcde(::Type{<:MonicLegendre})  = (-1,0,1,-2,0)

julia> SP.Bn(P::Type{<:MonicLegendre}, ::Val{0}) =  0

julia> 𝐐  =  Rational{Int}
Rational{Int64}

julia> x = variable(Polynomial{𝐐})
Polynomial(x)

julia> [basis(MonicLegendre{𝐐}, i)(x) for i  in 0:5]
6-element Array{Polynomial{Rational{Int64}},1}:
 Polynomial(1//1)
 Polynomial(x)
 Polynomial(-1//3 + x^2)
 Polynomial(-3//5*x + x^3)
 Polynomial(3//35 - 6//7*x^2 + x^4)
 Polynomial(5//21*x - 10//9*x^3 + x^5)
```

"""
abstract type AbstractCCOP{T,N} <: AbstractCOP{T,N} end


##
## -----
##
## interface for a  given type

basis_symbol(::Type{<:AbstractCOP}) = "P"
Polynomials.domain(::Type{<:AbstractCOP}) = Polynomials.Interval(-Inf,  Inf)
Base.extrema(P::Type{<:AbstractCOP}) = (first(domain(P)), last(domain(P)))

"""
   abcde

A named tuple returning  the  constants a,b,c,d,e  for a CCOP type with
(a⋅x²+b⋅x+c)*P₍ᵢ₊₂₎'' + (d⋅x + e) * P₍ᵢ₊₁₎ + λᵢ Pᵢ = 0.
"""
abcde(::Type{<:AbstractCOP}) = throw(MethodError())

# kn is the leading term (Section 3 table)
leading_term(P::Type{<:AbstractCOP},  n::Int) =  kn(P, n)

# Set defaults to be monic
kn(::Type{P},  n) where{P <: AbstractCOP} = one(eltype(one(P))) #  need one(P), as eltype(Polynomial{T}) != T

# k₍ᵢ₊₁₎/kᵢ
k1k0(::Type{P},  i) where {P<:AbstractCOP} =  kn(P,i+1)/kn(P,i)
# k₍ᵢ₊₁₎/k₍ᵢ₋₁₎
k1k_1(::Type{P},  i) where {P<:AbstractCOP} =  kn(P,i+1)/kn(P,i-1)

monic(p::P) where {T,N,P <: AbstractCOP{T,N}} = N == 0 ? p : p/(kn(P,degree(p))*p[end])


# subtypes  to keep track of number of parameters
# passed to  @registerN macros
abstract type AbstractCCOP0{T,N} <: AbstractCCOP{T,N} end
abstract type AbstractCCOP1{α,T,N} <: AbstractCCOP{T,N} end
abstract type AbstractCCOP2{α, β,T,N} <: AbstractCCOP{T,N}  end

## Display
# show parameters in constructor's name
function Base.show(io::IO, mimetype::MIME"text/plain", p::P) where {α,P<:AbstractCCOP1{α}}
    print(io,"$(P.name){$(α)}(")
    printpoly(io, p, mimetype)
    print(io,")")
end
function Base.show(io::IO, mimetype::MIME"text/plain", p::P) where {α,β,P<:AbstractCCOP2{α,β}}
    print(io,"$(P.name){$(α),$(β)}(")
    printpoly(io, p, mimetype)
    print(io,")")
end

# We want to  be able to strip  off T,N  or α,...,T,N
# * constructorof(P{α,..,T,N}) =  P
# * ⟒(P(α,...,T,N)) =  P(α...)
⟒(P::Type{<:AbstractCCOP1{α}}) where {α} = constructorof(P){α}
⟒(P::Type{<:AbstractCCOP2{α,β}}) where {α, β} = constructorof(P){α,  β}

# Can't  change  N  here
function Base.setindex!(p::AbstractCOP{T,N}, value::Number, idx::Int) where {T, N}

    ## widen size...
    idx < 0 &&  throw(ArgumentError("Negative index"))
    val = T(value)
    d = N - 1
    if idx > d || (idx == d && iszero(value))
        throw(ArgumentError("Polynomials of AbstractCCOP type have fixed size parameter, N, which  can't  be changed through  assignment. Make new polynomial  instance?"))
    end
    setindex!(p.coeffs, val,  idx+1)
end

Polynomials.degree(p::AbstractCOP{T,N})  where {T,N} = N-1
Polynomials.isconstant(p::AbstractCOP) = degree(p) <=  0

Polynomials.zero(::Type{P},  var::Polynomials.SymbolLike=:x) where {P<:AbstractCOP} = ⟒(P)(eltype(P)[], var)
Polynomials.one(::Type{P},  var::Polynomials.SymbolLike=:x) where {P<:AbstractCOP} = ⟒(P)(ones(eltype(P),1), var)
Polynomials.variable(P::Type{<:AbstractCOP},  var::Polynomials.SymbolLike=:x) = (basis(P,1,var) - Bn(P,0)) / An(P,0)

#  XXX var is keyword here (odd)
function Polynomials.basis(::Type{P}, n::Int; var::Polynomials.SymbolLike=:x) where {P  <: AbstractCCOP}
    T = eltype(P)
    cs = zeros(T, n+1)
    cs[end] = one(T)
    ⟒(P)(cs, var)
end


##
## -----
##
##
## -----
##

## abcde for each family
## Families       a   b   c   d       e
## Hermite        1   0   0  -2       0
## Hermiteₑ       0   0   1  -1       0
## ChebyshevT    -1   0   1  -1       0
## ChebyshevU    -1   0   1  -3       0
## Laguerre{α}    0   1   0  -1      α+1
## Bessel{α}      1   0   0  α        2
## Gegenbaeur{α} -1   0   1 -2(α+1)   0
## Jacobi{α,β}   -1   0   1 -(α+β+2) β-α
## standard


"""
    Clenshaw evaluation of an orthogonal polynomial 
"""
function eval_ccop(P::Type{<:AbstractCOP{T,N}}, cs, x::S) where {T,N,S}
    if @generated
        
        N == 0 && return zero(T) * zero(S)
        N == 1 && return cs[1] * one(S)
        #SS = eltype(one(S))
        Δ0 = :(cs[N-1])
        Δ1 = :(cs[N])
        for i in N-1:-1:2
            a = :(cs[i - 1] - c1 * Cn(P, i-1))
            b = :(Δ0 + Δ1 * muladd(x, An(P,i-1),Bn(P,i-1)))
            Δ0 = :(a)
            Δ1 = :(b)
        end
        c0 + c1* muladd(x, An(P,0), Bn(P,0))
    else
        clenshaw_eval(P, cs, x)
    end
end

function clenshaw_eval(P::Type{<:AbstractCOP{T,N}}, cs, x::S) where {N, T,S}


    N == 0 && return zero(T)*zero(S)
    N == 1 && return cs[1] * one(S)

    #SS = eltype(one(S))
    Δ0 = cs[end - 1]
    Δ1 = cs[end]
    @inbounds for i in N-1:-1:2
        Δ0, Δ1 = cs[i - 1] - Δ1 * Cn(P, i-1), Δ0 + Δ1 * muladd(x, An(P,i-1),Bn(P,i-1))
    end

    return Δ0 + Δ1 * muladd(x, An(P,0),  Bn(P,0))
end
    

##
## Structural  equations
##
# An, Bn,  Cn
# p_{n+1} = (An*x + Bn)⋅p_n + Cn⋅p_{n-1}
function  An(P::Type{<:AbstractCOP}, n::Int)
    a,b,c,d,e = abcde(P)
    _An(P, a,b,c,d,e ,n) *  k1k0(P, n)
end

function Bn(P::Type{<:AbstractCOP}, n::Int)
    a,b,c,d,e = abcde(P)
    _Bn(P, a,b,c,d,e ,n) * k1k0(P,n)
end

function Cn(P::Type{<:AbstractCOP}, n::Int)
    a,b,c,d,e = abcde(P)
    _Cn(P, a,b,c,d,e,n) *  k1k_1(P, n)
end


function _An(P::Type{<:AbstractCCOP}, a,b,c,d,e, n::Int)
    one(eltype(P))
end


function _Bn(P::Type{<:AbstractCCOP}, a,b,c,d,e, n::Int)

    S = eltype(P)
    
    num = (2b*n*(a*n+d-a)-e*(-d+2a))
    den = (d+2a*n) * (d-2a+2a*n)

    iszero(den) && return Bn(P, Val(n))  
    
    val = one(S) * num / den
    
    val
end

function _Cn(P::Type{<:AbstractCCOP}, a,b,c,d,e, n::Int)

    S = eltype(P)
    
    numa = (a*n+d-2a) * n * (4c*a-b^2) + 4a^2*c -a*b^2 + a*e^2 - 4*a*c*d
    numa += d*b^2 - b*e*d + d^2*c
    num = -numa * (a*n + d - 2a) * n
    den = (d - 2a + 2a*n)^2 * (2a*n - 3a + d) * (2a*n - a + d)

    iszero(den) && return Cn(P, Val(n))
    
    val = one(S) * num  / den

    val
end

# an, bn, cn
# x⋅pn = [an,bn,cn] ⋅ [p_{n+1},p_n,p_{n-1}]
function an(P::Type{<:AbstractCOP}, n::Int)
    1/An(P,n)
end

function bn(P::Type{<:AbstractCOP}, n::Int)
    -Bn(P,n)/An(P,n)
end

function cn(P::Type{<:AbstractCOP}, n::Int)
    Cn(P,n)/An(P,n)
end

# αn, βn,γn
# σ⋅pn' = [αn, βn,γn] ⋅ [p_{n+1},p_n,p_{n-1}]
function αn(P::Type{<:AbstractCCOP}, n::Int)
    a,b,c,d,e = abcde(P)
    S = eltype(P)

    num = a * n
    val = one(S) *  num
    
    val /= k1k0(P,n) 
    return val
end

function βn(P::Type{<:AbstractCCOP}, n::Int)
    a,b,c,d,e = abcde(P)
    S = eltype(P)

    num = -n*(a*n+d-a)*(2e*a-d*b)
    den = (d+2*a*n)*(d-2a+2a*n)

    iszero(den) &&  return βn(P, Val(n))

    val = one(S) *  num  / den

    return val
end    


function γn(P::Type{<:AbstractCCOP}, n::Int)
    a,b,c,d,e = abcde(P)
    S = eltype(P)
    
    num = ((n-1) * (a*n + d - a) * (4c*a - b^2) + a*e^2 + d^2*c - b*e*d)
    num *= (a*n + d - a) * (a*n + d - 2a) * n
    den = (d - 2a + 2a*n)^2 * (2a*n - 3a + d) * (2a*n - a + d)
    iszero(den) &&  return γn(P, Val(n))
    
    val = one(S) * num / den
    val *= k1k0(P,n-1)
    
    return val
end



# for integration formulas
# ân, b̂n, ĉn
# pn = [ân, b̂n, ĉn] ⋅ [p'_{n+1},p'_n,p'_{n-1}]
function ân(P::Type{<:AbstractCCOP}, n::Int)
    a,b,c,d,e = abcde(P)
    S = eltype(P)
    
    val = one(S)
    val /= n+1
    val /= k1k0(P,n)
    return val
end

function b̂n(P::Type{<:AbstractCCOP}, n::Int)
    a,b,c,d,e = abcde(P)
    S = eltype(P)
    
    num = (2e*a - d*b)
    den = (d+2a*n)*(d-2a+2a*n)
    iszero(den) && return b̂n(P, Val(n))
        val = one(S) * num/ den
    return val
    
end

function ĉn(P::Type{<:AbstractCCOP}, n::Int)
    a,b,c,d,e = abcde(P)
    S = eltype(P)
    
    num = ((n-1)*(a*n+d-a)*(4c*a-b^2)+a*e^2+d^2*c-b*e*d)*a*n
    den =  (d-2a+2a*n) * (d-2a+2a*n) * (2a*n-3a+d) * (2a*n-a+d)
    iszero(den)  &&  return ĉn(P, Val(n))
    val = one(S)  *  num / den
    n > 0 && (val *= k1k0(P,n-1))
    return val
end

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


# αᴵn, βᴵn, γᴵn  (α^*,...)
# x⋅pn' = [αᴵn, βᴵn, γᴵn] ⋅  [p'_{n+1},p'_n,p'_{n-1}]
function αᴵn(P::Type{<:AbstractCCOP}, n::Int)
    n  = n - 1
    a,b,c,d,e = abcdeᴵ(P)
    Aᴵn = _An(P,a,b,c,d,e,n) * (n+2)/(n+1)*k1k0(P,n+1) # . *  k1k0(dP,n) = k(dP,n+1)/k(dP,n) = (n+2)kn(P,n+1)/((n+1)kn(P,n+1) = (n+2)/(n+1)*k1k0(P,n+1)
    1/Aᴵn
end

function βᴵn(P::Type{<:AbstractCCOP}, n::Int)
    n = n - 1
    a,b,c,d,e = abcdeᴵ(P)
    Aᴵn = _An(P,a,b,c,d,e,n) # * k1k0(dP,n)
    Bᴵn = _Bn(P,a,b,c,d,e,n) # * k1k0(dP,n)
   
    - Bᴵn / Aᴵn 
end

function γᴵn(P::Type{<:AbstractCCOP}, n::Int)
    n = n - 1
    a,b,c,d,e = abcdeᴵ(P)
    Aᴵn = _An(P,a,b,c,d,e,n)  # * k(dP,n+1)/k(dP,n)
    Cᴵn = _Cn(P,a,b,c,d,e,n)  # * k(dP,n+1)/k(dP,n-1)
    Cᴵn/Aᴵn * (n+1)/n * k1k0(P, n)     # k(dp,n)/k(dp,n-1) = (n+1)k(P,n+1)/(n k(P,n)) = (n+1)/n * k1k0(n)
end


##
## --------------------------------------------------
##
## Conversion

function Base.convert(::Type{Q},  p::P)  where {Q <: Polynomials.StandardBasisPolynomial, P <: AbstractCOP} 
    p(variable(Q, p.var))
end
function Base.convert(::Type{Q},  p::P)  where {Q <: AbstractCOP,  P <: Polynomials.StandardBasisPolynomial}
    _convert_ccop(Q, p)
end
function Base.convert(::Type{Q}, p::P)  where  {Q <: AbstractCCOP,  P <: AbstractCCOP}

    a,b,c,d,e = abcde(P)
    ā,b̄,c̄,d̄,ē = abcde(Q)
    
    if (a,b,c) == (ā,b̄,c̄)  &&  !(Q <: Hermite || P  <: Hermite)  # σ ==  σ̄   and not Hermite
        _convert_ccop(Q, p)
    else
        T = eltype(Q)
        convert(Q, convert(Polynomial{T}, p))
    end
end

##
## --------------------------------------------------
##
## multiply/addition/divrem with P{α…,T,N} =  Q{α...,T, M}
## We don't have (p::P{N},q::P{M}) where  {N,M, P<:AbstractCOP} as an available signature
## so we create these fall  backs, and direct +,  *, divrem to ⊕, ⊗, _divrrem  in the `register` macros
function ⊕(p::P, q::Q) where {P <: AbstractCOP, Q <: AbstractCOP}

    #@assert  ⟒(P) == ⟒(Q)
    #@assert eltype(p) == eltype(q)

    Polynomials.isconstant(p)  && return q + p[0]
    Polynomials.isconstant(q)  && return p + q[0]    
    p.var != q.var && throw(ArgumentError("Variables don't  match"))    

    d = max(degree(p), degree(q))
    as = [p[i]+q[i] for i in 0:d]
    return ⟒(P)(as, p.var)

end

function ⊗(p::P, q::Q) where {P <: AbstractCOP, Q <: AbstractCOP}

    #@assert  ⟒(P) == ⟒(Q)
    #@assert eltype(p) == eltype(q)

    Polynomials.isconstant(p)  && return q * p[0]
    Polynomials.isconstant(q)  && return p * q[0]    
    p.var != q.var && throw(ArgumentError("Variables don't  match"))    
        
    convert(⟒(P), convert(Polynomial, p) * convert(Polynomial, q))

#    R = eltype(p)
#    convert(⟒(P){R}, convert(Polynomial, p) * convert(Polynomial, q))
end

function _divrem(num::P, den::Q) where {P <: AbstractCOP, Q <: AbstractCOP}

    #@assert  ⟒(P) == ⟒(Q)
    #@assert eltype(num) == eltype(den)

    p1 = convert(Polynomial, num)
    p2 = convert(Polynomial, den)
    q,r = divrem(p1, p2)
    convert.(P, (q,r))

end

## Modifications needed due to `N` in the type parameter

function Polynomials.truncate(p::P;
                               rtol::Real = Base.rtoldefault(real(T)),
                               atol::Real = 0,) where {T,N,P<:AbstractCOP{T,N}}
    ps = coeffs(p)
    max_coeff = maximum(abs, ps)
    thresh = max_coeff * rtol + atol
    map!(c->abs(c) <= thresh ? zero(T) : c, ps,ps)
    ⟒(P){T}(ps, p.var)
end

Polynomials.truncate!(p::P;
                      rtol::Real = Base.rtoldefault(real(T)),
                      atol::Real = 0,) where {T,N,P<:AbstractCOP{T,N}} = error("`truncate!` not defined")

function Base.chop(p::P;
                   rtol::Real = Base.rtoldefault(real(T)),
                   atol::Real = 0,) where {T,N,P<:AbstractCOP{T,N}}
    
    N == 0 && return p
    i = N-1
    while i >= 0
        val = p[i]
        if !isapprox(val, zero(T); rtol = rtol, atol = atol)
            break
        end
        i -= 1
    end
    ⟒(P)(coeffs(p)[1:i+1], p.var)
end

Polynomials.chop!(p::P;
                  rtol::Real = Base.rtoldefault(real(T)),
                  atol::Real = 0,) where {T,N,P<:AbstractCOP{T,N}} = error("`chop!` not defined")


# use pn= [â,b̂,ĉ] ⋅ [p'_{n+1}, p'_n, p'_{n-1}] to
# find expression for p' in terms of p
function Polynomials.derivative(p::P, order::Integer=1) where {P <:AbstractCOP}

    R = eltype(one(eltype(p))/1)
    d = degree(p)
    order < 0 && throw(ArgumentError("order must  be non-negative"))
    order == 0 && return ⟒(P){R}(coeffs(p), p.var)
    d < order && return zero(⟒(P){R})
    
    as = zeros(R, d)
    ps = R.(coeffs(p))
    for n = d-1:-1:1
        a,b,c = ân(P,n),b̂n(P,n),ĉn(P,n)
        if !iszero(a)
            pn = ps[1+n+1]
            as[1+n] = pn/a
            ps[1+n] -= pn*b/a
            ps[1+n-1] -= pn*c/a
        end
    end
    a,b,c = ân(P,0),b̂n(P,0),ĉn(P,0)
    p1 = ps[1+1]
    as[1+0] = p1/a

    dp = ⟒(P)(as, p.var)
    order == 1 ? dp : derivative(dp, order - 1)
    
end

function Polynomials.integrate(p::P, C::Number=0) where {P <: AbstractCOP}
    
    T,S = eltype(p), typeof(C)
    R = promote_type(typeof(one(T) / 1), S)
    Q = ⟒(P){R}
    #if hasnan(p) || hasnan(C)
    #    error("XXX nan")
    #end
    n = degree(p)
    if n == 0
        return Q([C, p[0]], p.var)
    end
    
    as = zeros(R, n + 2)

    # case d=0 we do by hand,  as
    # P_0(x) = 1, so ∫P_o = x = variable(P)
    c₀,c₁ = coeffs(variable(p))
    pd = first(p.coeffs)
    as[1] = pd*c₀
    as[2] = pd*c₁
    @inbounds for d in 1:n
        pd = p.coeffs[d+1]
        as[1 + d + 1] += pd * ân(Q, d)
        as[1 + d]     += pd * b̂n(Q, d)
        if  d > 0
            as[1 + d - 1] += pd * ĉn(Q, d)
        end
    end

    # adjust constant
    ∫p = Q(as,  p.var)
    ∫p[0] = R(C) - ∫p(0)
    
    return  ∫p
end




#  Structs for connection, linearization (connection.jl)
struct Connection{P,Q}
    n::Int
    k::Int
end

struct Linearization{P,V}
    l::Int
    n::Int
    m::Int
end
