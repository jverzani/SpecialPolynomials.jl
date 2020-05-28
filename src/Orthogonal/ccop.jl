##
## --------------------------------------------------
##
## Classic Continuos Orthogonal Polynomials

"""
    AbstractCCOP{T,N}

Following [Koepf and Schmersau](https://arxiv.org/pdf/math/9703217.pdf), a family `y(x)=p_n(x)=k_x⋅x^n +  ...`  
for  `n  ∈  {0, 1,…}, k_n ≠ 0` of polynomials is a family of classic *continuous* orthogonal polynomials if each is  a
solution of the differential equation

`(a⋅x²+b⋅x+c) ⋅ y'' + (d⋅x + e) ⋅ y' + λᵢ⋅ y = 0.`

A family is characterized by the 5 coefficients: `a,b,c,d,e`.
Let `σ = (a⋅x²+b⋅x+c)`, `τ = (d⋅x + e)`.

From these  5  coefficients several structural  equations are represented. For example
the three-point recusion.

`P₍ᵢ₊₁) = (Aᵢ⋅x + Bᵢ) * Pᵢ - Cᵢ *  P₍ᵢ₋₁₎`,

where `Aᵢ,Bᵢ,Cᵢ` can be represented in formulas involving just  `a,b,c,d,e` and `i`.

Rearranging   gives the structural equation:

x⋅p_n   = [an, bn, cn] ⋅ [p_{n+1}, p_n, p_{n-1}]     #  Eqn (7)


The other structural equations are (equation  references are from Koepf and Schmersau):

σ⋅p'_n  = [αn, βn, γn]    ⋅  [p_{n+1}, p_n, p_{n-1}]    # Eqn (9), n ≥ 1
p_n    = [ân, b̂n, ĉn]    ⋅  [p'_{n+1}, p'_n, p'_{n-1}] # Eqn (19)
x⋅p'_n  = [αᴵn, βᴵn, γᴵn] ⋅  [p'_{n+1}, p'_n, p'_{n-1}] # Eqn  (14) with  α^*, β^*,  γ^* 

Using (7), Clenshaw polynomial evaluation using the three  point recursion is defined.

Using (19), expressions for derivatives are found.

Using  (19), expressions for integration are found (p7).

Using Thms 2,4, and 5, connection coefficients, `C(n,m)` satisfying 
`P_n(x) =  ∑  C(n,m)  Q_m(x) (n ≥ 0, 0 ≤  m ≤ n)` are  found. These 
allow  fallback  definitions for `convert(Polynomial,p)`,  `convert(P, p::Polynomial)`,
`convert(P{α…}, p::P(β…))` and through composition  `p*q`

Subtypes of `AbstractCCOP` are  created through  the `@register0` or  `@registerN` macros, where the  `N`  one  is used  if parameters are  needed to describe the family.

If non-monic versions are desired, then the  leading  term can be  specified through   `kn()`.  The `@register_monic` macro is useful  for creating  monic versios through  method delegation. Similarly, the `@register_shifted` macro is useful  to provide shifted versions.

Registering a system, defining an  `abcde` method, and  optionally defining `kn`, `k1k0`,  and `k1k_1` methods is
usually sufficient to define a new system, though the general  equations may need specializations when algebraic cancellation is required. 

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

[Koekoek and Swarttouw](https://arxiv.org/pdf/math/9602214.pdf)
present an encyclopedia of formula characterizing families of
orthogonal polynomials.


"""
abstract type AbstractCCOP{T,N} <: AbstractCOP{T,N} end


# subtypes  to keep track of number of parameters
# passed to  @registerN macros
abstract type AbstractCCOP0{T,N} <: AbstractCCOP{T,N} end
abstract type AbstractCCOP1{α,T,N} <: AbstractCCOP{T,N} end
abstract type AbstractCCOP2{α,β,T,N} <: AbstractCCOP{T,N}  end
abstract type AbstractCCOP3{α,β,γ,T,N} <: AbstractCCOP{T,N}  end

# We want to  be able to strip  off T,N  or α,...,T,N
# * constructorof(P{α,..,T,N}) =  P
# * ⟒(P(α,...,T,N)) =  P(α...)
⟒(P::Type{<:AbstractCCOP1{α}}) where {α} = constructorof(P){α}
⟒(P::Type{<:AbstractCCOP2{α,β}}) where {α, β} = constructorof(P){α, β}
⟒(P::Type{<:AbstractCCOP3{α,β,γ}}) where {α, β,γ} = constructorof(P){α, β, γ}


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



##
## Structural  equations
##
## Could move the AbstractCOP terms into `cop.jl`...
# An, Bn,  Cn
# p_{n+1} = (An*x + Bn)⋅p_n + Cn⋅p_{n-1}
An(P::Type{<:AbstractCOP}, n::Int) = Ãn(P,n) * k1k0(P, n)
function  Ãn(P::Type{<:AbstractCOP}, n::Int)
    a,b,c,d,e = abcde(P)
    Ãn(P, a,b,c,d,e ,n) 
end

Bn(P::Type{<:AbstractCOP}, n::Int) = B̃n(P,n) * k1k0(P,n)
function B̃n(P::Type{<:AbstractCOP}, n::Int)
    a,b,c,d,e = abcde(P)
    B̃n(P, a,b,c,d,e ,n) 
end

Cn(P::Type{<:AbstractCOP}, n::Int) = C̃n(P,n) * k1k_1(P, n)
function C̃n(P::Type{<:AbstractCOP}, n::Int)
    a,b,c,d,e = abcde(P)
    C̃n(P, a,b,c,d,e,n)
end


function Ãn(P::Type{<:AbstractCCOP}, a,b,c,d,e, n::Int)
    one(eltype(P))
end


function B̃n(P::Type{<:AbstractCCOP}, a,b,c,d,e, n::Int)

    S = eltype(P)
    
    num = (2b*n*(a*n+d-a)-e*(-d+2a))
    den = (d+2a*n) * (d-2a+2a*n)

    iszero(den) && return B̃n(P, Val(n))  
    
    val = (one(S) * num) / den
    
    val
end

function C̃n(P::Type{<:AbstractCCOP}, a,b,c,d,e, n::Int)

    S = eltype(P)
    
    numa = (a*n+d-2a) * n * (4c*a-b^2) + 4a^2*c -a*b^2 + a*e^2 - 4*a*c*d
    numa += d*b^2 - b*e*d + d^2*c
    num = -numa * (a*n + d - 2a) * n
    den = (d - 2a + 2a*n)^2 * (2a*n - 3a + d) * (2a*n - a + d)

    iszero(den) && return C̃n(P, Val(n))
    
    val = (one(S) * num)  / den

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
αn(P::Type{<:AbstractCOP}, n::Int) = α̃n(P,n) / k1k0(P,n) 
function α̃n(P::Type{<:AbstractCCOP}, n::Int)
    a,b,c,d,e = abcde(P)
    S = eltype(P)

    val = (one(S) *  a) * n
    
    return val
end

βn(P::Type{<:AbstractCOP}, n::Int) = β̃n(P, n)
function β̃n(P::Type{<:AbstractCCOP}, n::Int)
    a,b,c,d,e = abcde(P)
    S = eltype(P)

    num = -n*(a*n+d-a)*(2e*a-d*b)
    den = (d+2*a*n)*(d-2a+2a*n)

    iszero(den) &&  return β̃n(P, Val(n))

    val = (one(S) * num)  / den

    return val
end    

γn(P::Type{<:AbstractCCOP}, n::Int) = γ̃n(P,n) *  k1k0(P,n-1)
function γ̃n(P::Type{<:AbstractCCOP}, n::Int)
    a,b,c,d,e = abcde(P)
    S = eltype(P)

    num = ((n-1) * (a*n + d - a) * (4c*a - b^2) + a*e^2 + d^2*c - b*e*d)
    num *= (a*n + d - a) * (a*n + d - 2a) * n
    den = (d - 2a + 2a*n)^2 * (2a*n - 3a + d) * (2a*n - a + d)
    iszero(den) &&  return γ̃n(P, Val(n))
    
    val = (one(S) * num) / den
    
    return val
end



# for integration formulas
# ân, b̂n, ĉn
# pn = [ân, b̂n, ĉn] ⋅ [p'_{n+1},p'_n,p'_{n-1}]
ân(P::Type{<:AbstractCCOP}, n::Int) = ẫn(P, n) / k1k0(P,n)
function ẫn(P::Type{<:AbstractCCOP}, n::Int)
    a,b,c,d,e = abcde(P)
    S = eltype(P)
    
    val = one(S)
    val /= n+1
    return val
end

b̂n(P::Type{<:AbstractCCOP}, n::Int) =  b̂̃n(P,n) 
function b̂̃n(P::Type{<:AbstractCCOP}, n::Int)
    a,b,c,d,e = abcde(P)
    S = eltype(P)
    
    num = (2e*a - d*b)
    den = (d+2a*n)*(d-2a+2a*n)
    iszero(den) && return b̂̃n(P, Val(n))

    val = one(S) * num/ den
    return val
    
end

ĉn(P::Type{<:AbstractCCOP}, n::Int) = ĉ̃n(P,n) * k1k0(P,n-1)
function ĉ̃n(P::Type{<:AbstractCCOP}, n::Int)
    a,b,c,d,e = abcde(P)
    S = eltype(P)
    
    num = ((n-1) * (a*n  + d - a) * (4c*a - b^2) + a*e^2 +  d^2*c - b*e*d) * a* n
    den = (d  - 2a + 2a*n) * (d - 2a + 2a*n) * (2a*n - 3a + d) * (2a*n - a  + d)
    iszero(den)  &&  return ĉ̃n(P, Val(n))

    val = (one(S)  *  num) / den
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
αᴵn(P::Type{<:AbstractCCOP}, n::Int) = α̃ᴵn(P,n) / k1k0(P,n)
function α̃ᴵn(P::Type{<:AbstractCCOP}, n::Int)
    n  = n - 1
    a,b,c,d,e = abcdeᴵ(P)
    Aᴵn = Ãn(P,a,b,c,d,e,n) * (n+2)/(n+1)#*k1k0(P,n+1) # . *  k1k0(dP,n) = k(dP,n+1)/k(dP,n) = (n+2)kn(P,n+1)/((n+1)kn(P,n+1) = (n+2)/(n+1)*k1k0(P,n+1)
    1/Aᴵn
end

βᴵn(P::Type{<:AbstractCCOP}, n::Int) = β̃ᴵn(P::Type{<:AbstractCCOP}, n::Int)
function β̃ᴵn(P::Type{<:AbstractCCOP}, n::Int)
    n = n - 1
    a,b,c,d,e = abcdeᴵ(P)
    Ãᴵn = Ãn(P,a,b,c,d,e,n) 
    B̃ᴵn = B̃n(P,a,b,c,d,e,n) 
   
    - B̃ᴵn / Ãᴵn 
end

γᴵn(P::Type{<:AbstractCCOP}, n::Int)  = γ̃ᴵn(P,n) * k1k0(P, n-1) 
function γ̃ᴵn(P::Type{<:AbstractCCOP}, n::Int)
    n = n - 1
    a,b,c,d,e = abcdeᴵ(P)
    Ãᴵn = Ãn(P,a,b,c,d,e,n)  # * k(dP,n+1)/k(dP,n)
    C̃ᴵn = C̃n(P,a,b,c,d,e,n)  # * k(dP,n+1)/k(dP,n-1)
    C̃ᴵn / Ãᴵn * (n+1)/n #* k1k0(P, n)     # k(dp,n)/k(dp,n-1) = (n+1)k(P,n+1)/(n k(P,n)) = (n+1)/n * k1k0(n)
end

##
##  --------------------------------------------------
##

##
## --------------------------------------------------
##
## Conversion

function Base.convert(::Type{Q},  p::P)  where {Q <: Polynomials.StandardBasisPolynomial, P <: AbstractCOP} 
    p(variable(Q, p.var))
end
function Base.convert(::Type{Q},  p::P)  where {Q <: AbstractCCOP,  P <: Polynomials.StandardBasisPolynomial}
    _convert_cop(Q, p)
end
function Base.convert(::Type{Q}, p::P)  where  {Q <: AbstractCCOP,  P <: AbstractCCOP}

    a,b,c,d,e = abcde(P)
    ā,b̄,c̄,d̄,ē = abcde(Q)
    
    if (a,b,c) == (ā,b̄,c̄)  &&  !(Q <: Hermite || P  <: Hermite)  # σ ==  σ̄   and not Hermite
        _convert_cop(Q, p)
    else
        T = eltype(Q)
        convert(Q, convert(Polynomial{T}, p))
    end
end

##
## --------------------------------------------------
# # scalar ops
# function Base.:+(p::P, c::S) where {P <: AbstractCOP, S<:Number}
#     R = promote_type(eltype(p), S)
#     as = R[a for a in coeffs(p)]
#     if length(as) >= 1
#         as[1] += c  # *assume* basis(P,0) = 1
#     else
#         push!(as, c)
#     end
#     ⟒(P)(as, p.var)
# end

#  avoid dispatch when N is known
function Base.:+(p::P, c::S) where {T, N, P<:AbstractCOP{T,N}, S<:Number}
    R = promote_type(T,S)
    iszero(c) && return (N == 0 ? zero(⟒(P){R},p.var) :  ⟒(P){R,N}(R.(p.coeffs), p.var))
    N == 0 && return ⟒(P)(R[c], p.var)
    N == 1 && return ⟒(P)(R[p[0]+c], p.var)
    cs = R[iszero(i) ? p[i]+c : p[i] for i in 0:N-1]
    return ⟒(P){R,N}(cs, p.var)
end

# function Base.:*(p::P, c::S) where {T, N, P<:AbstractCOP{T,N}, S <: Number}
#     R = promote_type(T,S)
#     iszero(c)  && return zero(⟒(P){R})
#     ⟒(P){R,N}(p.coeffs .* c, p.var)
# end

function Base.:/(p::P, c::S) where {T,N,P<:AbstractCOP{T,N},S <: Number}
    R = eltype(one(T)/one(S))
    isinf(c)  && return zero(⟒(P){R})
    ⟒(P){R,N}(p.coeffs ./ c, p.var)
end


##
## multiply/addition/divrem with P{α…,T,N} =  Q{α...,T, M}
## We don't have (p::P{N},q::P{M}) where  {N,M, P<:AbstractCOP} as an available signature
## so we create these fall  backs, and direct +,  *, divrem to ⊕, ⊗, _divrrem  in the `register` macros
function ⊕(p::P, q::Q) where {T,N,S,M, P <: AbstractCOP{T,N}, Q <: AbstractCOP{S,M}}

    #@assert  ⟒(P) == ⟒(Q)
    #@assert eltype(p) == eltype(q)

    Polynomials.isconstant(p)  && return q + p[0]
    Polynomials.isconstant(q)  && return p + q[0]    
    p.var != q.var && throw(ArgumentError("Variables don't  match"))    

    if N==M
        as = [p[i]+q[i] for i in 0:N-1]
        return ⟒(P)(as, p.var) # might have trailing zeros here
    end

    N > M ? ⊕′(p,q) : ⊕′(q,p)

end

# assumes N >  M
@generated  function ⊕′(p1::P, p2::Q) where {T,N,S,M, P<:AbstractCOP{T,N}, Q<:AbstractCOP{S,M}}

    #@assert  ⟒(P) == ⟒(Q)
    R = promote_type(T,S)

    exprs = Any[nothing for i = 1:N]
    for i in  1:M
        exprs[i] = :(p1.coeffs[$i] + p2.coeffs[$i])
    end
    for i in M+1:N
        exprs[i] =:(p1.coeffs[$i])
    end

    return quote
        Base.@_inline_meta
        ⟒(P){$(R),$(N)}(tuple($(exprs...)), p1.var)
    end

end


function ⊗(p::P, q::Q) where {P <: AbstractCOP, Q <: AbstractCOP}

    Polynomials.isconstant(p)  && return q * p[0]
    Polynomials.isconstant(q)  && return p * q[0]    
    p.var != q.var && throw(ArgumentError("Variables don't  match"))    

    # use connection for linearization
    convert(⟒(P), convert(Polynomial, p) * convert(Polynomial, q))

end

function _divrem(num::P, den::Q) where {P <: AbstractCOP, Q <: AbstractCOP}

    #@assert  ⟒(P) == ⟒(Q)
    #@assert eltype(num) == eltype(den)

    p1 = convert(Polynomial, num)
    p2 = convert(Polynomial, den)
    q,r = divrem(p1, p2)
    convert.(⟒(P), (q,r))

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
    return ∫p + (R(C) - ∫p(0))

end


