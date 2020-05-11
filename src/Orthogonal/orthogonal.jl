abstract type AbstractOrthogonalPolynomial{T} <: AbstractSpecialPolynomial{T} end
abstract type OrthogonalPolynomial{T} <: AbstractOrthogonalPolynomial{T} end
abstract type DiscreteOrthogonalPolynomial{T} <: AbstractOrthogonalPolynomial{T} end


"""
    AbstractCCOP{T}


A classic continuous orthogonal polynomial (CCOP) is characterized by 5 coefficients: a,b,c,d,e where

(a⋅x²+b⋅x+c)*P₍ᵢ₊₂₎'' + (d⋅x + e) * P₍ᵢ₊₁₎ + λᵢ Pᵢ = 0.

Let σ = (a⋅x²+b⋅x+c), τ  =  (d⋅x + e).

From these several structural  equations are represented,  for  example
the three-point recusion.

P₍ᵢ₊₁) = (Aᵢ⋅x + Bᵢ) * Pᵢ - Cᵢ *  P₍ᵢ₋₁₎

comes from the structural equation

x⋅p_n   = [an, bn, cn]    ⋅ [p_{n+1}, p_n, p_{n-1}]     #  Eqn (7)


Following Koepf  and  Schmersau, Representations of Orthogonal Polynomials, https://arxiv.org/pdf/math/9703217.pdf
The other structure equations are

x⋅p_n   = [an, bn, cn]    ⋅ [p_{n+1}, p_n, p_{n-1}]     #  Eqn (7)
σ⋅p'_n  = [αn, βn, γn]    ⋅  [p_{n+1}, p_n, p_{n-1}]    # Eqn (9), n ≥ 1
p_n    = [ân, b̂n, ĉn]    ⋅  [p'_{n+1}, p'_n, p'_{n-1}] # Eqn (19)
x⋅p'_n  = [αᴵn, βᴵn, γᴵn] ⋅  [p'_{n+1}, p'_n, p'_{n-1}] # Eqn  (14) with  α^*, β^*,  γ^* 

Using (7) Clenshaw evaluation using the three  point recursion is defied

Using (19) expressions for derivatives are found.

Using  (19)  expressions   for  integration are  found  (p7).

Using Thms 2,4, and 5, connection coefficients,  C(n,m) satisfying 
P_n(x) =  ∑  C(n,m)  Q_m(x) (n ≥ 0, 0 ≤  m ≤ n) are  found. These 
allow  fallback  definitions for `convert(Polynomial,p)`,  `convert(P, p::Polynomial)`,
`convert(P{α…}, p::P(β…))` and through composition  `p*q`

"""
abstract type AbstractCCOP{T,N} <: AbstractOrthogonalPolynomial{T} end
ConvertibleTypes = Union{AbstractCCOP, Polynomials.StandardBasisPolynomial}


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
leading_term(P::Type{<:AbstractCCOP},  n::Int) =  kn(P, n)

# Set defaults to be monic
kn(::Type{<:ConvertibleTypes},  n, ::Type{S}=Float64) where{S} = one(S)

# k₍ᵢ₊₁₎/kᵢ
k1k0(::Type{<:ConvertibleTypes},  i, ::Type{S}=Float64) where {S} =  one(S) #throw(MethodError())
# k₍ᵢ₊₁₎/k₍ᵢ₋₁₎
k1k_1(::Type{<:ConvertibleTypes},  i, ::Type{S}=Float64) where {S} =  one(S) #throw(MethodError())

#kn(::Type{<:Polynomials.StandardBasisPolynomial}, n::Int, ::Type{S}=Float64) where {S} = one(S) # need this
#k1k0(::Type{<:Polynomials.StandardBasisPolynomial}, i, ::Type{S}=Float64) where {S} =  one(S)


Polynomials.domain(::Type{<:AbstractCCOP}) = Polynomials.Interval(-Inf,  Inf)



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


⟒(P::Type{<:AbstractCCOP1{α}}) where {α} = constructorof(P){α}
⟒(P::Type{<:AbstractCCOP2{α,β}}) where {α, β} = constructorof(P){α,  β}


function Base.setindex!(p::AbstractCCOP, value::Number, idx::Int)
    T = eltype(p)
    ## widen size...
    idx < 0 &&  throw(ArgumentError("Negative index"))
    val = T(value)
    d = length(coeffs(p)) - 1
    if idx > d
        append!(p.coeffs, zeros(T, idx-d))
    end
    setindex!(p.coeffs, val,  idx+1)
end


# variables  are  mutable

Polynomials.degree(p::AbstractCCOP0{T,N})  where {T,N} = N-1
Polynomials.degree(p::AbstractCCOP1{α,T,N})  where {α,T,N} = N-1
Polynomials.degree(p::AbstractCCOP2{α,β,T,N})  where {α,β,T,N} = N-1
Polynomials.isconstant(p::AbstractCCOP) = degree(p) <=  0

Polynomials.zero(P::Type{<:AbstractCCOP0{T}},  var::Polynomials.SymbolLike=:x)     where {T}     = ⟒(P)(T[], var)
Polynomials.zero(P::Type{<:AbstractCCOP1{α,T}},  var::Polynomials.SymbolLike=:x)   where {α,T}   = ⟒(P)(T[], var)
Polynomials.zero(P::Type{<:AbstractCCOP2{α,β,T}},  var::Polynomials.SymbolLike=:x) where {α,β,T} = ⟒(P)(T[], var)
Polynomials.zero(P::Type{<:AbstractCCOP},  var::Polynomials.SymbolLike=:x) where {α} = zero(⟒(P){Float64}, var)
Polynomials.zero(p::P) where {P <: AbstractCCOP} = zero(P, p.var)

Polynomials.one(P::Type{<:AbstractCCOP0{T}},  var::Polynomials.SymbolLike=:x)     where {T}     = ⟒(P)(T[1], var)
Polynomials.one(P::Type{<:AbstractCCOP1{α,T}},  var::Polynomials.SymbolLike=:x)   where {α,T}   = ⟒(P)(T[1], var)
Polynomials.one(P::Type{<:AbstractCCOP2{α,β,T}},  var::Polynomials.SymbolLike=:x) where {α,β,T} = ⟒(P)(T[1], var)
Polynomials.one(P::Type{<:AbstractCCOP},  var::Polynomials.SymbolLike=:x) where {α} = one(⟒(P){Float64}, var)
Polynomials.one(p::P) where {P <: AbstractCCOP} = one(P, p.var)

# This is not right!
Polynomials.variable(P::Type{<:AbstractCCOP0{T}},  var::Polynomials.SymbolLike=:x)     where {T}     = (basis(P,1,var) - Bn(P,0,T))/An(P,0,T)
Polynomials.variable(P::Type{<:AbstractCCOP1{α,T}},  var::Polynomials.SymbolLike=:x)   where {α,T}   = (basis(P,1,var) - Bn(P,0,T))/An(P,0,T)
Polynomials.variable(P::Type{<:AbstractCCOP2{α,β,T}},  var::Polynomials.SymbolLike=:x) where {α,β,T} = (basis(P,1,var) - Bn(P,0,T))/An(P,0,T)
Polynomials.variable(P::Type{<:AbstractCCOP},  var::Polynomials.SymbolLike=:x) where {α} = variable(⟒(P){Float64}, var)
Polynomials.variable(p::P) where {P <: AbstractCCOP} = variable(P, p.var)

function Polynomials.basis(::Type{P}, n::Int, var::Polynomials.SymbolLike=:x) where {P  <: AbstractCCOP}
    T = eltype(one(P))
    cs = zeros(T, n+1)
    cs[end] = one(T)
    ⟒(P)(cs, var)
end
Polynomials.basis(p::P, n::Int) where {P <: AbstractCCOP} = basis(P, n, p.var)

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
function eval_ccop(P::Type{<:AbstractCCOP}, ::Val{N}, cs::Vector{T}, x::S) where {T,N,S}
    if @generated
        
        N == 0 && return zero(T) * zero(S)
        N == 1 && return cs[1] * one(S)
        SS = eltype(one(S))
        Δ0 = :(cs[N-1])
        Δ1 = :(cs[N])
        for i in N-1:-1:2
            a = :(cs[i - 1] - c1 * Cn(P, i-1, eltype(SS)))
            b = :(Δ0 + Δ1 * muladd(x, An(P,i-1,eltype(SS)),Bn(P,i-1,eltype(SS))))
            Δ0 = :(a)
            Δ1 = :(b)
        end
        c0 + c1* muladd(x, An(P,0,eltype(SS)), Bn(P,0,eltype(SS)))
    else
        _eval_ccop(P,cs,x)
    end
end

function _eval_ccop(P::Type{<:AbstractCCOP}, cs::Vector{T}, x::S) where {T,S}

    N = length(cs)
    N == 0 && return zero(T)*zero(S)
    N == 1 && return cs[1] * one(S)

    SS = eltype(one(S))
    Δ0 = cs[end - 1]
    Δ1 = cs[end]
    @inbounds for i in N-1:-1:2
        Δ0, Δ1 = cs[i - 1] - Δ1 * Cn(P, i-1,eltype(SS)), Δ0 + Δ1 * muladd(x, An(P,i-1,eltype(SS)),Bn(P,i-1,eltype(SS)))
    end

    return Δ0 + Δ1 * muladd(x, An(P,0,eltype(SS)),  Bn(P,0,eltype(SS)))
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


function _Bn(P::Type{<:AbstractCCOP}, a,b,c,d,e, n::Int, ::Type{S}=Float64) where {S}
    num = (2b*n*(a*n+d-a)-e*(-d+2a))
    den = (d+2a*n) * (d-2a+2a*n)

    (iszero(num) && !iszero(den)) && return zero(S)
    iszero(den) && return Bn(P, Val(n), S)  
    
    val = one(S) * num / den
    val *=  k1k0(P,n,S)
    
    val
end

function Cn(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    a,b,c,d,e = abcde(P)
    _Cn(P, a,b,c,d,e,n,S)
end

function _Cn(P::Type{<:AbstractCCOP}, a,b,c,d,e, n::Int, ::Type{S}=Float64) where {S}

    numa = (a*n+d-2a) * n * (4c*a-b^2) + 4a^2*c -a*b^2 + a*e^2 - 4*a*c*d
    numa += d*b^2 - b*e*d + d^2*c
    num = -numa * (a*n+d-2a)*n
    den = (d - 2a + 2a*n)^2 * (2a*n - 3a + d) * (2a*n - a + d)

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

# an, bn, cn
# 0 = [an,bn,cn] ⋅ [p_{n+1},p_n,p_{n-1}]
function an(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    1/An(P,n,S)
end

function bn(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    -Bn(P,n,S)/An(P,n,S)
end

function cn(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    Cn(P,n,S)/An(P,n,S)
end

# αn, βn,γn
# σ⋅pn' = [αn, βn,γn] ⋅ [p_{n+1},p_n,p_{n-1}]
function αn(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    a,b,c,d,e = abcde(P)
    num = a * n
    val = one(S) *  num
    
    val /= k1k0(P,n, S) 
    return val
end

function βn(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    a,b,c,d,e = abcde(P)
    num = -n*(a*n+d-a)*(2e*a-d*b)
    den = (d+2*a*n)(d-2a+2a*n)
    iszero(den) &&  return βn(P, Val(n), S)
    val = one(S) *  num  / den
    return val
end    


function γn(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    a,b,c,d,e = abcde(P)
    num = ((n-1)*(a*n+d-a)*(4c*a-b^2)+a*e^2 + d^2*c - b*e*d)
    num *= (a*n+d-a)*(a*n+d-2a)*n
    den = (d-2a+2a*n)^2*(2a*n-3a+d)*(2a*n-a+d)
    iszero(den) &&  return γn(P, Val(n), S)
    
    val = one(S) * num / den
    val *= k1k0(P,n-1, S)
    
    return val
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


# for integration formulas
# ân, b̂n, ĉn
# pn = [ân, b̂n, ĉn] ⋅ [p'_{n+1},p'_n,p'_{n-1}]
function ân(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    a,b,c,d,e = abcde(P)
    val = one(S)
    val /= n+1
    val /= k1k0(P,n, S)
    return val
end

function b̂n(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    a,b,c,d,e = abcde(P)
    
    num = (2e*a - d*b)
    den = (d+2a*n)*(d-2a+2a*n)
    iszero(den) && return b̂n(P, Val(n), S)
        val = one(S) * num/ den
    return val
    
end

function ĉn(P::Type{<:AbstractCCOP}, n::Int, ::Type{S}=Float64) where {S}
    a,b,c,d,e = abcde(P)
    num = ((n-1)*(a*n+d-a)*(4c*a-b^2)+a*e^2+d^2*c-b*e*d)*a*n
    den =  (d-2a+2a*n)^2*(2a*n-3a+d)*(2a*n-a+d)
    iszero(den)  &&  return ĉn(P, Val(n),  S)
    
    val = one(S)  *  num / den        
    val *= k1k0(P,n-1,S)
    return val
end

##
## --------------------------------------------------
##
## Conversion

function Base.convert(::Type{Q},  p::P)  where {Q <: Polynomials.StandardBasisPolynomial, P <: AbstractCCOP} 
    p(variable(Q, p.var))
end
function Base.convert(::Type{Q},  p::P)  where {Q <: AbstractCCOP,  P <: Polynomials.StandardBasisPolynomial}
    _convert_ccop(Q, p)
end
function Base.convert(::Type{Q}, p::P)  where  {Q <: AbstractCCOP,  P <: AbstractCCOP}
    a,b,c,d,e = abcde(P)
    ā,b̄,c̄,d̄,ē = abcde(Q)
    
    if a == ā && b == b̄ && c == c̄  #  σ = σ̄
        _convert_ccop(Q, p)
    else
        T = eltype(one(Q))
        convert(Q, convert(Polynomial{T}, p))
    end
end

##
## --------------------------------------------------
##
## multiply/additioni with P{α…,T,N} =  Q{α...,T, M}
## We don't have (p::P{N},q::P{M}) where  {N,M, P<:AbstractCCOP} as an available signature
## so we create these fall  backs, and direct + and  * to ⊕ and  ⊗  in the `register` macros
function ⊕(p::P, q::Q) where {P <: AbstractCCOP, Q <: AbstractCCOP}

    @assert  ⟒(P) == ⟒(Q)
    @assert eltype(p) == eltype(q)

    Polynomials.isconstant(p)  && return q + p[0]
    Polynomials.isconstant(q)  && return p + q[0]    
    p.var != q.var && throw(ArgumentError("Variables don't  match"))    

    d = max(degree(p), degree(q))
    as = [p[i]+q[i] for i in 0:d]
    return ⟒(P)(as, p.var)

end

function ⊗(p::P, q::Q) where {P <: AbstractCCOP, Q <: AbstractCCOP}

    @assert  ⟒(P) == ⟒(Q)
    @assert eltype(p) == eltype(q)

    Polynomials.isconstant(p)  && return q * p[0]
    Polynomials.isconstant(q)  && return p * q[0]    
    p.var != q.var && throw(ArgumentError("Variables don't  match"))    
        
    convert(⟒(P), convert(Polynomial, p) * convert(Polynomial, q))

#    R = eltype(p)
#    convert(⟒(P){R}, convert(Polynomial, p) * convert(Polynomial, q))
end


# function Base.:*(p::P, q::Q) where {P <: AbstractCCOP, Q <: AbstractCCOP}

#     Polynomials.isconstant(p) && return q * p[0]
#     Polynomials.isconstant(q) && return p * q[0]
#     p.var != q.var && throw(ArgumentError("Variables don't  match"))

#     R = promote_type(P, Q)

#     if hasmethod(iterate, (Linearization{R, Val{:diagonal}}, ))
#         linearization_product(p, q)
#     elseif hasmethod(iterate, (Linearization{R, Val{:pq}}, ))
#         linearization_product_pq(p, q)        
#     else
#         # gives poor error message if ⟒(R) not available
#         convert(⟒(R), convert(Polynomial, p) * convert(Polynomial, q))
#     end


        
    #     R = promote_type(eltype(p), eltype(q))
    #     pq = convert(Polynomial, p) * convert(Polynomial, q)
    #     convert(⟒(P){R}, pq)
    # else
        
    #     p1,q1 = promote(p, q)
    #     p1 * q1
        
    # end
#end

# use pn= [â,b̂,ĉ] ⋅ [p'_{n+1}, p'_n, p'_{n-1}] to
# find expression for p' in terms of p
function Polynomials.derivative(p::P, order::Int=1) where {P <:AbstractCCOP}

    R = eltype(one(eltype(p))/1)
    d = degree(p)
    order < 0 && throw(ArgumentError("order must  be non-negative"))
    order == 0 && return ⟒(P){R}(coeffs(p), p.var)
    d < order && return zero(⟒(P){R})
    
    as = zeros(R, d)
    ps = R.(coeffs(p))
    for n = d-1:-1:1
        a,b,c = ân(P,n,R),b̂n(P,n,R),ĉn(P,n,R)
        if !iszero(a)
            pn = ps[1+n+1]
            as[1+n] = pn/a
            ps[1+n] -= pn*b/a
            ps[1+n-1] -= pn*c/a
        end
    end
    a,b,c = ân(P,0,R),b̂n(P,0,R),ĉn(P,0,R)
    p1 = ps[1+1]
    as[1+0] = p1/a

    dp = ⟒(P)(as, p.var)
    order == 1 ? dp : derivative(dp, order - 1)
    
end

function Polynomials.integrate(p::AbstractCCOP, C::Number)

    T,S = eltype(p), typeof(C)
    R = promote_type(eltype(one(T) / 1), S)
    P = ⟒(typeof(p)){R}
    #if hasnan(p) || hasnan(C)
    #    error("XXX nan")
    #end
    n = degree(p)
    if n == 0
        return P([C, p[0]], p.var)
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
        as[1 + d + 1] += pd * ân(P, d, S)
        as[1 + d]     += pd * b̂n(P, d, S)
        if  d > 0
            as[1 + d - 1] += pd * ĉn(P, d, S)
        end
    end

    # adjust constant
    ∫p = P(as,  p.var)
    ∫p[0] = R(C) - ∫p(0)
    
    return  ∫p
end










##           H   Hₑ  Lᵅ   Bᵅ   Gᵅ  Jᵅᵝ
## scalar    ✓   ✓   ✓   ✓    ✓   ✓
## p+q       ✓   ✓   ✓   ✓    ✓   ✓
## p*q
## Pᵅ -> Pᵝ  -   -    ✓   ✓    ✓   x (unless α, β > 0 
## pᵅ -> P   x   ✓   ✓   ✓    ✓   x (unless α, β > 0 
## P -> Pᵅ   x   ✓   ✓   ✓    ✓   x (unless α, β > 0  
## p'
## ∫p        ✓   ✓   ✓   ✓    ✓   x (unless  α,β > 0)  


using  Test
function testit()
    Ps  = (Hermite,
           ChebyshevHermite,
           Laguerre{0}, Laguerre{1},
           Gegenbauer{1/2}, Gegenbauer{2},
           ChebyshevFirst, ChebyshevSecond,
           Bessel{1/2}, Bessel{1/4},
           Jacobi{1,1/2},Jacobi{1/2,2},
           Jacobi{1/2,1/2},
           Jacobi{-1/2,1/2},
           Jacobi{1/2,-1/2}, Jacobi{-1/2,-1/2}
           )

    x = variable()
    @show :arithmetic
    for P in Ps
        # p(x), -p,p+s,p*s,p/s, p+q, p*q
        as,bs = [1,2,3], [1,2,0]
        asc = [2,2,3]
        s,c,p,q = 1, P(1,:θ), P(as), P(bs)
        p(1/2)
        p(x)
        @test -p == P(-as)
        @test p+s == P(asc)
        @test p*s == P(as .*
                       s)
        @test p/s == P(as ./ s)
        @test p+c == P(asc)
        @test p*c == P(as.*s)
        @test p+q == P(as+bs)
        #p*q
    end

    @show :integrate
    for P in Ps
        p = basis(P,3)
        out = integrate(p)(x) ≈ integrate(p(x))
        !out && @show  P
    end

    @show :conversion_within
    for P in Ps
        for Q in Ps
            a,b,c,d,e = abcde(P)
            ā,b̄,c̄,d̄,ē = abcde(Q)
    
            if a == ā && b == b̄ && c == c̄
                if P != Hermite && Q != Hermite
                    @show P,Q
                    p = P(rand(1:4, 5))
                    @test _convert_ccop(P, _convert_ccop(Q, p)) ≈ p
                end
            end
        end
    end

    
    for  (P,Q) in  (
        (Hermite, Hermite),
        (Laguerre{1/2}, Laguerre{3/2}),
        (Gegenbauer{1/4}, Gegenbauer{3/4}),
        (Bessel{1/4},  Bessel{3/4}),
        (Jacobi{1/2,1/2}, Jacobi{3/4, 3/4}),
        (Jacobi{1/2,1/2}, Jacobi{-1/2, -1/2}),
                    )
        @show  P,Q
        for i in 2:6
            p,q = basis(P,i), basis(Q,i)
            @test convert(P,q)(x) ≈ q(x)
            @test convert(Q,p)(x) ≈ p(x)
        end
    end

    
    @show :conversion_Polynomial
    Q = Polynomial
    for P in Ps
        P ∈ (Hermite, # exclude Hermite, as ĉs are an issue
             #ChebyshevHermite,
             #Laguerre{0}, Laguerre{1},
             #Gegenbauer{1/2}, Gegenbauer{2},
             #Bessel{1}, Bessel{2},
             #Jacobi{1,1/2},
             #Jacobi{1/2,2},
             #Jacobi{1/2,1/2},
             #Jacobi{-1/2,1/2},
             #Jacobi{1/2,-1/2},
             #Jacobi{-1/2,-1/2},
             ) && continue
        @show P, Q
        p = P([1,2,3,4,5])
        q = Q([1,2,3,4,5])
        @test convert(Q, p) ≈ p(x) || hasnan(p(x))
        @test convert(P, q)(x) ≈ q
    end

    @show :Jacobi
    for i in 1:5
        α, β = 2*rand(2) .- 1
        P = Jacobi{α,β}
        Q = Polynomial
        x = Q(:x)
        
        for n in 3:6
            @show n
            p,q = basis(P,n), x^n
            @test convert(Q, convert(P, q)) ≈ q
            @test convert(P, convert(Q, p)) ≈ p            
            @test convert(Q, p) ≈ p(x)
            @test convert(P, q)(x) ≈ q
        end
    end

    @show :Jacobi_12
    for (α,β) = ((1/2, 1/2),
                 (1/2, -1/2),
                 (-1/2, 1/2),
                 (-1/2, -1/2)
                 )
        @show α,β

        P = Jacobi{α,β}
        Q = Polynomial
        x = Q(:x)
        
        for n in 3:6
            p,q = basis(P,n),x^n
            @test convert(Q, convert(P, q)) ≈ q
            @test convert(P, convert(Q, p)) ≈ p            
            @test convert(Q, p) ≈ p(x) || hasnan(p(x))
            @test convert(P, q)(x) ≈ q
        end
    end
    
end






