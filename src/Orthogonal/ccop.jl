##
## --------------------------------------------------
##
## Classic Continuous Orthogonal Polynomials

abstract type AbstractCCOPBasis <: AbstractCOPBasis end

"""
    AbstractCCOPBasis <:  AbstractCOPBasis

Following [Koepf and Schmersau](https://arxiv.org/pdf/math/9703217.pdf), a family `y(x)=p_n(x)=k_x⋅x^n +  ...`
for  `n  ∈  {0, 1,…}, k_n ≠ 0` of polynomials is a family of classic *continuous* orthogonal polynomials if each is  a
solution of the differential equation

`(a⋅x²+b⋅x+c) ⋅ yᵢ'' + (d⋅x + e) ⋅ yᵢ' + λᵢ⋅ yᵢ = 0.`

A family is characterized, up to choice of leading term, by the 5 coefficients: `a,b,c,d,e`.
Let `σ = (a⋅x²+b⋅x+c)`, `τ = (d⋅x + e)`.

From these  5  coefficients several structural  equations are represented.
For example the three-point recursion.

`P₍ᵢ₊₁₎ = (Aᵢ⋅x + Bᵢ) * Pᵢ - Cᵢ *  P₍ᵢ₋₁₎`,

where `Aᵢ,Bᵢ,Cᵢ` can be represented in formulas involving just  `a,b,c,d,e` and `i`.

Rearranging   gives the structural equation:

`x⋅p_n   = [an, bn, cn] ⋅ [p_{n+1}, p_n, p_{n-1}]`  (Eqn (7))


The other structural equations are (equation  references are from Koepf and Schmersau):

`σ⋅p'_n  = [αn, βn, γn] ⋅  [p_{n+1}, p_n, p_{n-1}]` (Eqn (9), n ≥ 1)

`p_n = [ân, b̂n, ĉn]  ⋅  [p'_{n+1}, p'_n, p'_{n-1}]` (Eqn (19))

`x⋅p'_n  = [αᴵn, βᴵn, γᴵn] ⋅  [p'_{n+1}, p'_n, p'_{n-1}]` (Eqn  (14))

Using (7), Clenshaw polynomial evaluation using the three  point recursion is defined.

Using (19), expressions for derivatives are found.

Using  (19), expressions for integration are found (p7).

Using their theorems 2,4, and 5, connection coefficients, `C(n,m)` satisfying
`P_n(x) =  ∑  C(n,m)  Q_m(x) (n ≥ 0, 0 ≤  m ≤ n)` are  found. These
allow  fallback  definitions for `convert(Polynomial,p)`,  `convert(P, p::Polynomial)`,
`convert(P{α…}, p::P(β…))` and through composition polynomial  multiplication,  `p*q`.


If non-monic versions are desired, then the  leading  term can be  specified through `kn()` (which by default is defined by the  method `k1k0(P,i)`, the ratio of  `kᵢ₊₁/kᵢ`).  The `@register_monic` macro is useful  for creating  monic versions through  method delegation from the common non-monic systems. Similarly, the `@register_shifted` macro is useful  to provide shifted versions (cf. [`ShiftedLegendre`](@ref)).

Registering a system, defining an `abcde` method, and optionally
defining `k1k0` is usually sufficient to define a new system, though
the general equations may need specializations when algebraic
cancellation is required.

The defaults for evaluation and multplication are generally an order
of magnitude slower than a directly defined function. For some
families this is done (e.g. `Chebyshev`,`ChebyshevU`, `Hermite`, `Laguerre`), but not all.

[Koekoek and Swarttouw](https://arxiv.org/pdf/math/9602214.pdf)
present an encyclopedia of formula characterizing families of
orthogonal polynomials.


"""
AbstractCCOPBasis
# type for dispatch
const AbstractCCOPPolynomial = AbstractUnivariatePolynomial{<:AbstractCCOPBasis,T,X} where {T,X}

# for conversion to base case
# \upstigma[tab]
ϛ(P::Type{<:AbstractUnivariatePolynomial{<:AbstractCCOPBasis}}) = Polynomial

## Display
# show parameters in constructor's name
# function Base.show(io::IO, mimetype::MIME"text/plain", p::P) where {α,P<:AbstractCCOP1{α}}
#     print(io, "$(P.name){$(α)}(")
#     printpoly(io, p, mimetype)
#     print(io, ")")
# end
# function Base.show(
#     io::IO,
#     mimetype::MIME"text/plain",
#     p::P,
# ) where {α,β,P<:AbstractCCOP2{α,β}}
#     print(io, "$(P.name){$(α),$(β)}(")
#     printpoly(io, p, mimetype)
#     print(io, ")")
# end
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
An(B::Type{<:AbstractCOPBasis}, n::Int) = Ãn(B, n) * k1k0(B, n)
function Ãn(B::Type{<:AbstractCOPBasis}, n::Int)
    a, b, c, d, e = abcde(B)
    Ãn(B, a, b, c, d, e, n)
end

Bn(B::Type{<:AbstractCOPBasis}, n::Int) = B̃n(B, n) * k1k0(B, n)
function B̃n(B::Type{<:AbstractCOPBasis}, n::Int)
    a, b, c, d, e = abcde(B)
    B̃n(B, a, b, c, d, e, n)
end

Cn(B::Type{<:AbstractCOPBasis}, n::Int) = C̃n(B, n) * k1k_1(B, n)
function C̃n(B::Type{<:AbstractCOPBasis}, n::Int)
    a, b, c, d, e = abcde(B)
    C̃n(B, a, b, c, d, e, n)
end
C̃n(B::Type{<:AbstractCOPBasis}, n::Val) = throw(ArgumentError("not defined"))

# may be overridden to speed up eval_cop()
function ABCₙ(B::Type{<:AbstractCOPBasis}, n::Int)
    a, b = An(B,n), Bn(B,n)
    n > 0 && return (A=a, B=b, C=Cn(B,n))
    n <= 0 && return (A=a, B=b, C=0*b)
end
## ---> CCOPBasis now

function Ãn(B::Type{<:AbstractCCOPBasis}, a, b, c, d, e, n::Int)
    1
    #one(eltype(B))
end

function B̃n(B::Type{<:AbstractCCOPBasis}, a, b, c, d, e, n::Int)
    #S = eltype(B)

    num = (2b * n * (a * n + d - a) - e * (-d + 2a))
    den = (d + 2a * n) * (d - 2a + 2a * n)

    iszero(den) && return B̃n(B, Val(n))

    val =  num / den

    val
end

function C̃n(B::Type{<:AbstractCCOPBasis}, a, b, c, d, e, n::Int)
    S = Int # eltype(B)

    numa =
        (a * n + d - 2a) * n * (4c * a - b^2) + 4a^2 * c - a * b^2 + a * e^2 - 4 * a * c * d
    numa += d * b^2 - b * e * d + d^2 * c
    num = -numa * (a * n + d - 2a) * n
    den = (d - 2a + 2a * n)^2 * (2a * n - 3a + d) * (2a * n - a + d)

    iszero(den) && return C̃n(B, Val(n))

    val = ( num) / den

    val
end

# an, bn, cn
# x⋅pn = [an,bn,cn] ⋅ [p_{n+1},p_n,p_{n-1}]
function an(B::Type{<:AbstractCOPBasis}, n::Int)
    1 / An(B, n)
end

function bn(B::Type{<:AbstractCOPBasis}, n::Int)
    -Bn(B, n) / An(B, n)
end

function cn(B::Type{<:AbstractCOPBasis}, n::Int)
    Cn(B, n) / An(B, n)
end

# αn, βn,γn
# σ⋅pn' = [αn, βn,γn] ⋅ [p_{n+1},p_n,p_{n-1}]
αn(B::Type{<:AbstractCOPBasis}, n::Int) = α̃n(B, n) / k1k0(B, n)
function α̃n(B::Type{<:AbstractCCOPBasis}, n::Int)
    a, b, c, d, e = abcde(B)
    #    S = eltype(B)

    val = (1 * a) * n

    return val
end

βn(B::Type{<:AbstractCOPBasis}, n::Int) = β̃n(B, n)
function β̃n(B::Type{<:AbstractCCOPBasis}, n::Int)
    a, b, c, d, e = abcde(B)
    S = Int # eltype(P)

    num = -n * (a * n + d - a) * (2e * a - d * b)
    den = (d + 2 * a * n) * (d - 2a + 2a * n)

    iszero(den) && return β̃n(B, n)

    val = (one(S) * num) / den

    return val
end

γn(B::Type{<:AbstractCOPBasis}, n::Int) = γ̃n(B, n) * k1k0(B, n - 1)
function γ̃n(B::Type{<:AbstractCCOPBasis}, n::Int)
    a, b, c, d, e = abcde(B)
    S = Int ##eltype(P)

    num = ((n - 1) * (a * n + d - a) * (4c * a - b^2) + a * e^2 + d^2 * c - b * e * d)
    num *= (a * n + d - a) * (a * n + d - 2a) * n
    den = (d - 2a + 2a * n)^2 * (2a * n - 3a + d) * (2a * n - a + d)
    iszero(den) && return γ̃n(B, n)

    val = (one(S) * num) / den

    return val
end

# for integration formulas
# ân, b̂n, ĉn
# pn = [ân, b̂n, ĉn] ⋅ [p'_{n+1},p'_n,p'_{n-1}]
ân(B::Type{<:AbstractCOPBasis}, n::Int) = ẫn(B, n) / k1k0(B, n)
function ẫn(B::Type{<:AbstractCCOPBasis}, n::Int)
    a, b, c, d, e = abcde(B)
    S = Int # XXeltype(B)

    val = one(S)
    val /= n + 1
    return val
end

b̂n(B::Type{<:AbstractCOPBasis}, n::Int) = b̂̃n(B, n)
function b̂̃n(B::Type{<:AbstractCCOPBasis}, n::Int)
    a, b, c, d, e = abcde(B)
    S = Int #eltype(B)

    num = (2e * a - d * b)
    den = (d + 2a * n) * (d - 2a + 2a * n)
    iszero(den) && return b̂̃n(B, Val(n))

    val = one(S) * num / den
    return val
end
b̂̃n(B::Type{<:AbstractCCOPBasis}, n::Val) = throw(ArgumentError("Not defined"))

ĉn(B::Type{<:AbstractCOPBasis}, n::Int) = ĉ̃n(B, n) * k1k0(B, n - 1)
function ĉ̃n(B::Type{<:AbstractCCOPBasis}, n::Int)
    a, b, c, d, e = abcde(B)
    S = Int #eltype(B)

    num =
        ((n - 1) * (a * n + d - a) * (4c * a - b^2) + a * e^2 + d^2 * c - b * e * d) * a * n
    den = (d - 2a + 2a * n) * (d - 2a + 2a * n) * (2a * n - 3a + d) * (2a * n - a + d)

    iszero(den) && return ĉ̃n(B, Val(n))

    val = (one(S) * num) / den
    return val
end
ĉ̃n(B::Type{<:AbstractCCOPBasis}, n::Val) = throw(ArgumentError("Not defined"))

"""
    abcdeᴵ

If P_n is a CCOP, then

`x pᵢ = αᵢ p₍ᵢ₊₁₎ + βᵢ pᵢ + γᵢ p₍ᵢ₋₁₎`.

`Sp` will `pᵢ'` with this  choice of `a`,`b`,`d`,`e` derived from  those of  `pᵢ`. (Eqn  (13))
"""
function abcdeᴵ(B::Type{<:AbstractCCOPBasis})
    a, b, c, d, e = abcde(B).a, abcde(B).b, abcde(B).c, abcde(B).d, abcde(B).e
    (a=a, b=b, c=c, d=d + 2a, e=e + b)
end

# αᴵn, βᴵn, γᴵn  (α^*,...)
# x⋅pn' = [αᴵn, βᴵn, γᴵn] ⋅  [p'_{n+1},p'_n,p'_{n-1}]
αᴵn(B::Type{<:AbstractCOPBasis}, n::Int) = α̃ᴵn(B, n) / k1k0(B, n)
function α̃ᴵn(B::Type{<:AbstractCOPBasis}, n::Int)
    n = n - 1
    a, b, c, d, e = abcdeᴵ(B)
    Aᴵn = Ãn(B, a, b, c, d, e, n) * (n + 2) / (n + 1)#*k1k0(P,n+1) # . *  k1k0(dP,n) = k(dP,n+1)/k(dP,n) = (n+2)kn(P,n+1)/((n+1)kn(P,n+1) = (n+2)/(n+1)*k1k0(P,n+1)
    1 / Aᴵn
end

βᴵn(B::Type{<:AbstractCOPBasis}, n::Int) = β̃ᴵn(B, n)
function β̃ᴵn(B::Type{<:AbstractCOPBasis}, n::Int)
    n = n - 1
    a, b, c, d, e = abcdeᴵ(B)
    Ãᴵn = Ãn(B, a, b, c, d, e, n)
    B̃ᴵn = B̃n(B, a, b, c, d, e, n)

    -B̃ᴵn / Ãᴵn
end

γᴵn(B::Type{<:AbstractCOPBasis}, n::Int) = γ̃ᴵn(B, n) * k1k0(B, n - 1)
function γ̃ᴵn(B::Type{<:AbstractCOPBasis}, n::Int)
    n = n - 1
    a, b, c, d, e = abcdeᴵ(B)
    Ãᴵn = Ãn(B, a, b, c, d, e, n)  # * k(dP,n+1)/k(dP,n)
    C̃ᴵn = C̃n(B, a, b, c, d, e, n)  # * k(dP,n+1)/k(dP,n-1)
    C̃ᴵn / Ãᴵn * (n + 1) / n #* k1k0(P, n)     # k(dp,n)/k(dp,n-1) = (n+1)k(P,n+1)/(n k(P,n)) = (n+1)/n * k1k0(n)
end

##
##  --------------------------------------------------
##
# delegate P{B} -> B for tests
# preferred usage is B (not P)
abcde(::Type{P}) where {B<:AbstractCOPBasis, P<:AbstractUnivariatePolynomial{B}} = abcde(B)
abcdeᴵ(::Type{P}) where {B<:AbstractCOPBasis, P<:AbstractUnivariatePolynomial{B}} = abcdeᴵ(B)
for fn ∈ (:An, :an, :αn, :ân,:αᴵn,
          :Bn, :bn, :βn, :b̂n,:βᴵn,
          :Cn, :cn, :γn, :ĉn, :γᴵn,
          :ABCₙ, :k0, :kn, :k1k0, :k1k_1,
          :monic,
          :gauss_nodes_weights, :weight_function, :innerproduct)
    @eval begin
        $(fn)(::Type{P},i::Int) where {B<:AbstractCOPBasis, P<:AbstractUnivariatePolynomial{B}} = $(fn)(B,i)
    end
end

##
## --------------------------------------------------
##
## Conversion

function Base.convert(
    ::Type{Q},
    p::P,
) where {Q<:Polynomials.StandardBasisPolynomial,
         P<:AbstractCOPBasis}
    X = Polynomials.indeterminate(Q, p)
    T = eltype(Q)
    x = variable(⟒(Q){T,X})
    p(x)
end

function Base.convert(
    ::Type{Q},
    p::P,
) where {Q<:AbstractOrthogonalPolynomial{<:AbstractCCOPBasis},
         P<:Polynomials.StandardBasisPolynomial}
    _convert_cop(Q, p)
end

## Conversion
## * use FastTransforms, when available for T <: AbstractFloat; see connection.jl
## * use  _convert_cop when possible (needs to match σ)
## * use  conversion  through Polynomial type
function Base.convert(::Type{Q}, p::P) where {
    Q<:AbstractOrthogonalPolynomial{<:AbstractCCOPBasis},
    P<:AbstractOrthogonalPolynomial{<:AbstractCCOPBasis}}
    _convert(Q, p)
end

# work around method ambiguity introduced in abstract
# dispatch  to  specific FastTransform  method  (defined in `connection.jl`) or
# use this default
function _convert(::Type{Q}, p::P) where {
    Q<:AbstractOrthogonalPolynomial{<:AbstractCCOPBasis},
    P<:AbstractOrthogonalPolynomial{<:AbstractCCOPBasis}}
    a, b, c, d, e = abcde(basistype(P))
    ā, b̄, c̄, d̄, ē = abcde(basistype(Q))

    same_σ = a == ā && b == b̄ && c == c̄

    if same_σ  && !(Q <: Hermite || P <: Hermite)  # σ ==  σ̄   and not Hermite
        _convert_cop(Q, p)
    else
        T = eltype(Q)
        convert(Q, convert(Polynomial{T}, p))
    end
end

##
## --------------------------------------------------
# # scalar ops

# Polynomials.scalar_add (c,p)
#function Base.:+(p::P, c::S) where {B<:AbstractCCOPBasis, T,X,P<:AbstractUnivariatePolynomial{B,T,X}, S<:Number}
function Polynomials.scalar_add(c::S, p::P) where {B<:AbstractCCOPBasis, T,X,P<:AbstractUnivariatePolynomial{B,T,X}, S<:Number}
    c′ = c / k0(B) #one(T) * ⟒(P)(c)[0] #  c / k0(P), needed to add to a coefficient
    R = promote_type(T, typeof(c′))
    N = length(p)
    iszero(c) && return (N == 0 ? zero(⟒(P){R,X}) : ⟒(P){R,X}(R.(p.coeffs)))
    N == 0 && return ⟒(P){R,X}(R(c))
    N == 1 && return ⟒(P){R,X}(R[p[0] + c′])
    cs = R[iszero(i) ? p[0] + c′ : p[i] for i in 0:(N - 1)]
    return ⟒(P){R,X}(Val(false), cs)
end

# ##
# ## multiply/addition/divrem with P{α…,T} =  Q{α...,T, M}
# ## We don't have (p::P{N},q::P{M}) where  {N,M, P<:AbstractCOP} as an available signature
# ## so we create these fall  backs, and direct +,  *, divrem to ⊕, ⊗, _divrrem  in the `register` macros
# function ⊕(p::P, q::Q) where {T,X,S,Y,M, P <: AbstractCOP{T,X}, Q <: AbstractCOP{S,Y,M}}

#     #@assert  ⟒(P) == ⟒(Q)
#     #@assert eltype(p) == eltype(q)
#     R = promote_type(T,S)
#     Polynomials.isconstant(p)  && return q + p(0)
#     Polynomials.isconstant(q)  && return p + q(0)

#     X != Y && throw(ArgumentError("Variables don't  match"))

#     if N==M
#         as = [p[i]+q[i] for i in 0:N-1]
#         return ⟒(P){R,X}(as) # might have trailing zeros here
#     end

#     N > M ? ⊕′(p,q) : ⊕′(q,p)

# end

# # assumes N >  M
# @generated  function ⊕′(p1::P, p2::Q) where {T,X,S,Y,M, P<:AbstractCOP{T,X}, Q<:AbstractCOP{S,Y,M}}

#     #@assert  ⟒(P) == ⟒(Q)
#     R = promote_type(T,S)

#     exprs = Any[nothing for i = 1:N]
#     for i in  1:M
#         exprs[i] = :(p1.coeffs[$i] + p2.coeffs[$i])
#     end
#     for i in M+1:N
#         exprs[i] =:(p1.coeffs[$i])
#     end

#     return quote
#         Base.@_inline_meta
#         ⟒(P){$(R),X,$(N)}(tuple($(exprs...)))
#     end

# end

# multiplication intercepts P{B,S,X} x P{B,T,X}
for P ∈ Polynomials.ZeroBasedDensePolynomialContainerTypes
    @eval begin
        function Base.:*(p::P, q::Q) where {B <: AbstractCCOPBasis,X,
                                            T, P<:$P{B,T,X},
                                            S, Q<:$P{B,S,X}}
            p′, q′ = _convert_cop.(Polynomial, (p, q))
            _convert_cop(⟒(P), p′ * q′)
        end
    end
end

# function ⊗(p::P, q::Q) where {P<:AbstractCOPPolynomial,Q<:AbstractCOPPolynomial}
#     error("⊗")
#     isconstant(p) && return q * constantterm(p)
#     isconstant(q) && return p * constantterm(q)
#     assert_same_variable(p, q) || throw(ArgumentError("`p` and `q` have different indeterminate"))

#     # use connection for linearization;  note:  evalauation  is  faster than _convert_cop
#     p′, q′ = _convert_cop.(Polynomial, (p, q))
#     _convert_cop(⟒(P), p′ * q′)
#     #    convert(⟒(P), convert(Polynomial, p) * convert(Polynomial, q))

# end
# use pn= [â,b̂,ĉ] ⋅ [p'_{n+1}, p'_n, p'_{n-1}] to
# find expression for p' in terms of p
function Polynomials.derivative(p::P) where {B<:AbstractCOPBasis,T,X, P<:AbstractUnivariatePolynomial{B,T,X}}
    R = eltype(one(eltype(p)) / 1)
    d = degree(p)
    d ≤ 0 && return zero(⟒(P){R,X})

    as = zeros(R, d)
    ps = R.(coeffs(p))
    for n in (d - 1):-1:1
        a, b, c = ân(B, n), b̂n(B, n), ĉn(B, n)
        if !iszero(a)
            pn = ps[1 + n + 1]
            as[1 + n] = pn / a
            ps[1 + n] -= pn * b / a
            ps[1 + n - 1] -= pn * c / a
        end
    end
    a, b = ân(B, 0), b̂n(B, 0)
    p1 = ps[1 + 1]
    as[1 + 0] = p1 / a

    dp = ⟒(P){R,X}(as)
    dp
end

function Polynomials.integrate(p::P) where {B<:AbstractCOPBasis, T, X, P<:AbstractUnivariatePolynomial{B,T,X}}

    R = typeof(one(T) / 1)
    Q = ⟒(P){R,X}
    hasnan(p) &&  return Q(NaN)

    n = degree(p)
    if n == -1
        return zero(Q)
    elseif n == 0
        return p(0) * variable(Q)
    end

    as = zeros(R, n + 2)

    # case d=0 we do by hand,  as
    # P_0(x) = c, so ∫P_o = x = c*variable(P)
    c₀, c₁ = coeffs(variable(p))
    pd = first(p.coeffs)
    as[1] = pd * c₀
    as[2] = pd * (c₁ * k0(B))
    @inbounds for d in 1:n
        pd            = p.coeffs[d + 1]
        as[1 + d + 1] += pd * ân(B, d)
        as[1 + d]     += pd * b̂n(B, d)
        if d > 1
            as[1 + d - 1] += pd * ĉn(B, d)
        end
    end

    ∫p = Q(as)
    return ∫p
end



##
## ------------------------------------------------
##

## Use EXPLICIT BARYCENTRIC WEIGHTS FOR POLYNOMIAL INTERPOLATION
## IN THE ROOTS OR EXTREMA OF CLASSICAL ORTHOGONAL POLYNOMIALS,
## By: HAIYONG WANG, DAAN HUYBRECHS, AND STEFAN VANDEWALLE
## MATHEMATICS OF COMPUTATION
## Volume 83, Number 290, November 2014, Pages 2893–2914 S 0025-5718(2014)02821-4
## https://www.jstor.org/stable/24488682

function lagrange_barycentric_nodes_weights(P::Type{<:PP}, n::Int) where {B<:AbstractCCOPBasis, PP<:AbstractUnivariatePolynomial{B}}
    ## formula (2.15) simplified as we use ratios in Lagrange
    xs, λs = gauss_nodes_weights(B, n+1)
    a, b, c, d, e = abcde(P)
    ws = [sqrt((a*xᵢ^2 + b*xᵢ + c) * λᵢ) for (xᵢ, λᵢ) ∈ zip(xs, λs)]
    N = length(ws)

    itr = isodd(n) ? (2:2:N) : (1:2:N)
    for i ∈ itr
        ws[i] = -ws[i]
    end

    xs, ws
end
