## common functionality for  COP classical orthogonal polynomials
## These are described by 5  values: a,.b.c.d.e

abstract type  AbstractCOPBasis <: AbstractOrthogonalBasis end
const AbstractCOPPolynomial   = AbstractUnivariatePolynomial{<:AbstractCOPBasis,T,X} where {T,X}

Polynomials.evalpoly(x, p::P) where {B<:AbstractCOPBasis,P<:AbstractUnivariatePolynomial{B}} = eval_cop(B, coeffs(p), x)

function Base.promote_rule(
    P::Type{<:Polynomials.AbstractPolynomial{T,X}},
    Q::Type{<:AbstractCOPPolynomial{S,X}},
) where {T,S,X}
    ϛ(Q){promote_type(T, S),X}
end

#Base.promote_rule(
#    P::Type{<:AbstractCOPPolynomial{S,X}},
#    Q::Type{<:Polynomials.AbstractPolynomial{T,X}},
#) where {T,S,X} = ϛ(P){promote_type(T, S),X}

# generic addition for polynomials
function ⊕(
    P::Type{<:AbstractCOPPolynomial},
    p::AbstractUnivariatePolynomial{B,T,X},
    q::AbstractUnivariatePolynomial{B,S,X},
) where {B<:AbstractCOPBasis,T,S,X}
    error("⊕")
    R = promote_type(T, S)
    P′ = ⟒(P){R,X}
    N,M = length(p), length(q)
    if N > M
        cs = Polynomials.:⊕(P, p.coeffs, q.coeffs)
        return P′(Val(false), cs)
    elseif N < M
        cs = Polynomials.:⊕(P, q.coeffs, p.coeffs)
        return P′(Val(false), cs)
    else
        cs = Polynomials.:⊕(P, p.coeffs, q.coeffs)
        return P′(cs)
    end
end

# # multiplication of polynomials of type P
# function ⊗(
#     P::Type{<:AbstractCOPPolynomial},
#     p::AbstractUnivariatePolynomial{B,T,X},
#     q::AbstractUnivariatePolynomial{B,S,X},
# ) where {B,T,S,X}
#     error("XXX")
#     isconstant(p) && return q * constantterm(p) # scalar mult is faster, well defined
#     isconstant(q) && return p * constantterm(q)
#     # use connection for linearization;  note:  evalauation  is  faster than _convert_cop
#     p′, q′ = _convert_cop.(Polynomial, (p, q))
#     _convert_cop(⟒(P), p′ * q′)
#     #    convert(⟒(P), convert(Polynomial, p) * convert(Polynomial, q))

# end

# function _divrem(num::P, den::Q) where {P<:AbstractCOPPolynomial,Q<:AbstractCOPPolynomial}
#     error("XXX")
#     #@assert  ⟒(P) == ⟒(Q)
#     #@assert eltype(num) == eltype(den)

#     p1 = convert(ϛ(P), num)
#     p2 = convert(ϛ(Q), den)
#     q, r = divrem(p1, p2)
#     convert.(⟒(P), (q, r))
# end

##
## -----
##
## interface for a  given type

Polynomials.basis_symbol(::Type{<:AbstractCOPPolynomial}) = "P"

Polynomials.domain(::Type{P}) where {P <:AbstractCOPPolynomial} = Polynomials.domain(basistype(P))
Polynomials.domain(::Type{<:AbstractCOPBasis}) = Polynomials.Interval(-Inf, Inf)

"""
    classical_hypergeometric(::Type{<:AbstractCOPBasis}, n, x)

Polynomials of this type may have a representation in terms of hypergeometric functions
"""
classical_hypergeometric(P::Type{<:AbstractCOPBasis}, n, x) # = throw(ArgumentError("No default method"))

"""
    abcde

A named tuple returning  the  constants a,b,c,d,e  for a CCOP type with
(a⋅x²+b⋅x+c)*P₍ᵢ₊₂₎'' + (d⋅x + e) * P₍ᵢ₊₁₎ + λᵢ Pᵢ = 0.
"""
abcde(::Type{<:AbstractCOPBasis}) # = throw(ArgumentError("No default method"))

"""
    k1k0

Let  `kᵢ` be the leading  coefficient of the  polynomial  in  the standard basis.
This  function implements  `k₍ᵢ₊₁₎/kᵢ` for `i ≥ 0`.

The values `kᵢ` and `k₍ᵢ₊₁₎/k₍ᵢ₋₁₎` can  be generated.

The default value  leads to monic polynomials.
"""
k1k0(::Type{B}, i::Int) where {B<:AbstractCOPBasis} = 1
k1k0(::Type{B}, i::Int) where {B<:Polynomials.StandardBasis} = 1

"""
      k0(::Type{P})

This is  `basis(P,0)`,  most  often  just `1`, but may be different with some normalizations, such as `Orthonormal`.
"""
k0(::Type{B}) where {B<:AbstractCOPBasis} = 1
k0(::Type{B}) where {B<:Polynomials.StandardBasis} = 1

# Set defaults to be monic
# For non-monic subtypes of `AbstractCOP` only need k1k0  to be defined, as `kn` and `k1k_1` come for free

# kn = prod(k1k0(i) for i in 0:n-1)  *  k0(P)
function kn(::Type{B}, n::Int) where {B<:AbstractCOPBasis}
    foldr(*, (k1k0(B, i) for i in 0:(n - 1)), init=float(k0(B)))
end

# k₍ᵢ₊₁₎/k₍ᵢ₋₁₎ =  (k₍ᵢ₊₁₎/kᵢ) ⋅ (kᵢ/k₍ᵢ₋₁₎)
function k1k_1(::Type{B}, i::Int) where {B<:AbstractCOPBasis}
    i > 0 || throw(ArgumentError("i needs to be positive"))
    k1k0(B, i) * k1k0(B, i - 1)
end

# leading term (Section 3 table)  is kn
leading_term(B::Type{<:AbstractCOPBasis}, n::Int) = kn(B, n)

# square root of ratio of norm2(P,n+1)/norm2(P,n)
# Let ωᵢ = √{∫ πᵢ² dw}, this is ωᵢ₊₁ ÷ ωᵢ
ω₁₀(::Type{B}, n) where {B<:AbstractCOPBasis} = sqrt(norm2(B, n + 1) / norm2(B, n))

# # Can't  change  N  here
# function Base.setindex!(p::AbstractCOPPolynomial{T,X}, value::Number, idx::Int) where {T,X}

#     ## widen size...
#     idx < 0 && throw(ArgumentError("Negative index"))
#     val = T(value)
#     N = length(p)
#     d = N - 1
#     if idx > d || (idx == d && iszero(value))
#         throw(
#             ArgumentError(
#                 "Polynomials of type `AbstractCOP` have a  fixed size parameter, N, which  can't  be changed through  assignment. Make new polynomial  instance?",
#             ),
#         )
#     end
#     setindex!(p.coeffs, val, idx + 1)
# end

#Polynomials.degree(p::AbstractCOPPolynomial)  = length(p) - 1
#Polynomials.isconstant(p::AbstractCOPPolynomial) = degree(p) <= 0

##  Evaluation

"""
    Clenshaw evaluation of an orthogonal polynomial
"""
function eval_cop(::Type{B}, cs, x::S) where {B<:AbstractCOPBasis,S}
    T, N = eltype(cs), length(cs)
    N == 0 && return zero(T) * zero(S)
    N == 1 && return (cs[1] * k0(B)) * one(S)
    _eval_cop(B, cs, x)
end

function _eval_cop(::Type{B}, cs, x::S) where {B<:AbstractCOPBasis,S}
    return clenshaw_eval(B, cs, x)
    # if @generated
    #     quote
    #         N = length(cs) # lastindex?
    #         Δ0 = cs[end - 1]
    #         Δ1 = cs[end]
    #         @inbounds for i in (Int(N)::Int - 1):-1:2
    #             Δ0, Δ1 = cs[i - 1] - Δ1 * Cn(P, i - 1),
    #             Δ0 + Δ1 * muladd(x, An(P, i - 1), Bn(P, i - 1))
    #         end
    #         p₀ = k0(P)
    #         p₁ = muladd(x, An(P, 0), Bn(P, 0)) * p₀
    #         Δ0 * p₀ + Δ1 * p₁
    #     end
    # else
    #     clenshaw_eval(P, cs, x)
    # end
end

#  evaluate  basis vector through hypergeometric formulation
(B::Basis{P})(x) where {P<:AbstractCOPPolynomial} = eval_basis(basistype(P), B.n, x)

# Evaluate a basis vector without realizing it and without using Clenshaw
eval_basis(::Type{P}, n, x) where {P<:AbstractCOPPolynomial} = classical_hypergeometric(basistype(P), n, x)

"""
    eval_hyper(::P, cs,  x::S)

Evaluate polynomial `P(cs)` by  computing value for each basis vector from its hypergeometric representation
"""
function eval_hyper(P::Type{<:AbstractCOPPolynomial}, cs, x::S) where {S}
    isempty(cs) && return zero(S)
    tot = cs[1] * one(S)
    for i in 2:length(cs)
        tot += cs[i] * Basis(P, i - 1)(x)
    end
    tot
end
