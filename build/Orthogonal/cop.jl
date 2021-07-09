## common functionality for  COP classical orthogonal polynomials
## These are described by 5  values: a,.b.c.d.e

abstract type AbstractCOP{T,X,N} <: AbstractOrthogonalPolynomial{T,X} end

function Base.promote_rule(P::Type{<:Polynomials.AbstractPolynomial{T,X}},
                           Q::Type{<:AbstractCOP{S,X}}) where {T,S,X}
    ϛ(Q){promote_type(T, S), X}
end

Base.promote_rule(P::Type{<:AbstractCOP{S,X}},
                  Q::Type{<:Polynomials.AbstractPolynomial{T,X}}) where {T,S,X} =
                      ϛ(P){promote_type(T, S), X}

# generic addition for polynomials
function ⊕(P::Type{<:AbstractCOP}, p::AbstractCOP{T,X,N} , q::AbstractCOP{S,X,M}) where {T,S,X,N,M}
    R = promote_type(T,S)
    P′ = ⟒(P){R,X}
    if N > M
        cs = Polynomials.:⊕(P, p.coeffs, q.coeffs)
        return P′{N}(cs)
    elseif N < M
        cs = Polynomials.:⊕(P, q.coeffs, p.coeffs)
        return P′{M}(cs)
    else
        cs = Polynomials.:⊕(P, p.coeffs, q.coeffs)
        return P′(cs)
    end
end    

# multiplication of polynomials of type P
function ⊗(P::Type{<:AbstractCOP}, p::AbstractCOP{T,X,N} , q::AbstractCOP{S,X,M}) where {T,S,X,N,M}
    isconstant(p)  && return q * constantterm(p) # scalar mult is faster, well defined
    isconstant(q)  && return p * constantterm(q)
    # use connection for linearization;  note:  evalauation  is  faster than _convert_cop
    p′,q′ = _convert_cop.(Polynomial, (p,q))
    _convert_cop(⟒(P), p′ * q′)
#    convert(⟒(P), convert(Polynomial, p) * convert(Polynomial, q))

end


function _divrem(num::P, den::Q) where {P <: AbstractCOP, Q <: AbstractCOP}

    #@assert  ⟒(P) == ⟒(Q)
    #@assert eltype(num) == eltype(den)

    p1 = convert(ϛ(P), num)
    p2 = convert(ϛ(Q), den)
    q,r = divrem(p1, p2)
    convert.(⟒(P), (q,r))

end

##
## -----
##
## interface for a  given type

basis_symbol(::Type{<:AbstractCOP}) = "P"
Polynomials.domain(::Type{<:AbstractCOP}) = Polynomials.Interval(-Inf, Inf)
classical_hypergeometric(P::Type{<:AbstractCOP}, n, x) = throw(ArgumentError("No default method"))

"""
   abcde

A named tuple returning  the  constants a,b,c,d,e  for a CCOP type with
(a⋅x²+b⋅x+c)*P₍ᵢ₊₂₎'' + (d⋅x + e) * P₍ᵢ₊₁₎ + λᵢ Pᵢ = 0.
"""
abcde(::Type{<:AbstractCOP}) = throw(ArgumentError("No default method"))


"""
    k1k0

Let  `kᵢ` be the leading  coeffiecient of the  polynomial  in  the standard basis. 
This  function implements  `k₍ᵢ₊₁₎/kᵢ` for `i ≥ 0`.

The values `kᵢ` and `k₍ᵢ₊₁₎/k₍ᵢ₋₁₎` can  be generated.

The default value  leads to monic polynomials.
"""
k1k0(::Type{P},  i::Int) where {P<:AbstractCOP} =  one(eltype(P))
k1k0(::Type{P},  i::Int) where {P<:Polynomials.StandardBasisPolynomial} =  one(eltype(P))

"""
      k0(::Type{P})

This is  `basis(P,0)`,  most  often  just `1`, but may be different with some normalizations, such as `Orthonormal`.
"""
k0(::Type{P})  where{P <: AbstractCOP}  =  one(eltype(P))
k0(::Type{P}) where {P <: Polynomials.StandardBasisPolynomial} = one(eltype(P))

# Set defaults to be monic
# For non-monic subtypes of `AbstractCOP` only need k1k0  to be defined, as `kn` and `k1k_1` come for free

# kn = prod(k1k0(i) for i in 0:n-1)  *  k0(P)
kn(::Type{P},  n::Int) where{P <: AbstractCOP} = foldr(*, (k1k0(P,i) for i in 0:n-1), init=k0(P)) 

# k₍ᵢ₊₁₎/k₍ᵢ₋₁₎ =  (k₍ᵢ₊₁₎/kᵢ) ⋅ (kᵢ/k₍ᵢ₋₁₎)
function k1k_1(::Type{P},  i::Int) where {P<:AbstractCOP}
    @assert i > 0
    k1k0(P,i)*k1k0(P,i-1)
end

# leading term (Section 3 table)  is kn
leading_term(P::Type{<:AbstractCOP},  n::Int) =  kn(P, n)

# square root of ratio of norm2(P,n+1)/norm2(P,n)
# Let ωᵢ = √{∫ πᵢ² dw}, this is ωᵢ₊₁ ÷ ωᵢ
ω₁₀(::Type{P},n)  where {P <:  AbstractCOP} = sqrt(norm2(P,n+1)/norm2(P,n))


# Can't  change  N  here
function Base.setindex!(p::AbstractCOP{T,X,N}, value::Number, idx::Int) where {T, X,N}

    ## widen size...
    idx < 0 &&  throw(ArgumentError("Negative index"))
    val = T(value)
    d = N - 1
    if idx > d || (idx == d && iszero(value))
        throw(ArgumentError("Polynomials of type `AbstractCOP` have a  fixed size parameter, N, which  can't  be changed through  assignment. Make new polynomial  instance?"))
    end
    setindex!(p.coeffs, val,  idx+1)
end


Polynomials.degree(p::AbstractCOP{T,X,N})  where {T,X,N} = N-1
Polynomials.isconstant(p::AbstractCOP) = degree(p) <=  0

##  Evaluation

"""
    Clenshaw evaluation of an orthogonal polynomial 
"""
function eval_cop(P::Type{<:AbstractCOP{T,X,N}}, cs, x::S) where {T,X,N,S}
    N == 0 && return zero(T) * zero(S)
    N == 1 && return (cs[1] * k0(P)) *  one(S)
    _eval_cop(P,cs, x)
end

function _eval_cop(P::Type{<:AbstractCOP{T,X,N}}, cs, x::S) where {T,X,N,S}
    if @generated
        quote
            Δ0 = cs[end - 1]
            Δ1 = cs[end]
            @inbounds for i in N-1:-1:2
                Δ0, Δ1 = cs[i - 1] - Δ1 * Cn(P, i-1), Δ0 + Δ1 * muladd(x, An(P,i-1), Bn(P,i-1))
            end
            p₀ = k0(P)
            p₁ =  muladd(x, An(P,0),  Bn(P,0)) * p₀
            Δ0 * p₀  + Δ1 * p₁
        end
    else
        clenshaw_eval(P, cs, x)
    end
end

#  evaluate  basis vector through hypergeometric formulation
(B::Basis{P})(x) where  {P <: AbstractCOP} = eval_basis(P, B.n, x)

# Evaluate a basis vector without realizing it and without using Clenshaw
eval_basis(::Type{P}, n, x) where {P <: AbstractCOP} = classical_hypergeometric(P, n, x)

"""
    eval_hyper(::P, cs,  x::S)

Evaluate polynomial `P(cs)` by  computing value for each basis vector from its hypergeometric representation
"""
function eval_hyper(P::Type{<:AbstractCOP}, cs, x::S) where {S}
    isempty(cs) &&  return zero(S)
    tot = cs[1] * one(S)
    for i in 2:length(cs)
        tot += cs[i] * Basis(P, i-1)(x)
    end
    tot
end
