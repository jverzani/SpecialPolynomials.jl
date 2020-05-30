## common functionality for  COP classical orthogonal polynomials
## These are described by 5  values: a,.b.c.d.e

abstract type AbstractCOP{T,N} <: AbstractOrthogonalPolynomial{T} end

function Base.promote_rule(P::Type{<:Polynomials.AbstractPolynomial{T}},
                           Q::Type{<:AbstractCOP{S}}) where {T,S}
    ϛ(Q){promote_type(T, S)}
end

Base.promote_rule(P::Type{<:AbstractCOP{S}},
                  Q::Type{<:Polynomials.AbstractPolynomial{T}}) where {T,S} =
                      ϛ(P){promote_type(T, S)}

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
Polynomials.domain(::Type{<:AbstractCOP}) = Polynomials.Interval(-Inf,  Inf)
classical_hypergeometric(P::Type{<:AbstractCOP}, n, x) = throw(ArgumentError("No default method"))

"""
   abcde

A named tuple returning  the  constants a,b,c,d,e  for a CCOP type with
(a⋅x²+b⋅x+c)*P₍ᵢ₊₂₎'' + (d⋅x + e) * P₍ᵢ₊₁₎ + λᵢ Pᵢ = 0.
"""
abcde(::Type{<:AbstractCOP}) = throw(ArgumentError("No default method"))


"""
    k1k0

Let  `kᵢ` be the leading  coeffiecient of the  polynomial  in  the standard basis. This  function implement  `k₍ᵢ₊₁₎/kᵢ` for `i ≥  0`.

With  the assumption that  `k₀ = 1`, the values `kᵢ` and `k₍ᵢ₊₁₎/k₍ᵢ₋₁₎` can  be generated.

The default value  leads to monic polynomials.
"""
k1k0(::Type{P},  i::Int) where {P<:AbstractCOP} =  one(eltype(P))

# Set defaults to be monic
# For non-monic subtypes of `AbstractCOP` only need k1k0  to be defined, as `kn` and `k1k_1` come for free

# leading term (Section 3 table)  is kn
leading_term(P::Type{<:AbstractCOP},  n::Int) =  kn(P, n)

# kn = prod(k1k0(i) for i in 0:n-1)  *  kn(P,0)
# **Assume** basis(P,0) = 1, so kn(P,0)  = 1
kn(::Type{P},  n::Int) where{P <: AbstractCOP} = iszero(n) ?  one(eltype(P)) : prod(k1k0(P,i) for i in  0:n-1)

# k₍ᵢ₊₁₎/k₍ᵢ₋₁₎ =  (k₍ᵢ₊₁₎/kᵢ) ⋅ (kᵢ/k₍ᵢ₋₁₎)
function k1k_1(::Type{P},  i::Int) where {P<:AbstractCOP}
    @assert i > 0
    k1k0(P,i)*k1k0(P,i-1)
end


# Can't  change  N  here
function Base.setindex!(p::AbstractCOP{T,N}, value::Number, idx::Int) where {T, N}

    ## widen size...
    idx < 0 &&  throw(ArgumentError("Negative index"))
    val = T(value)
    d = N - 1
    if idx > d || (idx == d && iszero(value))
        throw(ArgumentError("Polynomials of type `AbstractCOP` have a  fixed size parameter, N, which  can't  be changed through  assignment. Make new polynomial  instance?"))
    end
    setindex!(p.coeffs, val,  idx+1)
end


Polynomials.degree(p::AbstractCOP{T,N})  where {T,N} = N-1
Polynomials.isconstant(p::AbstractCOP) = degree(p) <=  0


"""
    Clenshaw evaluation of an orthogonal polynomial 
"""
function eval_cop(P::Type{<:AbstractCOP{T,N}}, cs, x::S) where {T,N,S}
    if @generated
        N == 0 && return zero(T) * zero(S)
        N == 1 && return cs[1] * one(S)
        #SS = eltype(one(S))
        Δ0 = :(cs[N-1])
        Δ1 = :(cs[N])
        for i in N-1:-1:2
            a = :(cs[i - 1] - Δ1 * Cn(P, i-1))
            b = :(Δ0 + Δ1 * muladd(x, An(P,i-1),Bn(P,i-1)))
            Δ0 = :(a)
            Δ1 = :(b)
        end
        Δ0 + Δ1* muladd(x, An(P,0), Bn(P,0))
    else
        clenshaw_eval(P, cs, x)
    end
end

