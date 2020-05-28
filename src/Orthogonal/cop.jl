## common functionality for  COP classical orthogonal polynomials
## These are described by 5  values: a,.b.c.d.e

abstract type AbstractCOP{T,N} <: AbstractOrthogonalPolynomial{T} end
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

# kn is the leading term (Section 3 table)
leading_term(P::Type{<:AbstractCOP},  n::Int) =  kn(P, n)

# Set defaults to be monic
kn(::Type{P},  n::Int) where{P <: AbstractCOP} = one(eltype(one(P))) #  need one(P), as eltype(Polynomial{T}) != T

# k₍ᵢ₊₁₎/kᵢ
k1k0(::Type{P},  i::Int) where {P<:AbstractCOP} =  kn(P,i+1)/kn(P,i)
# k₍ᵢ₊₁₎/k₍ᵢ₋₁₎
k1k_1(::Type{P},  i::Int) where {P<:AbstractCOP} =  kn(P,i+1)/kn(P,i-1)

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

