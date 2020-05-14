"""
    AbstractSpecialPolynomial{T}

An abstract type to distinguish the different polynomial types in this package.

The concrete types specify different bases for the space of polynomials of degree `n` or less. 

This package includes: 

* several classic orthogonal polynomials.
* Newton and Lagrange interpolating polynomials
* Bernstein polynomials

As many of the methods for the base `Polynomials` class are directly coded if possible, but quite a few
depend on conversion to the base `Polynomial` type (which uses the standard polynomial basis).

"""
abstract type AbstractSpecialPolynomial{T} <: Polynomials.AbstractPolynomial{T} end

# polynomial like Vector{T} with variable
Base.eltype(::Type{<:AbstractSpecialPolynomial{T}}) where {T} = T
Base.eltype(::Type{<:AbstractSpecialPolynomial}) = Float64

# to strip off type parameters. See Polynomials.constructorof
# We want ⟒(P{α,T,N}) = P{α} where {T, N}x
@generated function constructorof(::Type{T}) where T
      getfield(parentmodule(T), nameof(T))
 end
⟒(P::Type{<:AbstractSpecialPolynomial}) = constructorof(P) 
⟒(::Type{<:Polynomial}) = Polynomial

## Display

function unicode_subscript(io, j)
    a = ("⁻","","","₀","₁","₂","₃","₄","₅","₆","₇","₈","₉")
    for i in string(j)
        print(io, a[Int(i)-44])
    end
end

function Polynomials.showterm(io::IO, ::Type{P}, pj::T, var, j, first::Bool, mimetype) where {N, T, P <: AbstractSpecialPolynomial}
    iszero(pj) && return false
    !first &&  print(io, " ")
    print(io, Polynomials.hasneg(T)  && Polynomials.isneg(pj) ? "- " :  (!first ? "+ " : ""))
    print(io, "$(abs(pj))⋅$(basis_symbol(P))")
    unicode_subscript(io, j)
    print(io,"($(var))")
    return true
end

# function Polynomials.showterm(io::IO, ::Type{P}, pj::T, var, j, first::Bool, mimetype) where {P <: AbstractSpecialPolynomial, T}
#     iszero(pj) && return false
#     !first &&  print(io, " ")
#     print(io, Polynomials.hasneg(T) && pj < 0 ? "- " :  (!first ? "+ " : ""))
#     print(io, "$(abs(pj))⋅e_$j($var)")
#     return true
# end

# Conversion is *possible* when multiplication  in `P` is defined
# but otherwise,  it must be handled specially
# function Base.convert(::Type{P}, q::Q) where {P <: AbstractSpecialPolynomial, Q <: AbstractSpecialPolynomial}
#     convert(P, convert(Polynomial, q))
# end

## Polynomial operations

# polynomial operations -p, p+c, p*c, p/c, p+q, p*q

# faster than +(promote(p,c)...); as that allocates twice
# useful to define p + P(:θ) (mismatched symbols
function Base.:+(p::P, c::S) where {P <: AbstractSpecialPolynomial, S<:Number}
    R = promote_type(eltype(p), S)
    as = R[a for a in coeffs(p)]
    if length(as) >= 1
        as[1] += c  # *assume* basis(P,0) = 1
    else
        push!(as, c)
    end
    ⟒(P)(as, p.var)
end
   
function Base.:*(p::P, c::S) where {P<:AbstractSpecialPolynomial, S<:Number}
    R = promote_type(eltype(p), S)
    as = R[a for a in coeffs(p)]
     ⟒(P)(as * c, p.var)
end

#function Base.:/(p::P, c::S) where {P <: AbstractSpecialPolynomial, S<:Number}
#    as = copy(coeffs(p))
#    ⟒(P)(as / c, p.var)
#end


# function Base.:+(p1::P, p2::P) where {P <: AbstractSpecialPolynomial}
#     p1.var != p2.var && throw(ArgumentError("Polynomials must have same variable"))
#     n = max(length(p1), length(p2))
#     c = [p1[i] + p2[i] for i = 0:n]
#     return P(c, p1.var)
# end

# function Base.:+(p::P, q::Q) where {P <: AbstractSpecialPolynomial, Q <: AbstractSpecialPolynomial}

#     isconstant(p) && return q + p[0]
#     isconstant(q) && return p + q[0]
#     p.var != q.var && throw(ArgumentError("Variables don't  match"))
    
#     if ⟒(P) == ⟒(Q)
#         d = max(degree(p), degree(q))
#         as = [p[i]+q[i] for i in 0:d]
#         return ⟒(P)(as, p.var)
#     else
#         p1, q1 = promote(p, q)
#         p1 + q1
#     end
# end

#### XXX NOT needed,  should be  handled  by  promote_type
#function Base.:*(p1::P, p2::Q) where {P <: AbstractSpecialPolynomial, Q <: AbstractSpecialPolynomial}
#    R = promote_type(P,Q)
#    convert(⟒(R), convert(Polynomial, p1) * convert(Polynomial, p2))
#end



##
## --------------------------------------------------
##

function Polynomials.derivative(p::P, order::Integer = 1) where {P <: AbstractSpecialPolynomial}
    T = eltype(one(P))
    q = convert(Polynomial{T}, p)
    convert(⟒(P), derivative(q, order))
end


function Polynomials.integrate(p::P, C::Number=0) where {P <: AbstractSpecialPolynomial}
    T = eltype(one(p))
    q = convert(Polynomial{T}, p)
    convert(⟒(P), integrate(q, C))
end

function Polynomials.integrate(p::P, a::Number, b::Number) where {P <: AbstractSpecialPolynomial}
    q = convert(Polynomial, p)
    integrate(q, a, b)
end


function Base.divrem(num::P, den::P) where {P <: AbstractSpecialPolynomial}

    p1 = convert(Polynomial, num)
    p2 = convert(Polynomial, den)
    q,r = divrem(p1, p2)
    convert.(P, (q,r))

end

function Polynomials.companion(p::P) where {P <: AbstractSpecialPolynomial}
    companion(convert(Polynomial, p))
end

function Polynomials.vander(::Type{P}, x::AbstractVector{T}, n::Integer) where {P <: AbstractSpecialPolynomial, T}
    N = length(x) - 1
    R = typeof(one(T)/one(T))
    V = zeros(R, N+1, n+1)

    for j in 0:n
        bj = Polynomials.basis(P, j)
        for i in 0:N
            V[i+1, j+1] = bj(x[i+1])
        end
    end

    V
end
