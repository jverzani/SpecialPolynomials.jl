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


# some shortcuts
basis =  Polynomials.basis
export basis

@generated function constructorof(::Type{T}) where T
    getfield(parentmodule(T), nameof(T))
end
⟒(P::Type{<:AbstractSpecialPolynomial}) = constructorof(P) # to strip off type parameters. See Polynomials.constructorof
# set up some defaults

function Polynomials.showterm(io::IO, ::Type{P}, pj::T, var, j, first::Bool, mimetype) where {P <: AbstractSpecialPolynomial, T}
    iszero(pj) && return false
    !first &&  print(io, " ")
    print(io, Polynomials.hasneg(T) && pj < 0 ? "- " :  (!first ? "+ " : ""))
    print(io, "$(abs(pj))⋅e_$j($var)")
    return true
end

# Conversion is *possible* when multiplication  in `P` is defined
# but otherwise,  it must be handled specially
# function Base.convert(::Type{P}, q::Q) where {P <: AbstractSpecialPolynomial, Q <: AbstractSpecialPolynomial}
#     convert(P, convert(Polynomial, q))
# end

function Base.:*(p1::P, p2::Q) where {P <: AbstractSpecialPolynomial, Q <: AbstractSpecialPolynomial}
    R = promote_type(P, Q)
    convert(⟒(R), convert(Polynomial, p1) * convert(Polynomial, p2))
end


function Base.:+(p1::P, p2::P) where {P <: AbstractSpecialPolynomial}
    p1.var != p2.var && throw(ArgumentError("Polynomials must have same variable"))
    n = max(length(p1), length(p2))
    c = [p1[i] + p2[i] for i = 0:n]
    return P(c, p1.var)
end



function Polynomials.derivative(p::P, order::Integer = 1) where {P <: AbstractSpecialPolynomial}
    q = convert(Polynomial, p)
    convert(⟒(P), derivative(q, order))
end


function Polynomials.integrate(p::P, C::Number=0) where {P <: AbstractSpecialPolynomial}
    q = convert(Polynomial, p)
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
