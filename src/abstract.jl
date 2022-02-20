"""
    AbstractSpecialPolynomial{T,X}

An abstract type to distinguish the different polynomial types in this package.

The concrete types specify different bases for the space of polynomials of degree `n` or less.

This package includes:

* several classic orthogonal polynomials.
* Newton and Lagrange interpolating polynomials
* Bernstein polynomials

As many of the methods for the base `Polynomials` class are directly coded if possible, but quite a few
depend on conversion to the base `Polynomial` type (which uses the standard polynomial basis).

"""
abstract type AbstractSpecialPolynomial{T,X} <: Polynomials.AbstractPolynomial{T,X} end

# polynomial like Vector{T} with variable
Base.eltype(::Type{<:AbstractSpecialPolynomial{T}}) where {T} = T
Base.eltype(::Type{<:AbstractSpecialPolynomial}) = Float64

# use ArgumentError to give message. Also could use hint.
Base.convert(::Type{Q}, p::P) where {P<:AbstractPolynomial,Q<:AbstractSpecialPolynomial} =
    throw(
        ArgumentError(
            "There is no `convert` method defined for a polynomial of type $P to one of type $Q. Maybe try converting through the `Polynomial` type (e.g `convert(Q, convert(Polynomial, p))`)",
        ),
    )

## Display

function unicode_subscript(io, j)
    a = ("⁻", "", "", "₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉")
    for i in string(j)
        print(io, a[Int(i) - 44])
    end
end

function Polynomials.showterm(
    io::IO,
    ::Type{P},
    pj::T,
    var,
    j,
    first::Bool,
    mimetype,
) where {N,T,P<:AbstractSpecialPolynomial}
    iszero(pj) && return false
    !first && print(io, " ")

    if Polynomials.hasneg(T) && Polynomials.isneg(pj)
        print(io, "- ")
        pj = -pj
    else
        !first && print(io, "+ ")
    end

    Polynomials.printcoefficient(io, pj, 1, mimetype) # j != 0; want to have ()
    print(io, "⋅")
    print(io, basis_symbol(P))

    unicode_subscript(io, j)
    print(io, "($(var))")
    return true
end

# Conversion is *possible* when multiplication  in `P` is defined
# but otherwise,  it must be handled specially
# function Base.convert(::Type{P}, q::Q) where {P <: AbstractSpecialPolynomial, Q <: AbstractSpecialPolynomial}
#     convert(P, convert(Polynomial, q))
# end

# return the domain as a tuple, not an interval object
function Base.extrema(::Type{P}) where {P<:AbstractSpecialPolynomial}
    dom = domain(P)
    first(dom), last(dom)
end
Base.extrema(p::P) where {P<:AbstractSpecialPolynomial} = extrema(P)

## Polynomial operations

##
## --------------------------------------------------
##

function Polynomials.derivative(p::P, order::Integer=1) where {P<:AbstractSpecialPolynomial}
    T = eltype(one(P))
    q = convert(Polynomial{T}, p)
    convert(⟒(P), derivative(q, order))
end

function Polynomials.integrate(p::P) where {P<:AbstractSpecialPolynomial}
    q = convert(Polynomial, p)
    integrate(q)
end

function Base.divrem(num::P, den::P) where {P<:AbstractSpecialPolynomial}
    p1 = convert(Polynomial, num)
    p2 = convert(Polynomial, den)
    q, r = divrem(p1, p2)
    convert.(P, (q, r))
end

function Polynomials.companion(p::P) where {P<:AbstractSpecialPolynomial}
    companion(convert(Polynomial, p))
end

function Polynomials.vander(
    ::Type{P},
    x::AbstractVector{T},
    n::Integer,
) where {P<:AbstractSpecialPolynomial,T}
    N = length(x) - 1
    R = typeof(one(T) / one(T))
    V = zeros(R, N + 1, n + 1)

    for j in 0:n
        bj = Polynomials.basis(P, j)
        for i in 0:N
            V[i + 1, j + 1] = bj(x[i + 1])
        end
    end

    V
end
