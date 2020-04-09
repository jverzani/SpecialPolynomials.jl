"""
    AbstractSpecialPolynomial{T}

An abstract type to distinguish the families of polynomials in this package.

The families consist of different bases for the space of polynomials of degree `n` or less. 

This package includes: 

* several orthgonal polynomial familes.
* Newton and Lagrange interpolating polynomials
* Bernstein polynomials

As many of the methods for the base `Polynomials` class are directly coded if possible, but quite a few
depend on conversion to the base `Polynomial` type (which uses the standard polynomial basis).

"""
abstract type AbstractSpecialPolynomial{T} <: Polynomials.AbstractPolynomial{T} end

# Macros to register POLY{α, T} and POLY{α, β, T}
# modified from `Polynomials.jl`
macro register1(name)
    poly = esc(name)
    quote
        Base.convert(::Type{P}, p::P) where {P<:$poly} = p
        Base.convert(P::Type{<:$poly}, p::$poly) where {T} = P(coeffs(p), p.var)
        Base.promote_rule(::Type{$poly{α,T}}, ::Type{$poly{α,S}}) where {α,T,S} =
            $poly{α,promote_type(T, S)}
        Base.promote_rule(::Type{$poly{α,T}}, ::Type{S}) where {α,T,S<:Number} = 
            $poly{α,promote_type(T, S)}
        $poly{α}(coeffs::AbstractVector{T}, var::Polynomials.SymbolLike = :x) where {α,T} =
            $poly{α,T}(coeffs, Symbol(var))
        $poly{α,T}(x::AbstractVector{S}, var = :x) where {α,T,S<:Number} = $poly{α}(T.(x), var)
        $poly{α,T}(n::Number, var = :x) where {α,T} = $poly{α}(T[n], var)        
        $poly{α}(n::Number, var = :x) where {α} = $poly{α}([n], var)
        $poly{α,T}(var=:x) where {α, T} = variable($poly{α,T}, var)
        $poly{α}(var=:x) where {α} = variable($poly{α}, var)
    end
end

# Macro to register POLY{α, β, T}
macro register2(name)
    poly = esc(name)
    quote
        Base.convert(::Type{P}, p::P) where {P<:$poly} = p
        Base.convert(P::Type{<:$poly}, p::$poly) where {T} = P(coeffs(p), p.var)
        Base.promote_rule(::Type{$poly{α,β,T}}, ::Type{$poly{α,β,S}}) where {α,β,T,S} =
            $poly{α,β,promote_type(T, S)}
        Base.promote_rule(::Type{$poly{α,β,T}}, ::Type{S}) where {α,β,T,S<:Number} =
            $poly{α,β,promote_type(T, S)}
        $poly{α,β}(coeffs::AbstractVector{T}, var::Polynomials.SymbolLike = :x) where {α,β,T} =
            $poly{α,β,T}(coeffs, Symbol(var))
        $poly{α,β,T}(x::AbstractVector{S}, var = :x) where {α,β,T,S<:Number} = $poly{α,β}(T.(x), var)
        $poly{α,β,T}(n::Number, var = :x) where {α,β,T} = $poly{α}(T(n), var)
        $poly{α,β}(n::Number, var = :x) where {α,β} = $poly{α,β}([n], var)
        $poly{α,β,T}(var=:x) where {α,β, T} = variable($poly{α,β,T}, var)
        $poly{α,β}(var=:x) where {α,β} = variable($poly{α,β}, var)
    end
end


# set up some defaults

function Polynomials.showterm(io::IO, ::Type{P}, pj::T, var, j, first::Bool, mimetype) where {P <: AbstractSpecialPolynomial, T}
    iszero(pj) && return false
    !first &&  print(io, " ")
    print(io, Polynomials.hasneg(T) && pj < 0 ? "- " :  (!first ? "+ " : ""))
    print(io, "$(abs(pj))⋅e_$j($var)")
    return true
end

## One of these next two *must* be defined for the polynomial type
## otherwise they will loop (e.g. x^2 needs convert, which needs x^2...)
function Base.convert(C::Type{P}, p::Polynomial) where {P <: AbstractSpecialPolynomial}
    x = variable(C)
    p(x)
end

function Base.:*(p1::P, p2::Q) where {P <: AbstractSpecialPolynomial, Q <: AbstractSpecialPolynomial}
    R = promote_type(P, Q)
    convert(R, convert(Polynomial, p1) * convert(Polynomial, p2))
end


function Base.:+(p1::P, p2::P) where {P <: AbstractSpecialPolynomial}
    p1.var != p2.var && throw(ArgumentError("Polynomials must have same variable"))
    n = max(length(p1), length(p2))
    c = [p1[i] + p2[i] for i = 0:n]
    return P(c, p1.var)
end



function Polynomials.derivative(p::P, order::Integer = 1) where {P <: AbstractSpecialPolynomial}
    q = convert(Polynomial, p)
    convert(P, derivative(q, order))
end


function Polynomials.integrate(p::P, C::Number=0) where {P <: AbstractSpecialPolynomial}
    q = convert(Polynomial, p)
    convert(P, integrate(q, C))
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
