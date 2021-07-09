"""
    AbstractInterpolatingPolynomial{T,X}

Abstract type for interpolating polynomials.

These are polynomial representations of `p(x)` satisfying `p(x_i) =
y_i` for a specified set of `x` values and `y` values.

For a collection of points `(x_0,y_0), ..., (x_n, y_n)` there is a
*unique* polynomial of degree `n` or less satisfying
`p(x_i)=y_i`. This fact allows the specification of `p(x)` using a
vector of coefficients relative to some set of basis vectors.

The two main types, `Lagrange` and `Newton`, store the nodes within
the instance. In particular, the type does not contain all the
information needed to describe the family. So methods like
`convert(::Type, p)` will not work. Use `fit(Type, xs, p)`, as
appropriate, instead.

"""
abstract type AbstractInterpolatingPolynomial{T,X} <: AbstractSpecialPolynomial{T,X} end


function Polynomials.showterm(io::IO, ::Type{P}, pj::T, var, j, first::Bool, mimetype) where {N, T, P <: AbstractInterpolatingPolynomial}
    iszero(pj) && return false
    !first &&  print(io, " ")
    print(io, Polynomials.hasneg(T)  && Polynomials.isneg(pj) ? "- " :  (!first ? "+ " : ""))
    print(io, "$(abs(pj))⋅$(basis_symbol(P))_$(j)($var)")
    return true
end

basis_symbol(::Type{<:AbstractInterpolatingPolynomial}) = "P"

function Polynomials.variable(::Type{P}, var::Polynomials.SymbolLike=:x) where {P <: AbstractInterpolatingPolynomial}
    @warn "No `variable` function is defined for a type, just an instance"
    throw(MethodError())
end


# pick new nodes for product
function _new_nodes(xs::Vector{S}, ys::Vector{T}) where {S, T}
    sort!(ys)
    n = length(ys)
    out = setdiff(ys, xs)
    length(out) >= n-1 && return out[1:n-1]
    out  = out/1
    tmp = intersect(xs, ys)
    sort!(tmp)
    for i in 2:length(tmp)
        a,b = tmp[i-1], tmp[i]
        c = (a+b)/2
        while c in xs
            c = (a+c)/2
        end
        push!(out, c)
    end
    out
end
