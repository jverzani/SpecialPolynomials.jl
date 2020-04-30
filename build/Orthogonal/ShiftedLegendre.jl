"""
    ShiftedLegendre{T}

The shifted [Legendre](https://en.wikipedia.org/wiki/Legendre_polynomials#Shifted_Legendre_polynomials) polynomials have `P̃(x) = P(2x-1)` and are defined over the domain `[0,1]`.

cf. [`Legendre`](@ref)

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> p = ShiftedLegendre([1,2,3])
ShiftedLegendre(1⋅L̃_0(x) + 2⋅L̃_1(x) + 3⋅L̃_2(x))

julia> convert(Polynomial, p)
Polynomial(2//1 - 14//1*x + 18//1*x^2)

julia> q = Legendre([1,2,3])
Legendre(1⋅L_0(x) + 2⋅L_1(x) + 3⋅L_2(x))

julia> p(1/4) == q(2*(1/4) - 1)
true
```
"""
struct ShiftedLegendre{T <: Number} <: AbstractLegendre{T}
    coeffs::Vector{T}
    var::Symbol
    function ShiftedLegendre{T}(coeffs::AbstractVector{T}, var::Symbol) where {T <: Number}
        length(coeffs) == 0 && return new{T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{T}(coeffs[1:last], var)
    end
end

export ShiftedLegendre

Polynomials.@register ShiftedLegendre

basis_symbol(::Type{<:ShiftedLegendre}) = "L̃"

Polynomials.domain(::Type{<:ShiftedLegendre}) = Polynomials.Interval(0, 1)
Polynomials.variable(::Type{P}, var::Polynomials.SymbolLike=:x) where {P <: ShiftedLegendre} = P([0, 1], var)

weight_function(::Type{<:ShiftedLegendre}) = x -> one(x)
generating_function(::Type{<:ShiftedLegendre}) = (t, x)  -> 1/sqrt(1 - 2*(2x-1)*t +t^2)

# 3 point recursion
An(::Type{ShiftedLegendre{T}}, n) where {T <: Integer} = (4n + 2)//(n + 1)
An(::Type{<:ShiftedLegendre}, n) = (4n + 2)/(n + 1)
Bn(::Type{ShiftedLegendre{T}}, n) where {T <: Integer} = -(2n + 1)//(n + 1)
Bn(::Type{<:ShiftedLegendre}, n) = -(2n + 1)/(n + 1)
Cn(::Type{ShiftedLegendre{T}}, n) where {T <: Integer} = -n //(n  + 1)
Cn(::Type{<:ShiftedLegendre}, n) = -n / (n  + 1)


(ch::ShiftedLegendre{T})(x::S) where {T,S} = orthogonal_polyval(ch, x)

function Base.convert(::Type{P}, p::Legendre{T}) where {P <: ShiftedLegendre, T}
    x = variable(Legendre{T})
    P(coeffs(p(2x-1)), p.var)
end

function Base.convert(::Type{P}, p::ShiftedLegendre{T}) where {P <: Legendre, T}
    q = convert(Polynomial{T}, p)
    convert(P, q)
end
