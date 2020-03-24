struct ShiftedLegendre{T <: Number} <: OrthogonalPolynomial{T}
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
weight_function(::Type{ShiftedLegendre{T}}) where {T} = x -> one(x)
Polynomials.variable(::Type{P}, var::Polynomials.SymbolLike=:x) where {P <: ShiftedLegendre} = P([0, 1], var)

# Bonnet's expresssion
An(::Type{<:ShiftedLegendre{T}}, n) where {T <: Integer} = (4n + 2)//(n + 1)
An(::Type{<:ShiftedLegendre}, n) = (4n + 2)/(n + 1)
Bn(::Type{<:ShiftedLegendre{T}}, n) where {T <: Integer} = -(2n + 1)//(n + 1)
Bn(::Type{<:ShiftedLegendre}, n) = -(2n + 1)/(n + 1)
Cn(::Type{<:ShiftedLegendre}, n) = -n //(n  + 1)


(ch::ShiftedLegendre{T})(x::S) where {T,S} = orthogonal_polyval(ch, x)

function Base.convert(::Type{P}, p::Legendre{T}) where {P <: ShiftedLegendre, T}
    x = variable(Legendre{T})
    P(coeffs(p(2x-1)), p.var)
end

function Base.convert(::Type{P}, p::ShiftedLegendre{T}) where {P <: Legendre, T}
    q = convert(Polynomial, p)
    convert(P, q)
end
