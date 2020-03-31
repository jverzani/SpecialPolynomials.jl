
struct Newton{S<:Number, T <: Number} <: AbstractInterpolatingPolynomial{T}
    xs::Vector{S}
    coeffs::Vector{T}
    var::Symbol
    function  Newton(xs::Vector{S}, coeffs::Vector{T}, var::Symbol=:x) where {S,T}
        xs  = unique(xs)
        length(coeffs) == length(xs) || throw(ArgumentError("xs and coeffs must have the same length"))
        new{S, T}(xs, coeffs, :x)
    end
end

export Newton

basis_symbol(::Type{<:Newton}) = "p"

## Boilerplate code reproduced here, as there are two type parameters
Base.convert(::Type{P}, p::P) where {P <: Newton} =  p
Base.promote_rule(::Type{Newton{S, T}}, ::Type{R}) where {S, T, R <: Number} = Newton{S, promote_type(T,R)}

Polynomials.domain(::Type{<:Newton}) = Polynomials.Interval(-Inf, Inf)
Polynomials.variable(p::Newton{S, T}, var::Polynomials.SymbolLike=:x) where {S, T} = IMPLEMENTME()


function Base.:+(p1::P, p2::P) where {P <: Newton}
    p1.var == p2.var || throw(ArgumentError("p1 and p2 must share the same variable"))
    cs =  coeffs(p1) +  coeffs(p2)
    Newton(p1.xs, cs, p1.var)
end

## XXX
function Base.:*(p1::P, p2::P) where {P <: Newton}
    throw(ArgumentError("Multiplication not defined for Newton interpolants."))
end

function Base.:*(c::Number, p::P) where {P <: Newton}
    Newton(p.xs, c*coeffs(p), p.var)
end
Base.:*(p::P, c::Number) where {P <: Newton} = c*p

Base.:-(p::P) where {P<:Newton} = (-1)*p


function  Base.convert(Q::Type{<:Polynomial},  p::Newton{S,T})  where  {S, T}
    q = zero(Q)/1
    x =  variable(q)
    cs = coeffs(p)
    xs = p.xs
    for i in  eachindex(xs)
        ci = cs[i]
        li = prod((x - xs[j]) / (xs[i] - xs[j]) for j in eachindex(xs) if j != i)
        q += ci * li
    end
    q
end


function Polynomials.fit(P::Type{<:Newton},
                         xs::AbstractVector{S},
                         ys::AbstractVector{T};
                         var = :x,) where {S, T}
    Newton(xs, ys, var)
end
