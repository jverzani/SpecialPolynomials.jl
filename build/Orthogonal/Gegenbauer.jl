"""
   Gegenbauer{α, T <: Number}

The Gegenbauer family of orthogonal polynomials has weight function `(1-x^2)^(α-1/2)`. The parameter `α` is  specified  in the constructor.

```jldoctest
julia> p =  Gegenbauer{1/2}([1,2,3])
Gegenbauer(1⋅C^(0.5)_0(x) + 2⋅C^(0.5)_1(x) + 3⋅C^(0.5)_2(x))

julia> convert(Polynomial, p)
Polynomial(-0.5 + 2.0*x + 4.5*x^2)
```

"""
struct Gegenbauer{α, T <: Number} <: OrthogonalPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function Gegenbauer{α, T}(coeffs::AbstractVector{T}, var::Symbol=:x) where {α, T <: Number}
        length(coeffs) == 0 && return new{α, T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{α, T}(coeffs[1:last], var)
    end
    function Gegenbauer{α}(coeffs::AbstractVector{T}, var::Symbol=:x) where {α, T <: Number}
        Gegenbauer{α, T}(coeffs, var)
    end
end

export Gegenbauer

# boilerplate code needed, as two type parameters
Base.convert(::Type{P}, p::P) where {P <: Gegenbauer} =  p
Base.promote_rule(::Type{Gegenbauer{α, T}}, ::Type{Gegenbauer{α, S}}) where {α, T, S} = Gegenbauer{α, promote_type(T,S)}
Base.promote_rule(::Type{Gegenbauer{α, T}}, ::Type{S}) where {α, T, S <: Number} = Gegenbauer{α, promote_type(T,S)}
Gegenbauer{α}(n::Number, var=:x) where  {α}= Gegenbauer{α}([n], var)
Gegenbauer{α, T}(n::S, var=:x) where {α, T, S <: Number} = Gegenbauer{α, T}(T[n], var)
function Gegenbauer{α, T}(xs::AbstractVector{S}, var=:x) where {α, T, S <: Number}
    R = promote_type(T, S)
    Gegenbauer{α, R}(R.(xs), var)
end


Polynomials.domain(::Type{<:Gegenbauer}) = Polynomials.Interval(-1, 1)
basis_symbol(::Type{Gegenbauer{α, T}}) where {α, T} = "C^($α)"

weight_function(::Type{Gegenbauer{α, T}}) where {α, T} = x -> (1-x^2)^(α-1/2)
generating_function(::Type{Gegenbauer{α, T}}) where {α, T} = (t,x) -> begin
    1/(1-2x*t +t^2)^α
end


An(::Type{Gegenbauer{α,T}}, n) where {α,T} = 2(n+α)/(n+1)
An(::Type{Gegenbauer{α}}, n) where {α} = 2(n+α)/(n+1)
Bn(::Type{Gegenbauer{α,T}}, n) where {α,T} = zero(T)/1
Bn(::Type{Gegenbauer{α}}, n) where {α} = 0/1
Cn(::Type{Gegenbauer{α,T}}, n) where {α,T} = - (n - 1 + 2α)/(n+1)
Cn(::Type{Gegenbauer{α}}, n) where {α} = - (n - 1 + 2α)/(n+1)

# compute <Jn, Jn>
function norm2(::Type{Gegenbauer{α, T}}, n) where{α, T}
    pi * 2^(1-2α) * gamma(n + 2α) / (gamma(n+1) * (n+α) * gamma(α)^2)

    gamma(n + α + 1) / gamma(n+1)
end

(ch::Gegenbauer)(x::S) where {T,S} = orthogonal_polyval(ch, x)



# XXX could tidy up by doing higher orders at once
# XXX this returns answer in a different basis.
function Polynomials.derivative(p::Gegenbauer{α, T}, order::Int=1) where {α,T}
    order == 0 && return p
    order < 0 && throw(ArgumentError("order must be ≥ 0"))

    xs = coeffs(p)
    qs = 2 * α * xs[2:end]
    q = Gegenbauer{α+1}(qs, p.var)

    if order > 1
        return derivative(q, order-1)
    else
        return q
    end

end


# integral: https://mathworld.wolfram.com/GegenbauerPolynomial.html
# \int Pn = 1/(2(n+α)) * [P_{n+1} - P_{n-1}]
function Polynomials.integrate(p::Gegenbauer{α, T}, C::Number=0) where {α,T}

    R = eltype(one(T)/1)
    d = degree(p)
    qs = zeros(R, d+2)
    qs[1] = -p[1]/(2*(1+α))

    for i in 1:(d-1)
        qs[i+1] = p[i-1]/(2(i-1+α)) - p[i+1]/(2(i+1+α))
    end
    qs[d+1] = p[d-1] / (2(d-1+α))
    qs[d+2] = p[d] / (2(d+α))

    q = Gegenbauer{α}(qs, p.var)
    q = q - q(0) + R(C)

    return q
end
