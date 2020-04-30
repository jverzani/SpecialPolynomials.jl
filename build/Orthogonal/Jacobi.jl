"""
    Jacobi{α,  β, T}

Implements the [Jacobi](https://en.wikipedia.org/wiki/Jacobi_polynomials) polynomials. These have weight function `w(x) = (1-x)^α ⋅ (1+x)^β` over the domain `[-1,1]`. Many orthogonal families are special cases. The parameters are  specified to the constructors:

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> p = Jacobi{-1/2, -1/2}([0,0,1])
Jacobi(1⋅J^(α, β)_2(x))

julia> convert(Polynomial, p)
Polynomial(-0.375 + 0.75*x^2)

julia> monic(p::Polynomial) = p/p[end];

julia> (monic ∘ convert).(Polynomial, (p, basis(Chebyshev, 2))) # scaled version of each other
ERROR: UndefVarError: ChebyshevTT not defined
Stacktrace:
 [1] top-level scope at none:1
```

"""
struct Jacobi{α,  β, T <: Number} <: OrthogonalPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function Jacobi{α, β, T}(coeffs::AbstractVector{T}, var::Symbol=:x) where {α, β, T <: Number}
        α <= -1 && throw(ArgumentError("α > -1 is necessary"))
        β <= -1 && throw(ArgumentError("β > -1 is necessary"))
        length(coeffs) == 0 && return new{α, β, T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{α, β, T}(coeffs[1:last], var)
    end
end

export Jacobi

Polynomials.@register2 Jacobi

basis_symbol(::Type{<: Jacobi{α, β}}) where {α, β} = "J^(α, β)"

Polynomials.domain(::Type{<:Jacobi{α, β}}) where {α, β} = Polynomials.Interval(-1, 1, β < 0, α < 0)
weight_function(::Type{<:Jacobi{α, β}}) where {α, β} = x -> (1-x)^α *  (1+x)^β
generating_function(::Type{<:Jacobi{α, β}}) where {α, β} = (t,x) -> begin
    R = sqrt(1 - 2x*t+t^2)
    2^(α + β) * 1/R  * (1 - t + R)^(-α) * (1 + t + R)^(-β)
end


An(::Type{<:Jacobi{α, β}}, n) where {α, β} = iszero(n) ? (α+β+2)/2  : (2n+α+β+1) * (2n+α+β+2) / (2(n+1)*(n+α+β+1))
Bn(::Type{<:Jacobi{α, β}}, n) where {α, β} = iszero(n) ? (α-β)/2    : (α^2-β^2)  * (2n+α+β+1) / (2(n+1)*(n+α+β+1)*(2n+α+β))
Cn(::Type{<:Jacobi{α, β}}, n) where {α, β} = iszero(n) ? α*β*(α+β+2)/((α+β)*(α+β+1)) : -(n+α)*(n+β)*(2n+α+β+2)/((n+1)*(n+α+β+1)*(2n+α+β))

(ch::Jacobi)(x::S) where {T,S} = orthogonal_polyval(ch, x)

# compute <Jn, Jn>
function norm2(::Type{<:Jacobi{α, β}}, n) where{α, β}
    α > -1 && β > -1 || throw(ArgumentError("α, β > -1 is necessary"))
    2^(α+β+1)/(2n+α+β+1) * (gamma(n+α+1) *  gamma(n + β +1))/(gamma(n+α+β+1)*gamma(n+1))
end
