struct Jacobi{α,  β, T <: Number} <: OrthogonalPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function Jacobi{α, β, T}(coeffs::AbstractVector{T}, var::Symbol=:x) where {α, β, T <: Number}
        length(coeffs) == 0 && return new{α, β, T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{α, β, T}(coeffs[1:last], var)
    end
    function Jacobi{α, β}(coeffs::AbstractVector{T}, var::Symbol=:x) where {α, β, T <: Number}
        Jacobi{α, β, T}(coeffs, var)
    end
end

export Jacobi

# boilerplate code, as we have 3 type parameters
Base.convert(::Type{P}, p::P) where {P <: Jacobi} =  p
Base.convert(P::Type{<:Jacobi}, p::Jacobi) =  P(p.coeffs, p.var)
Base.promote_rule(::Type{Jacobi{α, β, T}}, ::Type{Jacobi{α, β, S}}) where {α, β, T, S} = Jacobi{α, β, promote_type(T,S)}
Base.promote_rule(::Type{Jacobi{α, β, T}}, ::Type{S}) where {α, β, T, S <: Number} = Jacobi{α, β, promote_type(T,S)}
Jacobi{α, β}(n::Number, var=:x) where  {α, β}= Jacobi{α, β}([n], var)
Jacobi{α, β, T}(n::S, var=:x) where {α, β, T, S <: Number} = Jacobi{α, β, T}(T[n], var)
Jacobi{α, β, T}(xs::AbstractVector{S}, var=:x) where {α, β, T, S <: Number} = Jacobi{α, β, T}(T.(xs), var)


basis_symbol(::Type{Jacobi{α, β, T}}) where {α, β, T} = "J^(α, β)"

Polynomials.domain(::Type{<:Jacobi}) = Polynomials.Interval(-1, 1)


weight_function(::Type{Jacobi{α, β, T}}) where {α, β, T} = x -> (1-x)^α *  (1+x)^β
generating_function(::Type{Jacobi{α, β, T}}) where {α, β, T} = (t,x) -> begin
    R = sqrt(1 - 2x*t+t^2)
    2^(α + β) * 1/R  * (1 - t + R)^(-α) * (1 + t + R)^(-β)
end


An(::Type{Jacobi{α, β, T}}, n) where {α, β, T} = iszero(n) ? (α+β+2)/2  : (2n+α+β+1) * (2n+α+β+2) / (2(n+1)*(n+α+β+1))
Bn(::Type{Jacobi{α, β, T}}, n) where {α, β, T} = iszero(n) ? (α-β)/2    : (α^2-β^2)  * (2n+α+β+1) / (2(n+1)*(n+α+β+1)*(2n+α+β))
Cn(::Type{Jacobi{α, β, T}}, n) where {α, β, T} = iszero(n) ? α*β*(α+β+2)/((α+β)*(α+β+1)) : -(n+α)*(n+β)*(2n+α+β+2)/((n+1)*(n+α+β+1)*(2n+α+β))

(ch::Jacobi)(x::S) where {T,S} = orthogonal_polyval(ch, x)

# compute <Jn, Jn>
function norm2(::Type{Jacobi{α, β, T}}, n) where{α, β, T}
    α > -1 && β > -1 || throw(ArgumentError("α, β > -1 is necessary"))
    2^(α+β+1)/(2n+α+β+1) * (gamma(n+α+1) *  gamma(n + β +1))/(gamma(n+α+β+1)*gamma(n+1))
end
