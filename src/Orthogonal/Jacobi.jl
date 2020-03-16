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

Base.convert(::Type{P}, p::P) where {P <: Jacobi} =  p
Base.convert(P::Type{<:Jacobi}, p::Jacobi) =  P(p.coeffs, p.var)
Base.promote_rule(::Type{Jacobi{α, β, T}}, ::Type{Jacobi{α, β, S}}) where {α, β, T, S} = Jacobi{α, β, promote_type(T,S)}
Base.promote_rule(::Type{Jacobi{α, β, T}}, ::Type{S}) where {α, β, T, S <: Number} = Jacobi{α, β, promote_type(T,S)}
Jacobi{α, β}(n::Number, var=:x) where  {α, β}= Jacobi{α, β}([n], var)
Jacobi{α, β, T}(n::S, var=:x) where {α, β, T, S <: Number} = Jacobi{α, β, T}(T[n], var)
Jacobi{α, β, T}(xs::AbstractVector{S}, var=:x) where {α, β, T, S <: Number} = Jacobi{α, β, T}(T.(xs), var)



function Polynomials.showterm(io::IO, ::Type{Jacobi{α, β, T}}, pj::T, var, j, first::Bool, mimetype) where {α, β, T}
    iszero(pj) && return false
    !first &&  print(io, " ")
    print(io, Polynomials.hasneg(T) && pj < 0 ? "- " :  (!first ? "+ " : ""))
    print(io, "$(abs(pj))⋅J^(α, β)_$j($var)")
    return true
end


weight_function(::Type{Jacobi{α, β, T}}) where {α, β, T} = x -> (1-x)^α *  (1+x)^β
generating_function(::Type{Jacobi{α, β, T}}) where {α, β, T} = (t,x) -> begin
    R = sqrt(1 - 2x*t+t^2)
    2^(α + β) * 1/R  * (1 - t + R)^(-α) * (1 + t + R)^(-β)
end

Polynomials.domain(::Type{<:Jacobi}) = Polynomials.Interval(-1, 1)
Polynomials.variable(::Type{Jacobi{α, β, T}}, var::Polynomials.SymbolLike=:x) where {α, β, T} =
    P([(α+1) - (α+β+2)/2, (α+β+2)/2], var)

An(::Type{Jacobi{α, β, T}}, n) where {α, β, T} = (2n+α+β-1)*(2n+α+β)/(2n*(n+α+β))
Bn(::Type{Jacobi{α, β, T}}, n) where {α, β, T} = (α^2-β^2)*(2n+α+β-1)/(2n*(n+α+β)*(2n+α+β-2))
Cn(::Type{Jacobi{α, β, T}}, n) where {α, β, T} = -(n+α-1)*(n+β-1)*(2n+α+β)/(n*(n+α+β)*(2n+α+β-2))
P0(::Type{Jacobi{α, β, T}}, x) where {α, β, T} = 1
P1(::Type{Jacobi{α, β, T}}, x) where {α, β, T} = (α+1) + (α+β+2)*(x-1)/2

# compute <Jn, Jn>
function Base.abs2(::Type{Jacobi{α, β, T}}, n) where{α, β, T}
    α > -1 && β > -1 || throw(ArgumentError("α, β > -1 is necessary"))
    2^(α+β+1)/(2n+α+β+1) * (gamma(n+α+1) *  gamma(n + β +1))/(gamma(n+α+β+1)*gamma(n+1))
end

(ch::Jacobi)(x::S) where {T,S} = orthogonal_polyval(ch, x)

## This is going  to be slow!
function Base.convert(P::Type{J}, p::Polynomial) where {J <: Jacobi}
    d = degree(p)
    R = eltype(one(eltype(p))/1)
    qs = zeros(R, d+1)
    for i in 0:d
        qs[i+1] = sum(p[j] * _jacobi_lambda(J, j, i) for j in i:d)
    end
    J(qs, p.var)
end

function _jacobi_lambda(J, n, j)
   p = convert(Polynomial, Polynomials.basis(J, j))
   innerproduct(J, x->x^n,  p) / abs2(J, j)
end
