struct GeneralizedLaguerre{α, T <: Number} <: OrthogonalPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function GeneralizedLaguerre{α, T}(coeffs::AbstractVector{T}, var::Symbol=:x) where {α, T <: Number}
        length(coeffs) == 0 && return new{α, T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{α, T}(coeffs[1:last], var)
    end
    function GeneralizedLaguerre{α}(coeffs::AbstractVector{T}, var::Symbol=:x) where {α, T <: Number}
        GeneralizedLaguerre{α, T}(coeffs, var)
    end
end

export GeneralizedLaguerre

Base.convert(::Type{P}, p::P) where {P <: GeneralizedLaguerre} =  p
Base.promote_rule(::Type{GeneralizedLaguerre{α, T}}, ::Type{GeneralizedLaguerre{α, S}}) where {α, T, S} = GeneralizedLaguerre{α, promote_type(T,S)}
Base.promote_rule(::Type{GeneralizedLaguerre{α, T}}, ::Type{S}) where {α, T, S <: Number} = GeneralizedLaguerre{α, promote_type(T,S)}
GeneralizedLaguerre{α}(n::Number, var=:x) where  {α}= GeneralizedLaguerre{α}([n], var)
GeneralizedLaguerre{α, T}(n::S, var=:x) where {α, T, S <: Number} = GeneralizedLaguerre{α, T}(T[n], var)
GeneralizedLaguerre{α, T}(xs::AbstractVector{S}, var=:x) where {α, T, S <: Number} = GeneralizedLaguerre{α, T}(T.(xs), var)

basis_symbol(::Type{GeneralizedLaguerre{α, T}}) where {α, T} = "L^($α)"

weight_function(::Type{GeneralizedLaguerre{α, T}}) where {α, T} = x^α * exp(-x)
generating_function(::Type{GeneralizedLaguerre{α, T}}) where {α, T} = (t,x) -> begin
    exp(-t*x/(1-t))/(1-t)^(α-1)
end

Polynomials.domain(::Type{<:GeneralizedLaguerre}) = Polynomials.Interval(0, Inf)
Polynomials.variable(::Type{GeneralizedLaguerre{α, T}}, var::Polynomials.SymbolLike=:x) where {α, T} =
    GeneralizedLaguerre{α, T}([(1+α), -1])

An(::Type{GeneralizedLaguerre{α, T}}, n) where {α, T} = -1/(n+1)
Bn(::Type{GeneralizedLaguerre{α, T}}, n) where {α, T} = 1 + α/(n+1)
Cn(::Type{GeneralizedLaguerre{α, T}}, n) where {α, T} = -(1 + (α - 1)/(k+1))

# compute <Jn, Jn>
function norm2(::Type{GeneralizedLaguerre{α, T}}, n) where{α, T}
    gamma(n + α + 1) / gamma(n+1)
end

(ch::GeneralizedLaguerre)(x::S) where {T,S} = orthogonal_polyval(ch, x)

function Base.convert(P::Type{GeneralizedLaguerre{β, T}}, p::GeneralizedLaguerre{α, S}) where {α, β, T, S}
    # use L(α,n) = sum(choose(α- β + n-1+1, n-1) L(β,i))
    d = degree(p)
    R = eltype(one(S)/1)
    qs = zeros(R, d+1)
    for i in 0:d
        for  j in i:d
            N = ( α - β + (j - i) - 1 )
            K = j - i
            qs[i+1] += p[j] * generalized_binomial(N, K)
        end
    end
    P(qs, p.var)
end

# compute generlized binomial coefficient
function  generalized_binomial(N,K)
    tot = 1/1
    for k in K:-1:1
        tot /= k
        tot *= N
        N -= 1
    end
    tot
end
