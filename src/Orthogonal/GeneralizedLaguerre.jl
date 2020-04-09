"""
   GeneralizedLaguerre{α, T <: Number}

The generalized Laguerre family have weight function `x^α * exp(-x)` over the domain `[0, oo)`. The parameter `α` is specified through the constructor.

```jldoctest
julia> p = GeneralizedLaguerre{1/2}([1,2,3])
GeneralizedLaguerre(1⋅L^(0.5)_0(x) + 2⋅L^(0.5)_1(x) + 3⋅L^(0.5)_2(x))

julia> convert(Polynomial, p)

Polynomial(9.625 - 9.5*x + 1.5*x^2)
```



"""
struct GeneralizedLaguerre{α, T <: Number} <: OrthogonalPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function GeneralizedLaguerre{α, T}(coeffs::AbstractVector{T}, var::Symbol=:x) where {α, T <: Number}
        α <= -1 && throw(ArgumentError("α > -1 is necessary"))
        length(coeffs) == 0 && return new{α, T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{α, T}(coeffs[1:last], var)
    end
end

export GeneralizedLaguerre

@register1 GeneralizedLaguerre


basis_symbol(::Type{<:GeneralizedLaguerre{α}}) where {α} = "L^($α)"

Polynomials.domain(::Type{<:GeneralizedLaguerre}) = Polynomials.Interval(0, Inf)



weight_function(::Type{<:GeneralizedLaguerre{α}}) where {α} = x  -> x^α * exp(-x)
generating_function(::Type{<:GeneralizedLaguerre{α}}) where {α} = (t,x) -> begin
    exp(-t*x/(1-t))/(1-t)^(α-1)
end


An(::Type{<:GeneralizedLaguerre{α}}, n) where {α} = -one(α)/(n+1)
Bn(::Type{<:GeneralizedLaguerre{α}}, n) where {α} = (2n + 1 +  α)/(n+1)
Cn(::Type{<:GeneralizedLaguerre{α}}, n) where {α} = -(n + α)/(n+1)

# compute <Jn, Jn>
function norm2(::Type{<:GeneralizedLaguerre{α}}, n) where{α}
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

# compute generalized binomial coefficient
function  generalized_binomial(N,K)
    tot = 1/1
    for k in K:-1:1
        tot /= k
        tot *= N
        N -= 1
    end
    tot
end
