## https://arxiv.org/pdf/1101.1798.pdf
"""
    Krawtchouk{N,T}

A family of discrete orthogonal polynomials for each N. 

We follow the parameterization  of  [Coleman](https://arxiv.org/pdf/1101.1798.pdf), notably with a parameter `s > 2`, not a parameter `0 < p  < 1`. With this  parameterization, the nodes are `xs=0,1,...,N` and the weight function is  `w_k=binomial(N,k)⋅(s-1)^k`.

"""
struct Krawtchouk{N, s, T<:Number} <: DiscreteOrthogonalPolynomial{T}
coeffs::Vector{T}
var::Symbol
    function Krawtchouk{N,s,T}(coeffs::AbstractVector{S}, var::Symbol=:x)  where {N, s, T <: Number, S <: Number}
        M = length(coeffs)
        M > N + 1 && throw(ArgumentError("Only degree $N or less polynomials may be specified"))
        new{N, s, T}(coeffs, var)
    end
end

@register2 Krawtchouk

export Krawtchouk

(ch::Krawtchouk)(x::S) where {S} = orthogonal_polyval(ch, x)

Polynomials.domain(::Type{<:Krawtchouk{N,s}}) where {N,s} = Polynomials.Interval(0, N)
basis_symbol(::Type{<:Krawtchouk}) = "K"
weight_function(w::Type{<:Krawtchouk{N,s}}) where {N,s} = k -> binomial(N,k) * (s-1)^k
###generating_function(w::Type{<:Krawtchouk{N,p}}) where {N,p} = (t,x) -> (1+t)^(N-x)*(1-t)^x


# Thm, 4.1
An(w::Type{<:Krawtchouk{N,s}}, n) where {N, s} = -s/(n+1)
Bn(w::Type{<:Krawtchouk{N,s}}, n) where {N, s} = (n + (s-1)*(N-n)) / (n+1)
Cn(w::Type{<:Krawtchouk{N,s}}, n) where {N, s} = -(s-1)*(N-n+1)/(n+1)

# prop 3.1
norm2(::Type{<:Krawtchouk{N,s}}, k) where {N, s} =  binomial(N,k) * s^N * (s-1)^k


function Base.:+(p1::Krawtchouk{N,s, T}, p2::Krawtchouk{M,s,S}) where {N,s,T,M,S}

    p1.var != p2.var && throw(ArgumentError("Polynomials must have same variable"))
    if N == M
        cs = [p1[i] + p2[i] for i = 0:N]
        return Krawtchouk{N,s}(cs, p1.var)
    else
        p, q = N >= M ? (p1, p2) : (p2, p1)
        P = typeof(p)
        pp = convert(P, convert(Polynomial,p) + convert(Polynomial, q))
        return P(coeffs(pp), p1.var)
    end
end

function Base.:*(p1::Krawtchouk{N,s,T}, p2::Krawtchouk{M,s,S}) where {N,T,M,S,s}
    R = promote_type(T,S)
    P = Krawtchouk{N+M-1, s, R}
    convert(P, convert(Polynomial, p1) * convert(Polynomial, p2))
end

