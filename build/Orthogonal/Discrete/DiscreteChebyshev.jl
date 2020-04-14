"""
    DiscreteChebyshev{N,T}

For each `N`, a family of discrete `N+1` orthogonal polynomials with `xs=0,1,...,N-1` and weight function `w_k=1`. 

# Example
```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> N = 8
8

julia> P = DiscreteChebyshev{N, Float64}
DiscreteChebyshev{8,Float64}

julia> p4, p5 = Polynomials.basis.(P, (4,5))
(DiscreteChebyshev(1.0⋅e_4(x)), DiscreteChebyshev(1.0⋅e_5(x)))

julia> p5(variable())
Polynomial(-2519.9999999999927 + 54228.0*x - 66150.0*x^2 + 26880.0*x^3 - 4410.0*x^4 + 252.0*x^5)

julia> abs(SpecialPolynomials.innerproduct(P, p4, p5)) <= sqrt(eps())
true

julia> SpecialPolynomials.innerproduct(P, p5, p5) ≈  N/(2*5+1)*prod(N^2-k^2 for k in 1:5)
true
```



"""
struct DiscreteChebyshev{N, T<:Number} <: DiscreteOrthogonalPolynomial{T}
coeffs::Vector{T}
var::Symbol
    function DiscreteChebyshev{N,T}(coeffs::AbstractVector{T}, var::Symbol=:x)  where {N, T <: Number}
        M = length(coeffs)
        # pad M to N
        M > N + 1 && throw(ArgumentError("Only degree $N or less polynomials can be specified"))
        new{N, T}(coeffs, var)
    end
end


export DiscreteChebyshev
@register1 DiscreteChebyshev

(ch::DiscreteChebyshev)(x::S) where {S} = orthogonal_polyval(ch, x)

Polynomials.domain(w::Type{<:DiscreteChebyshev{N}}) where {N} = Polynomials.Interval(0, N-1)
basis_symbol(::Type{<:DiscreteChebyshev}) = "C"
weight_function(w::Type{<:DiscreteChebyshev}) = x -> one(x)

## An,Bn,Cn for the orthogonal polynomials for the given weight function
An(w::Type{<:DiscreteChebyshev{N,T}}, n) where {N,T<:Integer} = 2(2n+1)//(n+1)
An(w::Type{<:DiscreteChebyshev{N}}, n) where {N} = 2(2n+1)/(n+1)
Bn(w::Type{<:DiscreteChebyshev{N,T}}, n) where {N,T<:Integer} = -(2n+1)*(N-1)//(n+1)
Bn(w::Type{<:DiscreteChebyshev{N}}, n) where {N} = -(2n+1)*(N-1)/(n+1)
Cn(w::Type{<:DiscreteChebyshev{N,T}}, n) where {N,T<:Integer} = -n*(N^2-n^2)//(n+1)
Cn(w::Type{<:DiscreteChebyshev{N}}, n) where {N} = -n*(N^2-n^2)/(n+1)

norm2(w::Type{<:DiscreteChebyshev{N}},n) where {N} = iszero(n) ? N/(2n+1) : N*prod(N^2-i^2 for i in 1:n)/(2n+1)


function Base.:+(p1::DiscreteChebyshev{N,T}, p2::DiscreteChebyshev{M,S}) where {N,T,M,S}
    p1.var != p2.var && throw(ArgumentError("Polynomials must have same variable"))
    if N == M
        cs = [p1[i] + p2[i] for i = 0:N]
        return DiscreteChebyshev(cs, p1.var)
    else
        p, q = N >= M ? (p1, p2) : (p2, p1)
        @show  convert(Polynomial,p) + convert(Polynomial, q)
        pp = convert(DiscreteChebyshev{N, promote_type(T,S)}, convert(Polynomial,p) + convert(Polynomial, q))
        return DiscreteChebyshev(coeffs(pp), p1.var)
    end
end

function Base.:*(p1::DiscreteChebyshev{N,T}, p2::DiscreteChebyshev{M,S}) where {N,T,M,S}
    R = promote_type(T,S)
    P = DiscreteChebyshev{N+M-1, R}
    convert(P, convert(Polynomial, p1) * convert(Polynomial, p2))
end
