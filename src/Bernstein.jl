abstract type AbstractBernstein{T} <: Polynomials.AbstractPolynomial{T} end

"""

    Bernstein{N, T}

A [Bernstein  polynomial](https://en.wikipedia.org/wiki/Bernstein_polynomial) is a polynomial expressed in terms of
Bernstein basic polynomials. For each degree, `n`, this is a set of `n+1` degree `n` polynomials of the form:
`beta_{ν, n} =  (ν choose n) x^ν  (1-x)^{n-ν}`, `0 ≤ x ≤ 1.`

The  `Bernstein{N,T}` type represents a polynomial of degree `N` or with a linear combination of the basis vectors using coefficients  of  type `T`.

```jldoctest
julia> p = Bernstein{3,Int}([0,0,1,0])
Bernstein(1⋅β(3, 2)(x))

julia> convert(Polynomial, p)
Polynomial(3*x^2 - 3*x^3)
```

"""
struct Bernstein{N, T<:Number} <: AbstractBernstein{T}
coeffs::Vector{T}
    var::Symbol
    function Bernstein{N, T}(coeffs::AbstractVector{T}, var::Symbol=:x)  where {N, T <: Number}
        M = length(coeffs)
        M > max(0,N)+1 && throw(ArgumentError("Only degree $N or less polynomials can be specified"))
        new{N, T}(T.(coeffs), var)
    end
end

export Bernstein
@register1 Bernstein

Bernstein(x::Number, var=:x) = convert(Bernstein, x)

# add default case with no N specified
function Bernstein(coeffs::AbstractVector{T}, var::Symbol=:x)  where {T <: Number}
    n = findlast(!iszero,coeffs)
    if n == nothing
        zero(Bernstein{0, T})
    else
        Bernstein{n-1, T}(coeffs[1:n], var)
    end
end

function Polynomials.showterm(io::IO, ::Type{Bernstein{N, T}}, pj::T, var, j, first::Bool, mimetype) where {N, T}
    iszero(pj) && return false
    !first &&  print(io, " ")
    print(io, Polynomials.hasneg(T) && Polynomials.isneg(pj) ? "- " :  (!first ? "+ " : ""))
    print(io, "$(abs(pj))⋅β($N, $j)($var)")
    return true
end


function (ch::AbstractBernstein{T})(x::S) where {T,S}
    # XXX make more efficient
    x ∉ domain(ch) && error("$x outside of domain")
    R = promote_type(T, S)
    length(ch) == 0 && return zero(R)
    convert(Polynomial{T}, ch)(x)
end



function Base.convert(P::Type{<:Polynomial}, ch::Bernstein{N,T}) where {N, T}
    N == -1 && return zero(Polynomial{T})
    out = P(zeros(T,1), ch.var)
    x = P([zero(T),one(T)], ch.var)
    @inbounds for (i,ai) in enumerate(coeffs(ch))
        nu = i - 1
        out += ai *  binomial(N, nu) * x^nu * (1-x)^(N-nu)
    end
    out

end

function Base.convert(C::Type{<:Bernstein{N, T}}, p::Polynomial) where {N, T}

    @assert degree(p) <= N

    R = eltype(one(T)/one(Int))
    cs = zeros(R, N+1)

    for (i, a_nu) in enumerate(coeffs(p))
        k = i - 1
        nk = binomial(N,k)
        for j in k:N
            cs[j+1] += a_nu/nk*binomial(j,k)
        end
    end
    Bernstein{N, R}(cs, p.var)
end

function Base.convert(::Type{<:AbstractBernstein}, p::Polynomial{T}) where {T}
    N = degree(p)
    convert(Bernstein{N, T}, p)
end

function Base.convert(P::Type{Bernstein{N,T}}, p::Bernstein{M, S}) where {N, T, M, S}
    @assert  N >=  M
    convert(P, convert(Polynomial{T},   p))
end

# create a basis vector
function Polynomials.basis(p::Bernstein{N,T}, k::Int) where {N, T}
    0 <= k <= N || throw(ArgumentError("k must be in 0 .. $N"))
    zs = zeros(T, N+1)
    zs[k+1] = one(T)
    Bernstein{N,T}(zs, p.var)
end

function bernstein(N, T, nu::Int, var::Symbol=:x)
    zs = zeros(T,  N+1)
    zs[nu+1] = one(T)
    Bernstein{N, T}(zs, var)
end

Polynomials.domain(::Type{<:AbstractBernstein}) = Polynomials.Interval(0, 1)
Polynomials.degree(p::Bernstein{N,T}) where {N, T} = degree(convert(Polynomial{T}, p))
Polynomials.variable(P::Type{<:Bernstein{N,T}}) where {N, T} = convert(P, Polynomial{T}([0,1]))
Base.one(P::Type{Bernstein{N,T}}) where {N,T} = Bernstein{N,T}(ones(T, N+1))
Base.one(P::Type{<:Bernstein{N}}) where {N} = Bernstein{N,Int}(ones(Int, N+1))
Base.one(P::Type{<:Bernstein})  = Bernstein{0,Int}([1])




function Base.:+(p1::Bernstein{N,T}, p2::Bernstein{M,S}) where {N, T, M, S}

    p1.var == p2.var || throw(ArgumentError("p1 and p2 must have the same symbol"))

    R = promote_type(T, S)

    if M == N
        return Bernstein{N, R}([p1[i] + p2[i] for i in 0:N], p1.var)
    else
        NN = max(M, N)
        q1 = convert(Polynomial{R}, p1)
        q2 = convert(Polynomial{R}, p2)
        return convert(Bernstein{NN, R}, q1+q2)
    end
end


function Base.:*(p1::Bernstein{N,T}, p2::Bernstein{M,S}) where {N,T, M,S}
    ## use b(n,k) * b(m,j) = choose(n,k)choose(m,j)/choose(n+m,k+j) b(n+m, k+j)

    p1.var == p2.var || throw(ArgumentError("p1 and p2 must have the same symbol"))

    R1 = promote_type(T,S)
    R = typeof(one(R1)/one(R1))

    c = zeros(R, N + M + 1)
    bnmi = 1 // (binomial(N+M,M))
    for i in 0:N
        for j in 0:M
            aij = binomial(N+M-(i+j),N-i) * binomial(i+j,i) * bnmi
            @inbounds c[i + j + 1] += aij * p1[i] * p2[j]
        end
    end

    Bernstein(c, p1.var)

end


function Polynomials.derivative(p::Bernstein{N, T}, order::Integer = 1) where {N, T}
    NN =  N  - 1
    cs = coeffs(p)
    NN < 0 &&  return p
    order ==  0 && return p
    csp = zeros(eltype(cs), NN+1)
    # nu = 0
    for (i,a_nu) in enumerate(cs)
        nu = i - 1
        nu > 0  && (csp[nu] += a_nu *  (NN+1))
        nu <= NN && (csp[nu+1] -= a_nu *  (NN+1))
    end
    csp
    pp = Bernstein{NN, T}(csp, p.var)
    derivative(pp, order-1)
end




function Polynomials.integrate(p::Bernstein{N, T}, C::S) where {N, T,S <: Number}
    R = promote_type(eltype(one(T) / 1), S)
    NN = N + 1
    cs = zeros(R, NN+1)

    @inbounds for (i,a_nu) in enumerate(coeffs(p))
        nu = i - 1
        for j in (nu+1):NN
            cs[j+1] += a_nu/(NN)
        end
    end

    Bernstein{NN, R}(cs, p.var) + C
end


function Base.divrem(num::Bernstein{N,T}, den::Bernstein{M,S}) where {N, T, M, S}
    p1 = convert(Polynomial{T}, num)
    p2 = convert(Polynomial{S}, den)
    q,r = divrem(p1, p2)
    convert.(Bernstein, (q,r))
end


function Polynomials.vander(P::Type{<:AbstractBernstein}, xs::AbstractVector{T}, n::Integer) where {T <: Number}
    N = length(xs) - 1 # xs = [x0, x1, ...]
    R = typeof(one(T)/one(T))
    V = zeros(R, N+1, n+1)
    for j in 0:n
        bnj = binomial(n,j)
        for i in 0:N
            x = xs[i+1]
            V[i+1, j+1] = bnj * x^j * (1-x)^(n-j)     #beta(n,j)(xj)
        end
    end
    V
end
