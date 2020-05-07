abstract type AbstractBernstein{T} <: Polynomials.AbstractPolynomial{T} end

"""

    Bernstein{N, T}

A [Bernstein  polynomial](https://en.wikipedia.org/wiki/Bernstein_polynomial) is a polynomial expressed in terms of
Bernstein basic polynomials. For each degree, `n`, this is a set of `n+1` degree `n` polynomials of the form:
`β_{ν, n} =  (ν choose n) x^ν  (1-x)^{n-ν}`, `0 ≤ x ≤ 1.`

The  `Bernstein{N,T}` type represents a polynomial of degree `N` or with a linear combination of the basis vectors using coefficients  of  type `T`.

```jldoctest Bernstein
julia> using Polynomials, SpecialPolynomials

julia> p = Bernstein{3,Int}([0,0,1,0])
Bernstein(1⋅β(3, 2)(x))

julia> convert(Polynomial, p)
Polynomials.Polynomial(3*x^2 - 3*x^3)
```

"""
struct Bernstein{N, T} <: AbstractBernstein{T}
    coeffs::Vector{T}
    var::Symbol
    function Bernstein{N, T}(coeffs::AbstractVector{T}, var::Polynomials.SymbolLike=:x)  where {N, T}
        M = length(coeffs)
        M > max(0,N)+1 && throw(ArgumentError("Only degree $N or less polynomials can be specified"))
        new{N,T}(coeffs, Symbol(var))
    end
end

export Bernstein
Polynomials.@register1 Bernstein

# add default case with no N specified
function Bernstein(coeffs::AbstractVector{T}, var::Symbol=:x)  where {T}
    N = length(coeffs)-1
    Bernstein{N, T}(coeffs, var)
end

function Polynomials.basis(P::Type{<:Bernstein{N,T}},i::Int,var=:x) where {N,T}
    cs = zeros(T, N+1)
    cs[i+1] = 1
    P(cs,var)
end
Polynomials.basis(P::Type{<:Bernstein{N}},i::Int,var=:x) where {N} = Polynomials.basis(Bernstein{N,Int},i,var)

function Base.getindex(p::Bernstein{N,T}, i::Int) where {N, T}
    i < 0 && throw(ArgumentError(""))
    i > N && return zero(eltype(T))
    coeffs(p)[i+1]
end

function Polynomials.showterm(io::IO, ::Type{Bernstein{N, T}}, pj::T, var, j, first::Bool, mimetype) where {N, T}
    iszero(pj) && return false
    !first &&  print(io, " ")
    print(io, Polynomials.hasneg(T) && Polynomials.isneg(pj) ? "- " :  (!first ? "+ " : ""))
    print(io, "$(abs.(pj))⋅β($N, $j)($var)")
    return true
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
Polynomials.variable(P::Type{<:Bernstein{N}}, var::Polynomials.SymbolLike=:x) where {N} = Bernstein([1 - i/N for i in N:-1:0], var)

## need to do zero and one differently,as coefficients have length N+1
Base.one(P::Type{Bernstein{N,T}}, var=:x) where {N,T} = Bernstein{N,T}(ones(T, N+1))
Base.one(P::Type{<:Bernstein{N}}, var=:x) where {N} = one(Bernstein{N,Int},var)
Base.one(P::Type{<:Bernstein}, var=:x)  = one(Bernstein{0,Int}, var)
Base.one(p::P) where {P <: Bernstein} = one(P,p.var)

Base.zero(::Type{<:Bernstein{N,T}}, var=:x) where {N,T} = Bernstein{N,T}(zeros(T, N+1), var)
Base.zero(::Type{<:Bernstein{N}}, var=:x) where {N} = zero(Bernstein{N,Int}, var)
Base.zero(::Type{<:Bernstein}, var=:x)  = zero(Bernstein{0,Int}, var)
Base.zero(p::P) where {P <: Bernstein} = zero(P,p.var)

# we use the compensated one below for increased accuracy, but is it really necessary?
function simple_deCasteljau_polyval(p::Bernstein{N,T}, t) where {N,T}
    bs = coeffs(p).*one(t)
    r = 1-t
    for j in 1:N
        for i = 0:N-j
            ii = i + 1
            bs[ii] = bs[ii]*r + bs[ii+1]*t
        end
    end
    bs[1]
end

## Evaluation https://doi.org/10.1016/j.camwa.2010.05.021
## Accurate evaluation of a polynomial and its derivative in Bernstein form☆
## Hao Jiang, Shengguo Li, Lizhi Cheng, and Fang Su
function deCasteljau_polyval(p::Bernstein{N,T}, t) where {N, T}
    b̂ = coeffs(p) .* one(t)
    errb̂ = 0 .* b̂ #zeros(eltype(b̂), N+1) 
    r̂,ρ = twoSum(1, -t)
    for j in 1:N
        for ii in 0:(N-j)
            i = ii+1
            s, π̂ = twoProd(b̂[i], r̂)
            v, σ̂ = twoProd(b̂[i+1], t)
            b̂[i],β̂ = twoSum(s, v)
            ŵ = π̂ .+ σ̂ .+ β̂ .+ (b̂[i] .* ρ)
            errb̂[i] = errb̂[i] * r̂ + (errb̂[i+1] * t + ŵ)
        end
    end
    b̂[1]+ errb̂[1]

end



# compensated form to  find p'(t)
function deCasteljau_derivative_val(p::Bernstein{N,T}, t) where {N,T}
    R = eltype(one(T)*one(eltype(t)))
    r̂,ρ  = twoSum(one(R), -t)
    Δb = zeros(R, N).*one(t)
    errΔb = similar(Δb) #zeros(R, N).*one(t)
    for ii=0:N-1
        i = ii+1
        Δb[i], errΔb[i] = twoSum(p[ii+1], -p[ii])
    end
    for j in 1:N-1
        for ii in 0:N-1-j
            i = ii+1
            s, π̂ = twoProd(Δb[i], r̂)
            v, σ̂ = twoProd(Δb[i+1], t)
            Δb[i], β̂ = twoSum(s, v)
            ŵ = π̂ + σ̂ + β̂ + (Δb[i] * ρ)
            errΔb[i] = errΔb[i] * r̂ + (errΔb[i+1] * t + ŵ)
        end
    end
    (Δb[1] + errΔb[1]) * N
end

# Moller-Knuth from Exact computations with approximate arithmetic
# http://perso.ens-lyon.fr/jean-michel.muller/
function twoSum(a::Number, b::Number)
    s = a + b
    â = s - b
    b̂ = s - â
    δ_a = a - â
    δ_b = b - b̂
    s, δ_a + δ_b
end
twoSum(x, y) = (x + y, zero(x)+zero(y))

function twoProd(a::Number, b::Number)
    ab = a * b
    ab, fma(a,b,-ab)
end
twoProd(x, y)  = x*y, 0


# poly evaluation
(p::AbstractBernstein{T})(x::S) where {T,S} = deCasteljau_polyval(p, x)



#
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
    bnmi = 1 / (binomial(N+M,M))
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

    #  use p'(t) = n ∑ (b_{i+1} - b_i) B_i,n-1
    dp =  Bernstein{NN, T}(N*diff(cs), p.var)
    
    order > 1 ? derivative(dp, order-1) : dp
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
