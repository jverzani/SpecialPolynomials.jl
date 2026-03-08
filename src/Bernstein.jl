#abstract type AbstractBernstein{T,X} <: AbstractSpecialPolynomial{T,X} end

# XXX TODO
# * the constraint 𝐍 < N - 1 isn't checked
# * methods written from Bernstein{𝐍,T,X} but should be {𝐍, B<:BernsteinBasis{𝐍}, AbstractUnivariatePolynomial{B}}
struct BernsteinBasis{𝐍} <: AbstractBasis end

"""

    Bernstein{N, T, X}

A [Bernstein  polynomial](https://en.wikipedia.org/wiki/Bernstein_polynomial) is a polynomial expressed in terms of
Bernstein basic polynomials. For each degree, `𝐍`, this is a set of `𝐍+1` degree `𝐍` polynomials of the form:
`β_{𝐍,ν} =  (ν choose 𝐍) x^ν  (1-x)^{𝐍-ν}`, `0 ≤ x ≤ 1.`

The `Bernstein{𝐍,T}` type represents a polynomial of degree `𝐍` or
less with a linear combination of the basis vectors using coefficients
of type `T`.

```jldoctest Bernstein
julia> using Polynomials, SpecialPolynomials

julia> p = basis(Bernstein{3},  2)
Bernstein(1⋅β₃,₂(x))

julia> convert(Polynomial, p)
Polynomial(3.0*x^2 - 3.0*x^3)
```

!!! note
    [StaticUnivariatePolynomials](https://github.com/tkoolen/StaticUnivariatePolynomials.jl) offers a  more  performant version.

"""
Bernstein{𝐍} = MutableDensePolynomial{BernsteinBasis{𝐍}} where {𝐍}
Polynomials._typealias(::Type{P}) where {P<:Bernstein{𝐍}} where {𝐍} = "Bernstein"

# struct Bernstein{𝐍,T,X} <: AbstractBernstein{T,X}
#     coeffs::Vector{T}
#     function Bernstein{𝐍,T,X}(coeffs::AbstractVector{T}) where {𝐍,T,X}
#         N = findlast(!iscoeffzero, coeffs)
#         N == nothing && return new{𝐍,T,X}(zeros(T, 0))
#         (N > 𝐍 + 1) && throw(ArgumentError("Wrong length for coefficients"))
#         return new{𝐍,T,X}(coeffs[1:N])
#     end
# end

export Bernstein
#Polynomials.@registerN Bernstein 𝐍
#Polynomials.:⟒(::Type{<:Bernstein{𝐍}}) where {𝐍} = Bernstein{𝐍}

function Polynomials.showterm(
    io::IO,
    ::Type{<:AbstractUnivariatePolynomial{BernsteinBasis{N},T,X}},
    pj::T,
    var,
    j,
    first::Bool,
    mimetype,
) where {N,T,X}
    iscoeffzero(pj) && return false
    !first && print(io, " ")
    print(io, Polynomials.hasneg(T) && Polynomials.isneg(pj) ? "- " : (!first ? "+ " : ""))
    print(io, "$(abs.(pj))⋅β")
    Polynomials.unicode_subscript(io, N)
    print(io, ",")
    Polynomials.unicode_subscript(io, j)
    print(io, "($var)")
    return true
end

# add default case with no N specified
function Bernstein(coeffs::AbstractVector{T}, var::Symbol=:x) where {T}
    l = findlast(!iscoeffzero, coeffs)
    𝐍 = l == nothing ? -1 : l - 1
    Bernstein{𝐍,T,Symbol(var)}(coeffs[1:(𝐍 + 1)])
end

# promotion between family
function Base.promote(
    p::P,
    q::Q,
) where {𝐍,T,X,P<:Bernstein{𝐍,T,X},𝐌,S,Y,Q<:Bernstein{𝐌,S,Y}}
    X == Y || throw(ArgumentError("different indeterminates"))
    MN = max(𝐌, 𝐍)
    R = promote_type(T, S)
    PQ = Bernstein{MN,R}
    convert(PQ, p), convert(PQ, q)
end

# promote up
function Base.convert(P::Type{<:Bernstein{𝐍}}, p::Bernstein{𝐌,S}) where {𝐍,𝐌,S}
    @assert 𝐍 >= 𝐌
    𝐑 = 𝐍 - 𝐌
    R = eltype(one(promote_type(eltype(P), S)) / 1)
    cs = zeros(R, 𝐍 + 1)
    for j in 0:𝐌
        pⱼ = p[j]
        iscoeffzero(pⱼ) && continue
        for i in j:(j + 𝐑)
            cs[1 + i] += pⱼ * binomial(𝐌, j) * binomial(𝐑, i - j) / binomial(𝐍, i)
        end
    end
    X = Polynomials.indeterminate(P, p)
    Bernstein{𝐍,R,X}(cs)
end

function Base.convert(
    P::Type{<:Polynomials.StandardBasisPolynomial},
    p::Bernstein{N,T},
) where {N,T}
    x = variable(P)
    p(x)
end

function Base.convert(C::Type{<:Bernstein{𝐍,T}}, p::Polynomial) where {𝐍,T}
    @assert degree(p) <= 𝐍

    R = eltype(one(T) / one(Int))
    cs = zeros(R, 𝐍 + 1)

    for (i, a_nu) in enumerate(coeffs(p))
        k = i - 1
        nk = binomial(𝐍, k)
        for j in k:𝐍
            cs[j + 1] += a_nu / nk * binomial(j, k)
        end
    end

    X = Polynomials.indeterminate(C, p)
    Bernstein{𝐍,R,X}(cs)
end

function Base.convert(P::Type{<:Bernstein}, p::Polynomial{T}) where {T}
    𝐍 = degree(p)
    X = Polynomials.indeterminate(P, p)
    convert(Bernstein{𝐍,T,X}, p)
end

Polynomials.domain(::Type{<:BernsteinBasis}) = Polynomials.Interval(0, 1)
Polynomials.domain(
    ::Type{P},
) where {𝐍,B<:BernsteinBasis{𝐍},P<:AbstractUnivariatePolynomial{B}} =
    Polynomials.Interval(0, 1)
Polynomials.domain(::Type{Bernstein}) = Polynomials.Interval(0, 1) # no N
Polynomials.degree(p::Bernstein{N,T,X}) where {N,T,X} = degree(convert(Polynomial{T,X}, p))

Polynomials.constantterm(p::Bernstein) = p[0] # β.,ᵢ = 0 + ... when i > 0

## need to do one  and variable differently,as coefficients have length N+1
Base.one(P::Type{Bernstein{N,T,X}}) where {N,T,X} = Bernstein{N,T,X}(ones(T, N + 1))
# Base.one(P::Type{Bernstein{N,T}}, var::Polynomials.SymbolLike) where {N,T} =
#     Bernstein{N,T,Symbol(var)}(ones(T, N + 1))
# function Base.one(P::Type{<:Bernstein{N}}, var::Polynomials.SymbolLike=:x) where {N}
#     N′ = max(N, 0)
#     one(Bernstein{N′,eltype(P),Symbol(var)})
# end

function Polynomials.variable(P::Type{Bernstein{𝐍,T,X}}) where {𝐍,T,X}
    𝐍 <= 0 && throw(ArgumentError("Need  𝐍 ≥ 1"))
    R = typeof(one(T) / 1)
    Bernstein{𝐍,R,X}([1 - (i * one(T)) / 𝐍 for i in 𝐍:-1:0])
end

function Polynomials.variable(
    P::Type{<:Bernstein{𝐍}},
    var::Polynomials.SymbolLike=:x,
) where {𝐍}
    S = eltype(P)
    variable(Bernstein{𝐍,S,Symbol(var)})
end
Polynomials.variable(P::Type{Bernstein}, var::Polynomials.SymbolLike=:x) =
    variable(Bernstein{1,eltype(P),Symbol(var)})

function basis(::Type{P}, k::Int) where {𝐍,T,X,P<:Bernstein{𝐍,T,X}}
    zs = zeros(Int, k + 1)
    zs[end] = 1
    P(zs)
end

function basis(P::Type{<:Bernstein}, k::Int, _var::Polynomials.SymbolLike=:x; var=_var)
    zs = zeros(Int, k + 1)
    zs[end] = 1
    P(zs, var)
end

# poly evaluation
Polynomials.evalpoly(x, p::Bernstein) = simple_deCasteljau_polyval(p, x)

# we could  use the compensated one
# from https://doi.org/10.1016/j.camwa.2010.05.021
# Accurate evaluation of a polynomial and its derivative in Bernstein form
# Hao Jiang, Shengguo Li, Lizhi Cheng, and Fang Su
#
# for increased accuracy, but is it really necessary here?
function simple_deCasteljau_polyval(p::Bernstein{𝐍,T,X}, t::S) where {𝐍,T,X,S}
    p₀ = p[0]
    R = eltype(p₀ * t)
    zR = p₀ * zero(t)

    R = promote_type(T, S)
    𝐍 == -1 && return zR
    𝐍 == 0 && return p[0] * one(t)
    n = length(coeffs(p))
    iscoeffzero(n) && return zR

    bs = [zR for _ in 1:(𝐍 + 1)]
    for i in eachindex(coeffs(p))
        bs[i] = p[i - 1]
    end

    r = 1 - t
    for j in 1:𝐍
        for i in 0:(𝐍 - j)
            ii = i + 1
            bs[ii] = bs[ii] * r + bs[ii + 1] * t
        end
    end
    bs[1]
end

# function deCasteljau_polyval(p::Bernstein{N,T}, t) where {N, T}
#     b̂ = coeffs(p) .* one(t)
#     errb̂ = 0 .* b̂ #zeros(eltype(b̂), N+1)
#     r̂,ρ = twoSum(1, -t)
#     for j in 1:N
#         for ii in 0:(N-j)
#             i = ii+1
#             s, π̂ = twoProd(b̂[i], r̂)
#             v, σ̂ = twoProd(b̂[i+1], t)
#             b̂[i],β̂ = twoSum(s, v)
#             ŵ = π̂ .+ σ̂ .+ β̂ .+ (b̂[i] .* ρ)
#             errb̂[i] = errb̂[i] * r̂ + (errb̂[i+1] * t + ŵ)
#         end
#     end
#     b̂[1]+ errb̂[1]

# end

# # compensated form to  find p'(t)
# function deCasteljau_derivative_val(p::Bernstein{N,T}, t) where {N,T}
#     R = eltype(one(T)*one(eltype(t)))
#     r̂,ρ  = twoSum(one(R), -t)
#     Δb = zeros(R, N).*one(t)
#     errΔb = similar(Δb) #zeros(R, N).*one(t)
#     for ii=0:N-1
#         i = ii+1
#         Δb[i], errΔb[i] = twoSum(p[ii+1], -p[ii])
#     end
#     for j in 1:N-1
#         for ii in 0:N-1-j
#             i = ii+1
#             s, π̂ = twoProd(Δb[i], r̂)
#             v, σ̂ = twoProd(Δb[i+1], t)
#             Δb[i], β̂ = twoSum(s, v)
#             ŵ = π̂ + σ̂ + β̂ + (Δb[i] * ρ)
#             errΔb[i] = errΔb[i] * r̂ + (errΔb[i+1] * t + ŵ)
#         end
#     end
#     (Δb[1] + errΔb[1]) * N
# end

# # Moller-Knuth from Exact computations with approximate arithmetic
# # http://perso.ens-lyon.fr/jean-michel.muller/
# function twoSum(a::Number, b::Number)
#     s = a + b
#     â = s - b
#     b̂ = s - â
#     δ_a = a - â
#     δ_b = b - b̂
#     s, δ_a + δ_b
# end
# twoSum(x, y) = (x + y, zero(x)+zero(y))

# function twoProd(a::Number, b::Number)
#     ab = a * b
#     ab, fma(a,b,-ab)
# end
# twoProd(x, y)  = x*y, 0

# not needed, but speeds things along
function Base.:+(p1::P, p2::Q) where {𝐍,T,X,P<:Bernstein{𝐍,T,X},𝐌,S,Q<:Bernstein{𝐌,S,X}}
    p, q = promote(p1, p2)
    p + q
end

function Base.:+(p::P, q::P) where {𝐍,T,X,P<:Bernstein{𝐍,T,X}}
    return Bernstein{𝐍}([p[i] + q[i] for i in 0:𝐍], X)
end

# no promote(p1,p2)  called here
function Base.:*(p::P, q::Q) where {𝐍,T,X,P<:Bernstein{𝐍,T,X},𝐌,S,Y,Q<:Bernstein{𝐌,S,Y}}
    ## use b(n,k) * b(m,j) = choose(n,k)choose(m,j)/choose(n+m,k+j) b(n+m, k+j)

    isconstant(p) && return q * constantterm(p)
    isconstant(q) && return p * constantterm(q)
    assert_same_variable(p, q) ||
        throw(ArgumentError("`p` and `q` have different indeterminate"))

    R = typeof(one(promote_type(T, S)) / 1)
    x = variable(Polynomial{R})
    return convert(Bernstein, p(x) * q(x))

    ## multiplication  through conversion  is as  efficient
    # n,  m =  degree(p), degree(q)
    # 𝐍 = n + m
    # cs = zeros(R, 𝐍 + 1)

    # for i in 0:𝐍
    #     for j in max(0,i-m):min(n, i)
    #         cs[1+i] += p[j] * q[i-j] * (one(R) * binomial(n,j)  *  binomial(m,i-j))
    #     end
    #     cs[1+i] /= binomial(m+n, i)
    # end

    # Bernstein{𝐍,R}(cs, p.var)

end

#  use p'(t) = n ∑ (b_{i+1} - b_i) B_i,n-1
function Polynomials.derivative(p::Bernstein{𝐍,T,X}) where {𝐍,T,X}
    cs = coeffs(p)
    𝐍 < -1 && return p
    iscoeffzero(𝐍) && return zero(p)

    cs = zeros(T, 𝐍)
    dp = Bernstein{𝐍 - 1,T,X}(𝐍 * diff(coeffs(p)))
    dp
end

# use ∫b_n,nu  = 1/(n+1)  * ∑_{nu+1}^{n+1} b_n+1, j
function Polynomials.integrate(p::Bernstein{𝐍,T,X}) where {𝐍,T,X}
    R = eltype(one(T) / 1)
    𝐍𝐍 = 𝐍 + 1

    cs = zeros(R, 𝐍𝐍 + 1)

    @inbounds for ν in eachindex(p)
        for j in (ν + 1):𝐍𝐍
            cs[1 + j] += p[ν] / 𝐍𝐍
        end
    end

    ∫p = Bernstein{𝐍𝐍,R,X}(cs)
    return ∫p
end

#  could do p 33 of http://mae.engr.ucdavis.edu/~farouki/bernstein.pdf
function Base.divrem(num::Bernstein{𝐍,T}, den::Bernstein{𝐌,S}) where {𝐍,T,𝐌,S}
    R = eltype((one(T) + one(S)) / 1)
    p1 = convert(Polynomial{R}, num)
    p2 = convert(Polynomial{R}, den)
    q, r = divrem(p1, p2)
    convert.(Bernstein, (q, r))
end

# Replace me. Doesn't seem necessary and poor stability of conversion.
# cf https://arxiv.org/pdf/1910.01998.pdf for one possible way, though for approximate GCD
function Base.gcd(num::Bernstein{𝐍,T}, den::Bernstein{𝐌,S}) where {𝐍,T,𝐌,S}
    ps = convert.(Polynomial, (num, den))
    qs = gcd(ps...)
    convert.(Bernstein, qs)
end

function Polynomials.vander(
    P::Type{<:Bernstein},
    xs::AbstractVector{T},
    n::Integer,
) where {T<:Number}
    N = length(xs) - 1 # xs = [x0, x1, ...]
    R = typeof(one(T) / one(T))
    V = zeros(R, N + 1, n + 1)
    for j in 0:n
        bnj = binomial(n, j)
        for i in 0:N
            x = xs[i + 1]
            V[i + 1, j + 1] = bnj * x^j * (1 - x)^(n - j)     #beta(n,j)(xj)
        end
    end
    V
end

# direct computation of roots using numerical stable method of
# Guðbjörn and Jónsson as
# noted by Corless, Sevyeri  in §2.1 of
# https://arxiv.org/pdf/1910.01998.pdf
function Polynomials.roots(p::Bernstein{𝐍,T}) where {𝐍,T}
    R = eltype(one(T) / 1)

    as = zeros(R, 𝐍 + 1)
    for (i, pᵢ) in pairs(coeffs(p))
        as[i] = pᵢ
    end

    Ap = diagm(-1 => ones(R, 𝐍 - 1))
    for j in 1:𝐍
        Ap[1, j] = -as[end - j]
    end

    Bp = diagm(-1 => ones(R, 𝐍 - 1))
    for j in 1:𝐍
        Bp[1, j] = -as[end - j]
    end

    for j in 2:𝐍
        Bp[j, j] = j / (𝐍 - j + 1)
    end
    Bp[1, 1] += as[end] / 𝐍

    eigvals(Ap, Bp)
end

# http://www.cecm.sfu.ca/personal/pborwein/MITACS/papers/LinMatPol.pdf
# Sec 4.1

function comrade_matrix(p::P) where {𝐍,T,P<:Bernstein{𝐍,T}}
    C₀, C₁ = comrade_pencil(p)
    C₀ * inv(C₁)
end

function comrade_pencil(p::P) where {𝐍,T,P<:Bernstein{𝐍,T}}
    n = 𝐍
    a, b = 0, 1
    C₀ = diagm(
        0 => a / (b - a) .* ((n:-1:1) ./ (1:n)) .* ones(T, n),
        -1 => b / (b - a) .* ones(T, n - 1),
    )
    C₁ = diagm(
        0 => 1 / (b - a) .* ((n:-1:1) ./ (1:n)) .* ones(T, n),
        -1 => 1 / (b - a) .* ones(T, n - 1),
    )

    C₀[end, end] *= p[end]
    C₁[end, end] *= p[end]

    for i in 0:(n - 1)
        C₀[1 + i, end] -= b / (b - a) * p[i]
        C₁[1 + i, end] -= 1 / (b - a) * p[i]
    end
    #    C₁[1,end] = -(n-1)/(b-a) * p[0]

    C₀, C₁
end

# block version
function comrade_pencil(p::P) where {𝐍,T,M<:AbstractMatrix{T},P<:Bernstein{𝐍,M}}
    m, m′ = size(p[0])
    @assert m == m′

    Iₘ = I(m)
    n = 𝐍

    a, b = 0, 1 # could generalize with Bernsteing
    R = eltype(one(T) / 1)
    C₀ = zeros(R, n * m, n * m)
    C₁ = zeros(R, n * m, n * m)
    Δ = 1:m

    for i in 1:(n - 1)
        C₀[(i - 1) * m .+ Δ, (i - 1) * m .+ Δ] .= (n + 1 - i) / i * a / (b - a) .* Iₘ
        C₀[i * m .+ Δ, (i - 1) * m .+ Δ] .= b / (b - a) .* Iₘ

        C₁[(i - 1) * m .+ Δ, (i - 1) * m .+ Δ] .= (n + 1 - i) / i * 1 / (b - a) .* Iₘ
        C₁[i * m .+ Δ, (i - 1) * m .+ Δ] .= 1 / (b - a) .* Iₘ
    end

    C₀[(n - 1) * m .+ Δ, (n - 1) * m .+ Δ] .= 1 / n * a / (b - a) .* p[n]
    C₁[(n - 1) * m .+ Δ, (n - 1) * m .+ Δ] .= 1 / n * 1 / (b - a) .* p[n]

    for i in 0:(n - 1)
        C₀[i * m .+ Δ, (n - 1) * m .+ Δ] .-= b / (b - a) .* p[i]
        C₁[i * m .+ Δ, (n - 1) * m .+ Δ] .-= 1 / (b - a) .* p[i]
    end

    #    C₁[Δ, (n-1)*Δ] .= -(n-1) * 1 / (b-a) .* p[0]

    C₀, C₁
end
