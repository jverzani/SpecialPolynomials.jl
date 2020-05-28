abstract type AbstractBernstein{T} <: AbstractSpecialPolynomial{T} end

"""

    Bernstein{N, T}

A [Bernstein  polynomial](https://en.wikipedia.org/wiki/Bernstein_polynomial) is a polynomial expressed in terms of
Bernstein basic polynomials. For each degree, `𝐍`, this is a set of `𝐍+1` degree `𝐍` polynomials of the form:
`β_{𝐍,ν} =  (ν choose 𝐍) x^ν  (1-x)^{𝐍-ν}`, `0 ≤ x ≤ 1.`

The `Bernstein{𝐍,T}` type represents a polynomial of degree `𝐍` or
less with a linear combination of the basis vectors using coefficients
of type `T`.

```jldoctest Bernstein
julia> using Polynomials, SpecialPolynomials

julia> p = basis(Bernstein{3},  2)
Bernstein(1.0⋅β₃,₂(x))

julia> convert(Polynomial, p)
Polynomial(3.0*x^2 - 3.0*x^3)
```

!!! Note
    [StaticUnivariatePolynomials](https://github.com/tkoolen/StaticUnivariatePolynomials.jl) Offers a  more  performant version.

"""
struct Bernstein{𝐍, T} <: AbstractBernstein{T}
    coeffs::Vector{T}
    var::Symbol
    function Bernstein{𝐍, T}(coeffs::AbstractVector{T}, var::Polynomials.SymbolLike=:x)  where {𝐍, T}
        N = findlast(!iszero, coeffs)
        N == nothing && return new{𝐍,T}(zeros(T,0), Symbol(var))
        (N > 𝐍 + 1) && throw(ArgumentError("Wrong length for coefficents"))
        return new{𝐍, T}(coeffs[1:N], Symbol(var))
    end
end

export Bernstein
Polynomials.@registerN Bernstein 𝐍
Polynomials.:⟒(::Type{<:Bernstein{𝐍}}) where  {𝐍} = Bernstein{𝐍}

function Polynomials.showterm(io::IO, ::Type{Bernstein{N, T}}, pj::T, var, j, first::Bool, mimetype) where {N, T}
    iszero(pj) && return false
    !first &&  print(io, " ")
    print(io, Polynomials.hasneg(T) && Polynomials.isneg(pj) ? "- " :  (!first ? "+ " : ""))
    print(io, "$(abs.(pj))⋅β")
    Polynomials.unicode_subscript(io, N)
    print(io, ",")
    Polynomials.unicode_subscript(io, j)
    print(io, "($var)")
    return true
end

# add default case with no N specified
function Bernstein(coeffs::AbstractVector{T}, var::Symbol=:x)  where {T}
    l = findlast(!iszero,  coeffs)
    𝐍 = l == nothing ?  -1 :  l - 1
    Bernstein{𝐍, T}(coeffs[1:𝐍+1], var)
end

# promotion between family
function Base.promote(p::P, q::Q) where {𝐍,T,P<:Bernstein{𝐍,T}, 𝐌,S,Q<:Bernstein{𝐌,S}}
    MN = max(𝐌,𝐍)
    R = promote_type(T,S)
    PQ = Bernstein{MN,R}
    convert(PQ,p), convert(PQ,q)
end

# promote up
function Base.convert(P::Type{<:Bernstein{𝐍}}, p::Bernstein{𝐌, S}) where {𝐍, 𝐌, S}
    @assert  𝐍 >=  𝐌
    𝐑 = 𝐍 -  𝐌
    R = eltype(one(promote_type(eltype(P),S))/1)
    cs = zeros(R, 𝐍 +  1)
    for  j  in 0:𝐌
        pⱼ = p[j]
        iszero(pⱼ) &&  continue
        for  i  in j:j+𝐑
            cs[1+i] +=  pⱼ * binomial(𝐌,j) * binomial(𝐑, i-j) /  binomial(𝐍, i)
        end
    end
    Bernstein{𝐍, R}(cs, p.var)
    
end


function Base.convert(P::Type{<:Polynomials.StandardBasisPolynomial}, p::Bernstein{N,T}) where {N, T}
    x = variable(P)
    p(x)
end

function Base.convert(C::Type{<:Bernstein{𝐍, T}}, p::Polynomial) where {𝐍, T}

    @assert degree(p) <= 𝐍

    R = eltype(one(T)/one(Int))
    cs = zeros(R, 𝐍+1)

    for (i, a_nu) in enumerate(coeffs(p))
        k = i - 1
        nk = binomial(𝐍,k)
        for j in k:𝐍
            cs[j+1] += a_nu/nk*binomial(j,k)
        end
    end
    
    Bernstein{𝐍, R}(cs, p.var)
end

function Base.convert(::Type{<:AbstractBernstein}, p::Polynomial{T}) where {T}
    𝐍 = degree(p)
    convert(Bernstein{𝐍, T}, p)
end



Polynomials.domain(::Type{<:AbstractBernstein}) = Polynomials.Interval(0, 1)
Polynomials.degree(p::Bernstein{N,T}) where {N, T} = degree(convert(Polynomial{T}, p))

## need to do one  and variable differently,as coefficients have length N+1
Base.one(P::Type{Bernstein{N,T}}, var::Polynomials.SymbolLike=:x) where {N,T} = Bernstein{N,T}(ones(T, N+1), var)
Base.one(P::Type{<:Bernstein{N}}, var::Polynomials.SymbolLike=:x) where {N} = one(Bernstein{N,eltype(P)}, var)

function Polynomials.variable(P::Type{<:Bernstein{𝐍}}, var::Polynomials.SymbolLike=:x) where {𝐍}
    S = eltype(P)
    𝐍 <= 0 && throw(ArgumentError("Need  𝐍 ≥ 1"))
    Bernstein{𝐍}([1 - (i*one(S))/𝐍 for i in 𝐍:-1:0], var)
end
Polynomials.variable(P::Type{Bernstein},  var::Polynomials.SymbolLike=:x) = variable(Bernstein{1, eltype(P)}, var)

function basis(::Type{Bernstein}, k::Int, _var::Polynomials.SymbolLike=:x; var=_var) 
    zs = zeros(Int, k+1)
    zs[end] = 1
    Bernstein(zs, var)
end

# poly evaluation
(p::AbstractBernstein{T})(x::S) where {T,S} = simple_deCasteljau_polyval(p, x)

# we could  use the compensated one
# from https://doi.org/10.1016/j.camwa.2010.05.021
# Accurate evaluation of a polynomial and its derivative in Bernstein form
# Hao Jiang, Shengguo Li, Lizhi Cheng, and Fang Su
#
# for increased accuracy, but is it really necessary here?
function simple_deCasteljau_polyval(p::Bernstein{𝐍,T}, t::S) where {𝐍,T, S}

    R = promote_type(T, S)
    𝐍 == -1 && return zero(R)
    𝐍 == 0  && return p[0] * one(R)
    n = length(coeffs(p))
    iszero(n) && return zero(R)


    bs = zeros(R, 𝐍+1)  #  pad  out
    for i in eachindex(coeffs(p))
        bs[i] = p[i-1]
    end

    r = 1-t
    for j in 1:𝐍
        for i = 0:𝐍-j
            ii = i + 1
            bs[ii] = bs[ii]*r + bs[ii+1]*t
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


# default in abstract.jl  assumes basis(P,0)=1
function Base.:+(p::Bernstein{𝐍,T}, c::Number) where {𝐍, T}
    u, v =  promote(p, c)
    u+v
end

# default `+` calls promote(p1,p2) then dispatches to this.
function Base.:+(p1::Bernstein{𝐍,T}, p2::Bernstein{𝐍,T}) where {𝐍, T}

    Polynomials.check_same_variable(p1, p2) || throw(ArgumentError("p1 and p2 must have the same symbol"))
    return Bernstein{𝐍, T}([p1[i] + p2[i] for i in 0:𝐍], p1.var)
end

# no promote(p1,p2)  called here
function Base.:*(p::Bernstein{N,T}, q::Bernstein{M,S}) where {N,T, M,S}
    ## use b(n,k) * b(m,j) = choose(n,k)choose(m,j)/choose(n+m,k+j) b(n+m, k+j)

    isconstant(p) && return q * p[0]
    isconstant(q) && return p * q[0]    
    p.var == q.var || throw(ArgumentError("p and q must have the same symbol"))

    
    R = typeof(one(promote_type(T,S))/1)
    x = variable(Polynomial{R})
    return convert(Bernstein, p(x)*q(x))

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
function Polynomials.derivative(p::Bernstein{𝐍, T}, order::Integer = 1) where {𝐍, T}
    cs = coeffs(p)
    𝐍 < -1 &&  return p
    iszero(𝐍) && return zero(p)
    order ==  0 && return p

    cs = zeros(T, 𝐍)
    dp = Bernstein{𝐍-1, T}(𝐍 * diff(coeffs(p)), p.var)
    
    order > 1 ? derivative(dp, order-1) : dp
end



# use ∫b_n,nu  = 1/(n+1)  * ∑_{nu+1}^{n+1} b_n+1, j
function Polynomials.integrate(p::Bernstein{𝐍, T}, C::S) where {𝐍, T,S <: Number}
    R = promote_type(eltype(one(T) / 1), S)
    𝐍𝐍 = 𝐍 + 1
    
    cs = zeros(R, 𝐍𝐍+1)

    @inbounds for ν in eachindex(p)
        for j in (ν+1):𝐍𝐍
            cs[1 + j] += p[ν] / 𝐍𝐍
        end
    end
    
    ∫p = Bernstein{𝐍𝐍, R}(cs, p.var)
    return ∫p - (∫p(0) - C)

end

#  could do p 33 of http://mae.engr.ucdavis.edu/~farouki/bernstein.pdf
function Base.divrem(num::Bernstein{𝐍,T}, den::Bernstein{𝐌,S}) where {𝐍, T, 𝐌, S}
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
