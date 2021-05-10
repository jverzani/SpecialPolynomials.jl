
# Chebyshev Polynomials of first and second kind
@register0 Chebyshev AbstractCCOP0
export Chebyshev
"""
   Chebyshev{<:Number}(coeffs::AbstractVector, var=:x)

Chebyshev polynomial of the first kind.

Construct a polynomial from its coefficients `a`, lowest order first, optionally in
terms of the given variable `x`. `x` can be a character, symbol, or string.

# Examples

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> Chebyshev([1, 0, 3, 4])
Chebyshev(1⋅T₀(x) + 3⋅T₂(x) + 4⋅T₃(x))

julia> Chebyshev([1, 2, 3, 0], :s)
Chebyshev(1⋅T₀(s) + 2⋅T₁(s) + 3⋅T₂(s))

julia> one(Chebyshev)
Chebyshev(1.0⋅T₀(x))
```

!!! note
    This is copied from the `ChebyshevT` example from the `Polynomials` package by Miles Lucas.


!!! note
    The sample chapter available online of [Numerical Methods for Special Functions" by Amparo Gil, Javier Segura, and Nico Temme](https://archive.siam.org/books/ot99/OT99SampleChapter.pdf) gives a very nice overview of these polynomials.
"""
Chebyshev



basis_symbol(::Type{<:Chebyshev})  = "T"
Polynomials.domain(::Type{<:Chebyshev}) = Polynomials.Interval{Open, Open}(-1, 1)
weight_function(::Type{<: Chebyshev}) = x -> one(x)/sqrt(one(x) - x^2)
generating_function(::Type{<: Chebyshev}) =  (t,x) -> (1-t*x)/(1-2*t*x - t^2)
function classical_hypergeometric(P::Type{<:Chebyshev}, n, x) where {α}
    as = (-n,n)
    bs = (one(eltype(P))/2, )
    pFq(as, bs, (1-x)/2)
end

function eval_basis(::Type{Chebyshev}, n, x)
    if abs(x) <= 1
        cos(n * acos(x))
    elseif x > 1
        cosh(n*acosh(x))
    else
        (-1)^n * cosh(n*acosh(-x))
    end
end

abcde(::Type{<:Chebyshev}) = NamedTuple{(:a,:b,:c,:d,:e)}((-1, 0, 1, -1, 0))

k0(P::Type{<:Chebyshev}) = one(eltype(P))
k1k0(P::Type{<:Chebyshev}, n::Int)  = iszero(n) ? one(eltype(P)) : 2*one(eltype(P))

norm2(P::Type{<:Chebyshev}, n::Int) = (iszero(n) ? 1 : 1/2) * pi
ω₁₀(P::Type{<:Chebyshev}, n::Int) = iszero(n) ? one(eltype(P))/sqrt(2) : one(eltype(P)) # √(norm2(n+1)/norm2(n)


leading_term(P::Type{<:Chebyshev}, n::Int)    = iszero(n) ? one(eltype(P)) : (2*one(eltype(P)))^(n-1)

# directly adding these gives a large (20x) speed up in polynomial evaluation
An(P::Type{<:Chebyshev}, n::Int) = iszero(n) ? one(eltype(P)) : 2*one(eltype(P))
Bn(P::Type{<:Chebyshev}, n::Int) = zero(eltype(P))
Cn(P::Type{<:Chebyshev}, n::Int) = one(eltype(P))

# 
C̃n(P::Type{<:Chebyshev}, ::Val{1}) = one(eltype(P))/2
ĉ̃n(P::Type{<:Chebyshev}, ::Val{0}) = one(eltype(P))/4
ĉ̃n(P::Type{<:Chebyshev}, ::Val{1}) = Inf
γ̃n(P::Type{<:Chebyshev}, n::Int) = (n==1) ? one(eltype(P))/2 : n*one(eltype(P))/4

function ⊗(::Type{<:Chebyshev}, p1::Chebyshev{T,X}, p2::Chebyshev{S,Y}) where {T,X,S,Y}

    isconstant(p1) &&  return p2 * p1[0]
    isconstant(p2) &&  return p1 * p2[0]
    assert_same_variable(X,Y)

    R = promote_type(T,S)
    z1 = _c_to_z(convert(Vector{R}, p1.coeffs))
    z2 = _c_to_z(convert(Vector{R}, p2.coeffs))
    prod = Polynomials.fastconv(z1, z2)
    ret = Chebyshev(_z_to_c(prod), X)
    return truncate(ret)
    
end

## Defining p' and ∫dp directly speeds things up, and works around an issue with ĉ
function Polynomials.derivative(p::Chebyshev{T,X,N}, order::Int = 1) where {T,X,N}
    order < 0 && throw(ArgumentError("Order of derivative must be non-negative"))
    R  = eltype(one(T)/1)
    P = Chebyshev{R,X}
    order == 0 && return P(coeffs(p))
    hasnan(p) && return P([NaN])
    order > length(p) && return zero(P)

    q =  P(copy(coeffs(p)))
    n = length(p)
    der = Vector{R}(undef, n)

    for j in n:-1:2
        der[j] = 2j * q[j]
        q[j - 2] += j * q[j] / (j - 2)
    end
    if n > 1
        der[2] = 4q[2]
    end
    der[1] = q[1]

    dp = P(der)
    return order > 1 ?  derivative(dp, order - 1) : dp

end

function Polynomials.integrate(p::P, C::S) where {T, X, P<:Chebyshev{T,X}, S <: Number}
    R = promote_type(eltype(one(T) / 1), S)
    if hasnan(p) || isnan(C)
        return ⟒{P}{R,X}([NaN])
    end
    n = length(p)
    if n == 1
        return ⟒{P}{R,X}([C, p[0]])
    end
    a2 = Vector{R}(undef, n + 1)
    a2[1] = zero(R)
    a2[2] = p[0]
    a2[3] = p[1] / 4
    @inbounds for i in 2:n - 1
        a2[i + 2] = p[i] / (2 * (i + 1))
        a2[i] -= p[i] / (2 * (i - 1))
    end
    a2[1] = C - sum(a2[1+i] for i in 2:2:n)

    return ⟒(P){R,X}(a2)
end


function Polynomials.companion(p::Chebyshev{T}) where {T}
    d = length(p) - 1
    d < 1 && throw(ArgumentError("Series must have degree greater than 1"))
    d == 1 && return diagm(0 => [-p[0] / p[1]])
    R = eltype(one(T) / one(T))

    scl = vcat(1.0, fill(R(√0.5), d - 1))

    diag = vcat(√0.5, fill(R(0.5), d - 2))
    comp = diagm(1 => diag,
                 -1 => diag)
    monics = p.coeffs ./ p.coeffs[end]
    comp[:, end] .-= monics[1:d] .* scl ./ scl[end] ./ 2
    return R.(comp)
end

function Base.divrem(num::Chebyshev{T,X}, den::Chebyshev{S,Y}) where {T,X,S,Y}
    X != Y && throw(ArgumentError("Polynomials must have same variable"))
    n = length(num) - 1
    m = length(den) - 1

    R = typeof(one(T) / one(S))
    P = Chebyshev{R,X}

    if n < m
        return zero(P), convert(P, num)
    elseif m == 0
        den[0] ≈ 0 && throw(DivideError())
        return num ./ den[end], zero(P)
    end

    znum = _c_to_z(num.coeffs)
    zden = _c_to_z(den.coeffs)
    quo, rem = _z_division(znum, zden)
    q_coeff = _z_to_c(quo)
    r_coeff = _z_to_c(rem)
    return P(q_coeff), P(r_coeff)
end

##
## --------------------------------------------------
##


# nodes/weights for fitting a polynomial using Lagrange polynomial
# cf. https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
function lagrange_barycentric_nodes_weights(::Type{<: Chebyshev}, n::Int)

    xs = [cospi((j+1/2)/(n+1)) for j in 0:n]
    ws = [(-1)^j*sinpi((j+1/2)/(n+1)) for j in 0:n]  # XX one loop

    xs, ws
end

# # noded for integrating against the weight function
# function gauss_nodes_weights(::Type{<:Chebyshev}, n::Int)
#     xs = cos.(pi/2n * (2*(n:-1:1).-1))
#     ws = pi/n * ones(n)
#     xs, ws
# end

## cf. fastgaussquadrature
#gauss_nodes_weights(p::Type{P}, n) where {P <: Chebyshev} =
#    FastGaussQuadrature.gausschebyshev(n)

##
## fitting coefficients
cks(val::Val{:interpolating},::Type{Chebyshev}, f, n::Int) = cks(val, Chebyshev{Float64}, f, n)
cks(val::Val{:lsq},::Type{Chebyshev}, f, n::Int) = cks(val, Chebyshev{Float64}, f, n)
cks(val::Val{:series},::Type{Chebyshev}, f) = cks(val, Chebyshev{Float64}, f)

# march forward to compute ck; c0 needs division by 1/2, as used
# ck = 2/(n+1) ∑ f(xⱼ) ⋅ T_j(xⱼ)
# The xs are the same as used by `lagrange_barycentric_nodes_weights`, but
# those weights refer to a different basis.
function cks(::Val{:interpolating}, P::Type{Chebyshev{T}}, f, n::Int) where {T}
    R = float(T)
    k = R(0):n
    ks = (k .+ 1/2)/(n+1)
    xks = cospi.(ks)
    fs = f.(xks)
    a = ones(R,n+1)
    b = copy(xks)
    cks = zeros(R, n+1)
    cks[0+1] = dot(a, fs)/2
    cks[1+1] = dot(b, fs)
    for r = 2:n
        for j in 1:n+1
            a[j], b[j] = b[j], muladd(xks[j], 2b[j],- a[j])
        end
        cks[r+1] = dot(b, fs)
    end
    (2/(n+1))*cks
end

## use discrete cosine  transformation to compute the ck.
## (slower than rrecursion, though maybe a bit more accurate, as implemented)
## https://archive.siam.org/books/ot99/OT99SampleChapter.pdf
# dct(P::Type{Chebyshev{T}}, f, n::Int) where  {T} =  [_dct(T,f,k,n) for k in 0:n]
# dct(::Type{Chebyshev}, f, n::Int) = dct(Chebyshev{Float64}, f,n)
# function _dct(T, f, k::Int, n::Int)
#     tot = zero(float(T))
#     for j in 0:n
#         θⱼ = (j+1/2)/(n+1)
#         tot += f(cospi(θⱼ)) * cospi(k * θⱼ)
#     end
#     ck = (2/(n+1)) * tot
#     ck
# end


# march forward to compute ck
# https://archive.siam.org/books/ot99/OT99SampleChapter.pdf
# ck ≈ 2/n ∑'' f(xⱼ) ⋅ T_k(xⱼ)
# xⱼ = cos(jπ/n)
function cks(::Val{:lsq}, P::Type{Chebyshev{T}}, f, n::Int) where {T}
    R = float(T)
    ks = (R(0):n)/n
    xks = cospi.(ks)
    fs = f.(xks)
    fs[1],fs[end] = fs[1]/2, fs[end]/2  # handle end point ''s
  
    a = ones(R,n+1)
    b = copy(xks)
    out = zeros(R, n+1)

    out[0+1] = dot(a, fs)/2
    out[1+1] = dot(b, fs)

    for r = 2:n
        for j in 1:n+1
            a[j], b[j] = b[j], muladd(xks[j], 2b[j],- a[j])
        end
        out[r+1] = dot(b, fs)
    end
    (2/n)*out
end

## Fit to a series
## use a heuristic to identify `n`
# ref: https://archive.siam.org/books/ot99/OT99SampleChapter.pdf
# Chebyshev interpolation can be interpreted as an approximation to Chebyshev series
# (or vice versa), provided that the coefficients decay fast and the discretization is accurate. In
# other words, Chebyshev series can be a good approximation to near minimax approximations
# (Chebyshev), which in turn are close to minimax approximations.
function cks(::Val{:series}, P::Type{Chebyshev{T}}, f) where {T}

    n = 3
    thresh = 2^n * 10 * eps(T)
    
    p = fit(Val(:interpolating), P, f, 2^n)

    while n < 10
        if maximum(abs, coeffs(p)[end-2^(n-1):end]) <= thresh
            break
        end
        p = fit(Val(:interpolating), P, f, 2^n)
        n += 1
        thresh *= 2
    end

    coeffs(chop(p, atol= thresh))

end

# ## Chebyshev interpolation of the second kind
# ## uses -1, 1 and roots of U_{n-1}, 
# ## (slower than recursion)
# ## used discrete cosine  transformation to compute the ck:
# ## https://archive.siam.org/books/ot99/OT99SampleChapter.pdf
# dct1(P::Type{ChebyshevU{T}}, f, n::Int) where  {T} =  [_dct1(T,f,k,n) for k in 0:n]
# dct1(::Type{Chebyshev},f, n::Int) = dct1(Chebyshev{Float64}, f, n)

# function _dct1(T, f, k::Int, n::Int)

#     tot = zero(float(T))
#     # j=0, n separate due to ∑''
#     j=0
#     θⱼ = j/n
#     tot += f(cospi(θⱼ)) * cospi(k*θⱼ)/2
    
#     j = n
#     θⱼ = j/n
#     tot += f(cospi(θⱼ)) * cospi(k*θⱼ)/2
    
#     for j in 1:n-1
#         θⱼ = j/n
#         tot += f(cospi(θⱼ)) * cospi(k*θⱼ)
#     end

#     ck = (2/n) * tot
#     ck
    
# end

##
## --------------------------------------------------
##
@register0 MonicChebyshev AbstractCCOP0
export MonicChebyshev
ϟ(::Type{<:MonicChebyshev}) = Chebyshev
ϟ(::Type{<:MonicChebyshev{T}}) where {T} = Chebyshev{T}
@register_monic(MonicChebyshev)
C̃n(P::Type{<:MonicChebyshev}, ::Val{1}) = one(eltype(P))


##

@register0 OrthonormalChebyshev AbstractCCOP0
export OrthonormalChebyshev
ϟ(::Type{<:OrthonormalChebyshev}) = Chebyshev
ϟ(::Type{<:OrthonormalChebyshev{T}}) where {T} = Chebyshev{T}
@register_orthonormal(OrthonormalChebyshev)

##
##  --------------------------------------------------
##
## Contrib
#=
zseries =#

function _c_to_z(cs::AbstractVector{T}) where {T}
    n = length(cs)
    U = typeof(one(T) / 2)
    zs = zeros(U, 2n - 1)
    zs[n:end] = cs ./ 2
    return zs .+ reverse(zs)
end

function _z_to_c(z::AbstractVector{T}) where {T}
    n = (length(z) + 1) ÷ 2
    cs = z[n:end]
    cs[2:n] *= 2
    return cs
end

function _z_division(z1::AbstractVector{T}, z2::AbstractVector{S}) where {T,S}
    R = eltype(one(T) / one(S))
    length(z1)
    length(z2)
    if length(z2) == 1
        z1 ./= z2
        return z1, zero(R)
    elseif length(z1) < length(z2)
        return zero(R), R.(z1)
    end
    dlen = length(z1) - length(z2)
    scl = z2[1]
    z2 ./= scl
    quo = Vector{R}(undef, dlen + 1)
    i = 1
    j = dlen + 1
    while i < j
        r = z1[i]
        quo[i] = z1[i]
        quo[end - i + 1] = r
        tmp = r .* z2
        z1[i:i + length(z2) - 1] .-= tmp
        z1[j:j + length(z2) - 1] .-= tmp
        i += 1
        j -= 1
    end

    r = z1[i]
    quo[i] = r
    tmp = r * z2
    z1[i:i + length(z2) - 1] .-= tmp
    quo ./= scl
    rem = z1[i + 1:i - 2 + length(z2)]
    return quo, rem
end

##
## --------------------------------------------------
##

@register0 ChebyshevU AbstractCCOP0
export ChebyshevU
"""
    ChebyshevU{T}

Implements the
[Chebyshev](https://en.wikipedia.org/wiki/Chebyshev_polynomials)
polynomials of the second kind. These have weight function 
`w(x) = sqrt(1-x^2)` over the domain `[-1,1]`.

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> p = ChebyshevU([1,2,3])
ChebyshevU(1⋅U₀(x) + 2⋅U₁(x) + 3⋅U₂(x))

julia> convert(Polynomial, p)
Polynomial(-2 + 4*x + 12*x^2)

julia> derivative(p)
ChebyshevU(4.0⋅U₀(x) + 12.0⋅U₁(x))

julia> roots(p)
2-element Array{Float64,1}:
 -0.6076252185107651
  0.27429188517743175
```
"""
ChebyshevU

basis_symbol(::Type{<:ChebyshevU})  = "U"
weight_function(::Type{<:ChebyshevU})  = x -> sqrt(one(x) - x^2)
generating_function(::Type{<: ChebyshevU}) = (t, x) -> 1 /  (1 - 2t*x + t^2)
function classical_hypergeometric(P::Type{<:ChebyshevU}, n, x) where {α}
    as = (-n,n+2)
    bs = ((3one(eltype(P)))/2, )
    (n+1)*pFq(as, bs, (1-x)/2)
end

Polynomials.domain(::Type{<:ChebyshevU}) = Polynomials.Interval(-1, 1)

abcde(::Type{<:ChebyshevU}) = NamedTuple{(:a,:b,:c,:d,:e)}((-1, 0, 1, -3, 0))

kn(P::Type{<:ChebyshevU}, n::Int) = (2 * one(eltype(P)))^n
k1k0(P::Type{<:ChebyshevU}, n::Int)  = 2 * one(eltype(P))
k1k_1(P::Type{<:ChebyshevU}, n::Int)  = 4 * one(eltype(P))

norm2(P::Type{<:ChebyshevU}, n::Int) = pi/2
ω₁₀(P::Type{<:ChebyshevU}, n::Int) = one(eltype(P))


# directly adding these gives a 5x speed up in polynomial evaluation
An(::Type{<:ChebyshevU}, n::Int) = 2
Bn(::Type{<:ChebyshevU}, n::Int) = 0
Cn(::Type{<:ChebyshevU}, n::Int) = 1
# work around cancellation
ĉ̃n(P::Type{<:ChebyshevU}, n::Int)  = -one(eltype(P)) /(4n+4)

function ⊗(::Type{<:ChebyshevU}, p::ChebyshevU{T,X}, q::ChebyshevU{S,Y}) where {T,X,S,Y}
    isconstant(p) &&  return q * p[0]
    isconstant(q) &&  return p * q[0]
    assert_same_variable(X,Y)

    M, N = degree(p), degree(q)
    R = promote_type(T, S)
    out = zeros(R, M + N + 1)

    # use Um * Un = sum(U_{m-n+2k} k in 0:n) (wikipedia page)
    for i in 0:M
        pᵢ = p[i]
        for j in 0:N
            qⱼ = q[j]
            pᵢqⱼ = pᵢ*qⱼ
            for k in max(i,j)-min(i,j):2:i+j
                out[k+1] += pᵢqⱼ
            end
        end
    end

    ChebyshevU{R,X}(out)
end



## σ_P = σ_Q:e
Base.convert(::Type{Q}, p::P) where {Q <: Chebyshev, P<: ChebyshevU} = _convert_cop(Q,p)
Base.convert(::Type{Q}, p::P) where {Q <: ChebyshevU, P<: Chebyshev} = _convert_cop(Q,p)


# function Polynomials.integrate(p::ChebyshevU{T}, C::Number=0) where {T}

#     # int Un dx = T_n+1/(n+1)
#     R = eltype(one(T)/1)
#     n = length(p)
#     ts = zeros(R, n+1)
#     for i in 1:n
#         ts[i+1] = p[i-1]/(i)
#     end

#     q = Chebyshev(ts, p.var)
#     q = q - q(0) + R(C)

#     return convert(ChebyshevU, q)
# end

##
## --------------------------------------------------
##
## fitting

## return xs ws
function lagrange_barycentric_nodes_weights(P::Type{<:ChebyshevU}, n::Int)
    xs = cospi.((0:n)/n)
    ws = ones(eltype(P), n+1)
    ws[1] /= 2
    ws[end] /= 2
    for i in  2:2:n+1
        ws[i] *= -1
    end
    xs, ws
end

@register0 MonicChebyshevU AbstractCCOP0
export MonicChebyshevU
ϟ(::Type{<:MonicChebyshevU}) = ChebyshevU
ϟ(::Type{<:MonicChebyshevU{T}}) where {T} = ChebyshevU{T}
@register_monic(MonicChebyshevU)


@register0 OrthonormalChebyshevU AbstractCCOP0
export OrthonormalChebyshevU
ϟ(::Type{<:OrthonormalChebyshevU}) = ChebyshevU
ϟ(::Type{<:OrthonormalChebyshevU{T}}) where {T} = ChebyshevU{T}
@register_orthonormal(OrthonormalChebyshevU)



