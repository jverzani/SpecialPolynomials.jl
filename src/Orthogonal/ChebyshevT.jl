export Chebyshev

## Contribued by @mileslucas to Polynomials; temporarily renamed while `ChebyshevT` is in `Polynomials`

"""
    Chebyshev{<:Number}(coeffs::AbstractVector, var=:x)

Chebyshev polynomial of the first kind.

Construct a polynomial from its coefficients `a`, lowest order first, optionally in
terms of the given variable `x`. `x` can be a character, symbol, or string.

# Examples

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> Chebyshev([1, 0, 3, 4])
Chebyshev(1⋅T_0(x) + 3⋅T_2(x) + 4⋅T_3(x))

julia> Chebyshev([1, 2, 3, 0], :s)
Chebyshev(1⋅T_0(s) + 2⋅T_1(s) + 3⋅T_2(s))

julia> one(Chebyshev)
Chebyshev(1.0⋅T_0(x))
```

!!! note
    This is copied from the `ChebyshevT` example from the `Polynomials` package by Miles Lucas.

"""
struct Chebyshev{T <: Number} <: OrthogonalPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function Chebyshev{T}(coeffs::AbstractVector{T}, var::Symbol) where {T <: Number}
        length(coeffs) == 0 && return new{T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{T}(coeffs[1:last], var)
    end
end

Polynomials.@register Chebyshev

basis_symbol(::Type{<:Chebyshev}) = "T"
Polynomials.domain(::Type{<:Chebyshev}) = Polynomials.Interval(-1, 1, false, false)
weight_function(::Type{<: Chebyshev}) = x -> one(x)/sqrt(one(x) - x^2)
generating_function(::Type{<: Chebyshev}) =  (t,x) -> (1-t*x)/(1-2*t*x - t^2)

Polynomials.variable(::Type{P}, var::Polynomials.SymbolLike=:x) where {P <: Chebyshev} = P([0, 1], var)

## [Numerical Methods for Special Functions" by Amparo Gil, Javier Segura, and Nico Temme](https://archive.siam.org/books/ot99/OT99SampleChapter.pdf) is one source, of many  for these.
An(::Type{<:Chebyshev}, n) = iszero(n) ? 1 : 2
Bn(::Type{<:Chebyshev}, n) = 0
Cn(::Type{<:Chebyshev}, n) = -1
P0(::Type{<:Chebyshev}, x) = one(x)
P1(::Type{<:Chebyshev}, x) = x

alpha(::Type{<:Chebyshev}, n) = 0.0
beta(::Type{<:Chebyshev}, n) = n <= 1 ? 1/2 : 1/4

norm2(P::Chebyshev{T}, n)  where {T} = iszero(n) ? pi*one(T)/1 : pi*one(T)/2
leading_term(::Type{<:Chebyshev},n::Int) = 2^(n-1)


"""
    (::Chebyshev)(x)

Evaluate the Chebyshev polynomial at `x`. If `x` is outside of the domain of [-1, 1], an error will be thrown. The evaluation uses Clenshaw Recursion.

# Examples
```jldoctest
julia> using SpecialPolynomials

julia> c = Chebyshev([2.5, 1.5, 1.0])
Chebyshev(2.5⋅T_0(x) + 1.5⋅T_1(x) + 1.0⋅T_2(x))

julia> c(0)
1.5

julia> c.(-1:0.5:1)
5-element Array{Float64,1}:
 2.0
 1.25
 1.5
 2.75
 5.0
```
"""
function (ch::Chebyshev{T})(x::S) where {T,S}
    ## if x ≈ 1 or -1 do something else
    ## As per: p29 of https://archive.siam.org/books/ot99/OT99SampleChapter.pdf
    if x ≈ 1
        n = degree(ch)
        n <= 0 && return ch[0]
        bn, bn1 = ch[n], zero(T)
        dn = bn
        for r = n-1:-1:1
            dn = 2*(x-1)*bn + dn + ch[r]
            bn, bn1 = dn + bn, bn
        end
        return x*bn - bn1 + ch[0]
    elseif x ≈ -1
        n = degree(ch)
        n <= 0 && return ch[0]
        bn, bn1 = ch[n], zero(T)
        dn = bn
        for r = n-1:-1:1
            dn = 2*(x+1)*bn - dn + ch[r]
            bn, bn1 = dn - bn, bn
        end
        return x*bn - bn1 + ch[0]
    end

    orthogonal_polyval(ch, x)

end


# TODO:
# Koepf, Wolfram. (1999). Efficient Computation of Chebyshev
# Polynomials in Computer Algebra. Computer Algebra Systems: A
# Practical Guide. 79-99.  Might have some gain asymptotically doing
# other tricks e.g computing `basis(Chebyshev,n)(x)` for large `n` can
# be done faster through A^(n-1)*[x,1] with A=[2x -1; 1 0]
# Or using divide and conquer, which for really large n is very efficient
#
function orthogonal_basis_polyval(P::Type{<:Chebyshev}, n, x)
    # no check x ∈ domain(P)
    iszero(n) && return one(x)
    n == 1 && return x
    a,r = divrem(n,2)
    r == 0 && return 2 * orthogonal_basis_polyval(P, a, x)^2 - 1
    2 * orthogonal_basis_polyval(P, a, x) * orthogonal_basis_polyval(P, a+1,x) - x
end
    

function Base.:*(p1::Chebyshev{T}, p2::Chebyshev{S}) where {T,S}
    p1.var != p2.var && throw(ArgumentError("Polynomials must have same variable"))
    R = promote_type(T,S)
    z1 = _c_to_z(convert(Vector{R}, p1.coeffs))
    z2 = _c_to_z(convert(Vector{R}, p2.coeffs))
    prod = Polynomials.fastconv(z1, z2)
    ret = Chebyshev(_z_to_c(prod), p1.var)
    return truncate!(ret)
end

function Polynomials.integrate(p::Chebyshev{T}, C::S) where {T,S <: Number}
    R = promote_type(eltype(one(T) / 1), S)
    if hasnan(p) || isnan(C)
        return Chebyshev([NaN])
    end
    n = length(p)
    if n == 1
        return Chebyshev{R}([C, p[0]])
    end
    a2 = Vector{R}(undef, n + 1)
    a2[1] = zero(R)
    a2[2] = p[0]
    a2[3] = p[1] / 4
    @inbounds for i in 2:n - 1
        a2[i + 2] = p[i] / (2 * (i + 1))
        a2[i] -= p[i] / (2 * (i - 1))
    end
    a2[1] += R(C) - Chebyshev(a2)(0)
    return Chebyshev(a2, p.var)
end


function Polynomials.derivative(p::Chebyshev{T}, order::Integer = 1) where {T}
    order < 0 && throw(ArgumentError("Order of derivative must be non-negative"))
    R  = eltype(one(T)/1)
    order == 0 && return convert(Chebyshev{R}, p)
    hasnan(p) && return Chebyshev(R[NaN], p.var)
    order > length(p) && return zero(Chebyshev{R})


    q =  convert(Chebyshev{R}, copy(p))
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

    pp = Chebyshev(der, p.var)
    return order > 1 ?  derivative(pp, order - 1) : pp

end

function Polynomials.companion(p::Chebyshev{T}) where T
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

function Base.divrem(num::Chebyshev{T}, den::Chebyshev{S}) where {T,S}
    num.var != den.var && throw(ArgumentError("Polynomials must have same variable"))
    n = length(num) - 1
    m = length(den) - 1

    R = typeof(one(T) / one(S))
    P = Chebyshev{R}

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
    return P(q_coeff, num.var), P(r_coeff, num.var)
end

#= Fitting =#

# nodes/weights for fitting a polynomial using Lagrange polynomial
# cf. https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
function lagrange_barycentric_nodes_weights(::Type{<: Chebyshev}, n::Int)

    xs = [cospi((j+1/2)/(n+1)) for j in 0:n]
    ws = [(-1)^j*sinpi((j+1/2)/(n+1)) for j in 0:n]  # XX one loop

    xs, ws
end

# noded for integrating against the weight function

function gauss_nodes_weights(::Type{<: Chebyshev}, n::Int)
    xs = cos.(pi/2n * (2*(1:n).-1))
    ws = pi/n * ones(n)
    xs, ws
end
has_fast_gauss_nodes_weights(::Type{<: Chebyshev}) = true

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

    coeffs(chop!(p, atol= thresh))

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
