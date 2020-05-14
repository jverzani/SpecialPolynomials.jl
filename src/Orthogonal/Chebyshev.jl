# Chebyshev Polynomials of first and second kind
@register0 Chebyshev
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
Chebyshev(1⋅T_0(x) + 3⋅T_2(x) + 4⋅T_3(x))

julia> Chebyshev([1, 2, 3, 0], :s)
Chebyshev(1⋅T_0(s) + 2⋅T_1(s) + 3⋅T_2(s))

julia> one(Chebyshev)
Chebyshev(1.0⋅T_0(x))
```

!!! note
    This is copied from the `ChebyshevT` example from the `Polynomials` package by Miles Lucas.


!!! note
The sample chapter available onlinie of [Numerical Methods for Special Functions" by Amparo Gil, Javier Segura, and Nico Temme](https://archive.siam.org/books/ot99/OT99SampleChapter.pdf) gives a very nice overview of these polynomials.
"""
Chebyshev



basis_symbol(::Type{<:Chebyshev})  = "T"
Polynomials.domain(::Type{<:Chebyshev}) = Polynomials.Interval(-1, 1, false, false)
weight_function(::Type{<: Chebyshev}) = x -> one(x)/sqrt(one(x) - x^2)
generating_function(::Type{<: Chebyshev}) =  (t,x) -> (1-t*x)/(1-2*t*x - t^2)


abcde(::Type{<:Chebyshev}) = NamedTuple{(:a,:b,:c,:d,:e)}((-1, 0, 1, -1, 0))

kn(::Type{<:Chebyshev}, n, ::Type{S}) where  {S} = iszero(n) ? one(S) : (2*one(S))^(n-1)
k1k0(::Type{<:Chebyshev}, n, ::Type{S}) where  {S} = iszero(n) ? one(S) : 2*one(S)
k1k_1(::Type{<:Chebyshev}, n, ::Type{S}) where  {S} = n==1 ? 2*one(S) : 4*one(S)

# directly adding these gives a large speed up in polynomial evaluation
#An(::Type{<:Chebyshev}, n::Int, ::Type{S}=Float64) where {S} = iszero(n) ? one(S) : 2*one(S)
#Bn(::Type{<:Chebyshev}, n::Int, ::Type{S}=Float64) where {S} = zero(S)
#Cn(::Type{<:Chebyshev}, n::Int, ::Type{S}=Float64) where {S} = one(S)

# 
Cn(::Type{<:Chebyshev}, ::Val{1}, ::Type{S}=Float64) where {S} = one(S)
ĉn(::Type{<:Chebyshev}, ::Val{0}, ::Type{S}) where {S} = one(S)/4
ĉn(::Type{<:Chebyshev}, ::Val{1}, ::Type{S}) where {S} = Inf
γn(P::Type{<:Chebyshev}, n::Int, ::Type{S}=Float64) where {S} = n==1 ? one(S)/2 : n*one(S)/4 *  k1k0(Chebyshev,n-1, S)

function ⊗(p1::Chebyshev{T}, p2::Chebyshev{S}) where {T,S}

    isconstant(p1) &&  return p2 * p1[0]
    isconstant(p2) &&  return p1 * p2[0]

    p1.var != p2.var && throw(ArgumentError("Polynomials must have same variable"))
    R = promote_type(T,S)
    z1 = _c_to_z(convert(Vector{R}, p1.coeffs))
    z2 = _c_to_z(convert(Vector{R}, p2.coeffs))
    prod = Polynomials.fastconv(z1, z2)
    ret = Chebyshev(_z_to_c(prod), p1.var)
    return truncate(ret)
    
end

## Defining p' and ∫dp directly speeds things up, and works around an issue with ĉ
function Polynomials.derivative(p::Chebyshev{T,N}, order::Int = 1) where {T, N}
    order < 0 && throw(ArgumentError("Order of derivative must be non-negative"))
    R  = eltype(one(T)/1)
    P = Chebyshev{R}
    order == 0 && return convert(P, p)
    hasnan(p) && return P([NaN], p.var)
    order > length(p) && return zero(P)


    q =  convert(P, copy(p))
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

    dp = P(der, p.var)
    return order > 1 ?  derivative(dp, order - 1) : dp

end

function Polynomials.integrate(p::P, C::S) where {T, P<:Chebyshev{T}, S <: Number}
    R = promote_type(eltype(one(T) / 1), S)
    if hasnan(p) || isnan(C)
        return ⟒{P}([NaN])
    end
    n = length(p)
    if n == 1
        return ⟒{P}{R}([C, p[0]])
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

    return ⟒(P)(a2, p.var)
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

## Fitting... still to take

##
## --------------------------------------------------
##

@register0 ChebyshevU
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
ChebyshevU(1⋅U_0(x) + 2⋅U_1(x) + 3⋅U_2(x))

julia> convert(Polynomial, p)
Polynomials.Polynomial(-2 + 4*x + 12*x^2)

julia> derivative(p)
ChebyshevU(4.0⋅U_0(x) + 12.0⋅U_1(x))

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
Polynomials.domain(::Type{<:ChebyshevU}) = Polynomials.Interval(-1, 1)

abcde(::Type{<:ChebyshevU}) = NamedTuple{(:a,:b,:c,:d,:e)}((-1, 0, 1, -3, 0))

kn(::Type{<:ChebyshevU}, n, ::Type{S}) where  {S} = (2*one(S))^n
k1k0(::Type{<:ChebyshevU}, n, ::Type{S}) where  {S} = 2*one(S)
k1k_1(::Type{<:ChebyshevU}, n, ::Type{S}) where  {S} = 4 * one(S)


# directly adding these gives a large speed up in polynomial evaluation
An(::Type{<:ChebyshevU}, n::Int, ::Type{S}=Float64) where {S} = 2
Bn(::Type{<:ChebyshevU}, n::Int, ::Type{S}=Float64) where {S} = 0
Cn(::Type{<:ChebyshevU}, n::Int, ::Type{S}=Float64) where {S} = 1
# work around cancellation
ĉn(::Type{<:ChebyshevU}, n::Int, ::Type{S}) where {S} = -one(S) /(4n+4) *  k1k0(ChebyshevU,n-1,S)

function ⊗(p1::ChebyshevU{T}, p2::ChebyshevU{S}) where {T,S}

    isconstant(p1) &&  return p2 * p1[0]
    isconstant(p2) &&  return p1 * p2[0]
    p1.var != p2.var && throw(ArgumentError("Polynomials must have same variable"))

    M, N = degree(p1), degree(p2)
    R = promote_type(T, S)
    out = zeros(R, M + N + 1)

    # use Um * Un = sum(U_{m-n+2k} k in 0:n) (wikipedia page)
    for m in 0:degree(p1)
        am = p1[m]
        for n in 0:degree(p2)
            bn = p2[n]
            ambn = am * bn
            for k in max(m,n)-min(m,n):2:m+n
                out[k+1] += ambn
            end
        end
    end

    ChebyshevU(out, p1.var)
end



## σ_P = σ_Q:
Base.convert(::Type{Q}, p::P) where {Q <: Chebyshev, P<: ChebyshevU} = _convert_ccop(Q,p)
Base.convert(::Type{Q}, p::P) where {Q <: ChebyshevU, P<: Chebyshev} = _convert_ccop(Q,p)


## Why does default fail?
# use T <--> U to take derivatives,  integrate
# function Polynomials.derivative(p::ChebyshevU{T}, order::Integer = 1) where {T}
#     order < 0 && throw(ArgumentError("Order of derivative must be non-negative"))
    
#     q = convert(Chebyshev, p)
#     return convert(ChebyshevU, derivative(q, order))

# end

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
