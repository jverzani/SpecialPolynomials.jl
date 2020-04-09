export ChebyshevTT

## Contribued by @mileslucas to Polynomials; temporarily renamed while `ChebyshevT` is in `Polynomials`

"""
    ChebyshevTT{<:Number}(coeffs::AbstractVector, var=:x)

Chebyshev polynomial of the first kind.

Construct a polynomial from its coefficients `a`, lowest order first, optionally in
terms of the given variable `x`. `x` can be a character, symbol, or string.

# Examples

```jldoctest
julia> ChebyshevTT([1, 0, 3, 4])
ChebyshevTT([1, 0, 3, 4])

julia> ChebyshevTT([1, 2, 3, 0], :s)
ChebyshevTT([1, 2, 3])

julia> one(ChebyshevTT)
ChebyshevTT([1.0])
```

!!! note
    The name `ChebyshevTT` will be replaced with `ChebyshevT` once that example is removed from the `Polynomials` package.

"""
struct ChebyshevTT{T <: Number} <: OrthogonalPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function ChebyshevTT{T}(coeffs::AbstractVector{T}, var::Symbol) where {T <: Number}
        length(coeffs) == 0 && return new{T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{T}(coeffs[1:last], var)
    end
end

Polynomials.@register ChebyshevTT

basis_symbol(::Type{<:ChebyshevTT}) = "T"
Polynomials.domain(::Type{<:ChebyshevTT}) = Polynomials.Interval(-1, 1, false, false)
weight_function(::Type{<: ChebyshevTT}) = x -> one(x)/sqrt(one(x) - x^2)
generating_function(::Type{<: ChebyshevTT}) =  (t,x) -> (1-t*x)/(1-2*t*x - t^2)

Polynomials.variable(::Type{P}, var::Polynomials.SymbolLike=:x) where {P <: ChebyshevTT} = P([0, 1], var)

An(::Type{<:ChebyshevTT}, n) = iszero(n) ? 1 : 2
Bn(::Type{<:ChebyshevTT}, n) = 0
Cn(::Type{<:ChebyshevTT}, n) = -1
P0(::Type{<:ChebyshevTT}, x) = one(x)
P1(::Type{<:ChebyshevTT}, x) = x

norm2(P::ChebyshevTT{T}, n)  where {T} = iszero(n) ? pi*one(T)/1 : pi*one(T)/2

function lagrange_barycentric_nodes_weights(::Type{<: ChebyshevTT}, n::Int)
    xs = [cos((2j+1)*pi/(2n+2)) for j in 0:n-1]
    ws = [(-1)^j*sin((2j+1)*pi/(2n+2)) for j in 0:n-1]  # XX one loop
    xs, ws
end

## used discrete cosine  transformation to compute the ck:
## https://archive.siam.org/books/ot99/OT99SampleChapter.pdf
cks(P::Type{ChebyshevTT{T}}, f, n::Int) where  {T} =  [ck(P,f,k,n) for k in 0:n]
function ck(P::Type{ChebyshevTT{T}}, f, k::Int, n::Int) where  {T}
    #n =  k
    tot = zero(T)/1
    an   =  pi/(n+1)
    for j in 0:n
        thetaj = (j+1/2)*an
        tot += f(cos(thetaj)) * cos(k * thetaj)
    end
    (iszero(k) ? 1 : 2) * tot  / (n+1)
end



"""
    (::ChebyshevTT)(x)

Evaluate the Chebyshev polynomial at `x`. If `x` is outside of the domain of [-1, 1], an error will be thrown. The evaluation uses Clenshaw Recursion.

# Examples
```jldoctest
julia> c = ChebyshevTT([2.5, 1.5, 1.0])
ChebyshevTT([2.5, 1.5, 1.0])

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
(ch::ChebyshevTT{T})(x::S) where {T,S} = orthogonal_polyval(ch, x)


function Base.:*(p1::ChebyshevTT{T}, p2::ChebyshevTT{S}) where {T,S}
    p1.var != p2.var && throw(ArgumentError("Polynomials must have same variable"))
    R = promote_type(T,S)
    z1 = _c_to_z(convert(Vector{R}, p1.coeffs))
    z2 = _c_to_z(convert(Vector{R}, p2.coeffs))
    prod = Polynomials.fastconv(z1, z2)
    ret = ChebyshevTT(_z_to_c(prod), p1.var)
    return truncate!(ret)
end

function Polynomials.integrate(p::ChebyshevTT{T}, C::S) where {T,S <: Number}
    R = promote_type(eltype(one(T) / 1), S)
    if hasnan(p) || isnan(C)
        return ChebyshevTT([NaN])
    end
    n = length(p)
    if n == 1
        return ChebyshevTT{R}([C, p[0]])
    end
    a2 = Vector{R}(undef, n + 1)
    a2[1] = zero(R)
    a2[2] = p[0]
    a2[3] = p[1] / 4
    @inbounds for i in 2:n - 1
        a2[i + 2] = p[i] / (2 * (i + 1))
        a2[i] -= p[i] / (2 * (i - 1))
    end
    a2[1] += R(C) - ChebyshevTT(a2)(0)
    return ChebyshevTT(a2, p.var)
end


function Polynomials.derivative(p::ChebyshevTT{T}, order::Integer = 1) where {T}
    order < 0 && throw(ArgumentError("Order of derivative must be non-negative"))
    R  = eltype(one(T)/1)
    order == 0 && return convert(ChebyshevTT{R}, p)
    hasnan(p) && return ChebyshevTT(R[NaN], p.var)
    order > length(p) && return zero(ChebyshevTT{R})


    q =  convert(ChebyshevTT{R}, copy(p))
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

    pp = ChebyshevTT(der, p.var)
    return order > 1 ?  derivative(pp, order - 1) : pp

end

function Polynomials.companion(p::ChebyshevTT{T}) where T
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

function Base.divrem(num::ChebyshevTT{T}, den::ChebyshevTT{S}) where {T,S}
    num.var != den.var && throw(ArgumentError("Polynomials must have same variable"))
    n = length(num) - 1
    m = length(den) - 1

    R = typeof(one(T) / one(S))
    P = ChebyshevTT{R}

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
