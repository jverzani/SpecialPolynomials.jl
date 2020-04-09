"""
    ChebyshevU{T}

Implements the
[Chebyshev](https://en.wikipedia.org/wiki/Chebyshev_polynomials)
polynomials of the second kind. These have weight function 
`w(x) = sqrt(1-x^2)` over the domain `[-1,1]`.

```jldoctest
julia> p = ChebyshevU([1,2,3])
ChebyshevU(1⋅U_0(x) + 2⋅U_1(x) + 3⋅U_2(x))

julia> convert(Polynomial, p)
Polynomial(-2 + 4*x + 12*x^2)

julia> derivative(p)
ChebyshevU(4.0⋅U_0(x) + 12.0⋅U_1(x))

julia> roots(p)
2-element Array{Float64,1}:
 -0.6076252185107651
  0.27429188517743175
```
"""
struct ChebyshevU{T <: Number} <: OrthogonalPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function ChebyshevU{T}(coeffs::AbstractVector{T}, var::Symbol) where {T <: Number}
        length(coeffs) == 0 && return new{T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{T}(coeffs[1:last], var)
    end
end

export ChebyshevU

Polynomials.@register ChebyshevU

basis_symbol(::Type{<:ChebyshevU}) = "U"

weight_function(::Type{<:ChebyshevU})  = x -> sqrt(one(x) - x^2)
generating_function(::Type{<: ChebyshevU}) = (t, x) -> 1 /  (1 - 2t*x + t^2)

Polynomials.domain(::Type{<:ChebyshevU}) = Polynomials.Interval(-1, 1)
Polynomials.variable(::Type{P}, var::Polynomials.SymbolLike=:x) where {P <: ChebyshevU} = P([0, 1], var)/2

An(::Type{<:ChebyshevU}, n) = 2
Bn(::Type{<:ChebyshevU}, n) = 0
Cn(::Type{<:ChebyshevU}, n) = -1

norm2(::Type{ChebyshevU{T}}, n) where {T} = pi/2 * one(T)

## return xs ws
function lagrange_barycentric_nodes_weights(::Type{<:ChebyshevU}, n::Int)
    xs = cos.((0:n-1)*pi/n)
    ws = ones(n)
    ws[1] = ws[end] = 1/2
    for i in  2:2:n
        ws[i] *= -1
    end
    xs, ws
end



(ch::ChebyshevU{T})(x::S) where {T,S} = orthogonal_polyval(ch, x)

# could do directly
function Base.convert(P::Type{<:ChebyshevU}, p::Polynomial)
    q = convert(ChebyshevTT,  p)
    convert(ChebyshevU, q )
end


# direct conversions between U and T types
function Base.convert(P::Type{<:ChebyshevTT}, p::ChebyshevU{T}) where {T}
    # Use
    # Un = 2sum(Tj, j in 1:2:n), n odd
    #    = 2sum(Tj, j in 0:n) - T_0  = 2sum(Tj j in 2:n) + T0, n even
    d = degree(p)
    d == -1 && return ChebysevTT(zeros(T, 1), p.var)
    ts = zeros(T, d + 1)
    q = ChebyshevTT(ts, p.var)
    q[0] = sum(p[i] for i in 0:2:d)
    tot = zero(T)
    for i in d:-2:1
        tot += 2p[i]
        q[i] = tot
    end
    tot = zero(T)
    for i in (d-1):-2:1
        tot += 2p[i]
        q[i] = tot
    end
    return q
end

function Base.convert(P::Type{<:ChebyshevU}, p::ChebyshevTT{T}) where {T <: Number}
    # use  a_i T_i = a_i/2 (U_i - U_i-2)
    d = degree(p)
    R = typeof(one(T)/2)
    d == -1 && return ChebysevU(zeros(R, 1), p.var)
    d == 0 && return ChebyshevU(p[0:0]/1, p.var)
    ts = zeros(R, d + 1)
    q = ChebyshevU(ts, p.var)
    q[d] = p[d]/2
    q[d-1] = p[d-1]/2
    for i in (d-2):-1:1
        q[i] = (p[i]-p[i+2])/2  # uses: a_i T_i = a_i/2 (U_i - U_i-2)
    end
    q[0] = p[0] - p[2]/2        # uses a_0 T_0 = a_0 U_0

    return q
end


function Base.:*(p1::ChebyshevU{T}, p2::ChebyshevU{S}) where {T,S}
    p1.var != p2.var && throw(ArgumentError("Polynomials must have same variable"))

    M, N = length(p1)-1, length(p2) - 1
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

# use T <--> U to take derivatives,  integrate
function Polynomials.derivative(p::ChebyshevU{T}, order::Integer = 1) where {T}
    order < 0 && throw(ArgumentError("Order of derivative must be non-negative"))

    q = convert(ChebyshevTT, p)
    return convert(ChebyshevU, derivative(q, order))

end

function Polynomials.integrate(p::ChebyshevU{T}, C::Number=0) where {T}

    # int Un dx = T_n+1/(n+1)
    R = eltype(one(T)/1)
    n = length(p)
    ts = zeros(R, n+1)
    for i in 1:n
        ts[i+1] = p[i-1]/(i)
    end

    q = ChebyshevTT(ts, p.var)
    q = q - q(0) + R(C)

    return convert(ChebyshevU, q)
end
