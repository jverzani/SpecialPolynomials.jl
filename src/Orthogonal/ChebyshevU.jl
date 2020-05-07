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

(ch::ChebyshevU{T})(x::S) where {T,S} = orthogonal_polyval(ch, x)

# could do directly
function Base.convert(P::Type{<:ChebyshevU}, p::Polynomial)
    q = convert(Chebyshev,  p)
    convert(ChebyshevU, q )
end

# convert Chebyshev (Q) ->  ChebyshevU (P)
# T_n = 1/2(U_n - U_n-2)
function Base.iterate(o::Connection{P, Q}, state=nothing) where {P <: ChebyshevU, Q <: Chebyshev}

    k, n = o.k, o.n

    if state == nothing
        i = k
        i > n && return nothing
        α = k == 0 ? 1 : 1/2
    elseif state  >= k + 2 # terminate
        return nothing
    else
        i = state
        i += 2
        α = -1/2
    end
    return(i, α), i
end

# convert ChebyshevU (Q) into Chebyshev (P)
# α(n,k) = k=0?-1:2
function Base.iterate(o::Connection{P, Q}, state=nothing) where {P <: Chebyshev, Q <: ChebyshevU}

    k, n = o.k, o.n

    if state == nothing
        i = k 
        i > n && return nothing
    elseif state + 1 > n # terminate
        return nothing
    else
        i = state
        i += 2
    end

    α = k == 0 ? 1 : 2
    return(i, α), i
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

    q = convert(Chebyshev, p)
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

    q = Chebyshev(ts, p.var)
    q = q - q(0) + R(C)

    return convert(ChebyshevU, q)
end

#fitting

## return xs ws
function lagrange_barycentric_nodes_weights(::Type{<:ChebyshevU}, n::Int)
    xs = cospi.((0:n)/n)
    ws = ones(n)
    ws[1] = ws[end] = 1/2
    for i in  2:2:n+1
        ws[i] *= -1
    end
    xs, ws
end
