struct Legendre{T <: Number} <: OrthogonalPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function Legendre{T}(coeffs::AbstractVector{T}, var::Symbol) where {T <: Number}
        length(coeffs) == 0 && return new{T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{T}(coeffs[1:last], var)
    end
end

export Legendre

Polynomials.@register Legendre

basis_symbol(::Type{<:Legendre}) = "L"

Polynomials.domain(::Type{<:Legendre}) = Polynomials.Interval(-1, 1)
weight_function(::Type{Legendre{T}}) where {T} = x -> one(x)
generating_function(::Type{<:Legendre}) = (t, x)  -> 1/sqrt(1 - 2x*t +t^2)
Polynomials.variable(::Type{P}, var::Polynomials.SymbolLike=:x) where {P <: Legendre} = P([0, 1], var)

# Bonnet's expresssion
An(::Type{Legendre{T}}, n) where {T <: Integer} = iszero(n) ? 1//1 : (2n-1)//n
An(::Type{<:Legendre}, n) = iszero(n) ? 1/1 : (2n-1)/n
Bn(::Type{<:Legendre}, n) = 0
Cn(::Type{Legendre{T}}, n) where {T <: Integer} = -(n-1)//n
Cn(::Type{<:Legendre}, n) = -(n-1)/n

norm2(::Type{<:Legendre}, n) = 2/(2n+1)

(ch::Legendre{T})(x::S) where {T,S} = orthogonal_polyval(ch, x)



# l = n - 2k ; k = (n-l)/2
# https://mathworld.wolfram.com/LegendrePolynomial.html
# compute (2l+1)n! / (2^(n-1)/2 * [1/2(n-l)]! * (l + n + 1)!!)
#      =  (2n+1-4k)n! / (2^k k! (2n+1-2k)!!)
function _legendre_lambda(n, l)
    k = div((n-l),2)
    N = n
    tot = 1/1
    tot *= (2n+1-4k)/1
    for i in 1:k
        tot *= N/(2i)  # n!/(2^k k!)
        N -= 1
    end

    for i in (2n+1-2k):-2:2
        tot *= N/i     # n!/(2n+1 - 2k)!!
        N -= 1
    end
    tot
end

function Base.convert(P::Type{<:Legendre}, p::Polynomial)
    d = degree(p)
    R = eltype(one(eltype(p))/1)
    qs = zeros(R, d+1)
    m = div(d, 2)
    for i in 0:d
        qs[i+1] = sum(p[j] * _legendre_lambda(j, i) for j in i:2:d)
    end
    Legendre(qs, p.var)
end



function Base.:*(p1::Legendre{T}, p2::Legendre{S}) where {T,S}
    p1.var != p2.var && throw(ArgumentError("Polynomials must have same variable"))

    p, q = length(p1)-1, length(p2) - 1
    R = eltype(one(T)  * one(S)/ 1)
    out = zeros(R, p + q + 1)

    ## https://www.cambridge.org/core/services/aop-cambridge-core/content/view/S2040618500035590
    ## P_p * P_q = A0 P_(p+q) + A2  P_(p+q-2)  + A4 P_(p+q-4)
    ## A_2k = (2p+2q-4k+1)/(2p+2q-2k+1) *  l_k l_{p-k} l_{q-k}/l_{p+q-k} wheree
    ## l = (2n)!/(2^n n! n!)
    for m in 0:p
        am = p1[m]
        iszero(am) && continue
        for n in 0:q
            bn = p2[n]
            iszero(bn) && continue
            ambn = am * bn
            for d in p+q:-2:0
                twok = p+q-d
                out[d+1] += ambn * _legendre_A(p,q,twok)
            end
        end
    end

    Legendre(out, p1.var)
end


function _legendre_A(p,q,twok)
    k = div(twok, 2)
    (2p + 2q - 4k+1)/(2p+2q-2k+1) * _legendre_l(k) * _legendre_l(p-k) * _legendre_l(q-k)  / _legendre_l(p+q-k)
end

function _legendre_l(n)
    n < 0 && return 0/1
    out  = 1/1
    for i  in 1:n
        out  *= (2i-1)/i
    end
    out
end


function Polynomials.derivative(p::Legendre{T}, order::Integer = 1) where {T}
    order < 0 && throw(ArgumentError("Order of derivative must be non-negative"))
    order == 0 && return convert(Legendre{T}, p)
    hasnan(p) && return Legendre(T[NaN], p.var)
    order > length(p) && return zero(Legendre{T})

    d = degree(p)
    qs = zeros(T, d)
    for i in 0:d-1
        gamma = 2*i+1
        qs[i+1] = gamma * sum(p[j] for j in (i+1):2:d)
    end

    q = Legendre(qs, p.var)

    if order > 1
        derivative(q, order-1)
    else
        q
    end

end



function Polynomials.integrate(p::Legendre{T}, C::Number=0) where {T}
    # (2n+1)Pn = d/dx (P_(n+1) - P_(n-1))
    R = eltype(one(T)/1)
    d = degree(p)
    qs = zeros(R, d+2)
    qs[1] = -p[0]/3
    for i in 1:(d-1)
        qs[i+1] = p[i-1]/(2*(i-1) + 1) - p[i+1]/(2*(i+1)+1)
    end
    qs[d+1] = p[d-1] / (2(d-1)+1)
    qs[d+2] = p[d] / (2d+1)

    q = Legendre(qs, p.var)
    q = q - q(0) + R(C)

    return q
end
