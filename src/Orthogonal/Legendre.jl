abstract type AbstractLegendre{T} <: OrthogonalPolynomial{T} end

"""
    Legendre{T}

Implements the [Legendre](https://en.wikipedia.org/wiki/Legendre_polynomials) polynomials. These have weight function `w(x) = 1` over the domain `[-1,1]`.

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> p = Legendre([1,2,3])
Legendre(1⋅L_0(x) + 2⋅L_1(x) + 3⋅L_2(x))

julia> convert(Polynomial, p)
Polynomial(-1//2 + 2//1*x + 9//2*x^2)

julia> p2m, p2m1 = basis.(Legendre, (8,9)) # evaluation P_{2m+k}(-1) =  (-1)^k
(Legendre(1⋅L_8(x)), Legendre(1⋅L_9(x)))

julia> p2m(-1) == 1
true

julia> p2m1(-1) == -1
true

julia> n = 5  # verify  Rodrigues' formula 
5

julia> x = Polynomial(:x)
Polynomial(x)

julia> derivative((x^2-1)^n, n) - 2^n *  factorial(n) * basis(Legendre, n)
Polynomial(0.0)

julia> p4, p5  =  basis.(Legendre, (4,5)) # verify  orthogonality  of  L4, L5
(Legendre(1⋅L_4(x)), Legendre(1⋅L_5(x)))

julia> SpecialPolynomials.innerproduct(Legendre, p4,  p5)
6.096184133406375e-16
```



"""
struct Legendre{T <: Number} <: AbstractLegendre{T}
    coeffs::Vector{T}
    var::Symbol
    function Legendre{T}(coeffs::AbstractVector{T}, var::Symbol) where {T <: Number}
        length(coeffs) == 0 && return new{T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        #return new{T}(coeffs[1:last], var)
        new{T}(resize!(coeffs,last), var)
    end
end

export Legendre

Polynomials.@register Legendre

basis_symbol(::Type{<:Legendre}) = "L"

Polynomials.domain(::Type{<:AbstractLegendre}) = Polynomials.Interval(-1, 1)
Polynomials.variable(::Type{P}, var::Polynomials.SymbolLike=:x) where {P <: Legendre} = P([0, 1], var)

weight_function(::Type{<:Legendre})  = x -> one(x)
generating_function(::Type{<:Legendre}) = (t, x)  -> 1/sqrt(1 - 2x*t +t^2)


# Bonnet's expresssion
An(::Type{Legendre{T}}, n) where {T <: Integer} = (2n+1)//(n+1)
An(::Type{<:Legendre}, n) = (2n+1)/(n+1)
Bn(::Type{<:Legendre}, n) = 0
Cn(::Type{Legendre{T}}, n) where {T <: Integer} = -n//(n+1)
Cn(::Type{<:Legendre}, n) = -n/(n+1)

norm2(::Type{<:Legendre}, n) = 2/(2n+1)


pqr(p::Legendre) = (x,n) -> (p=1-x^2, q=-2x, r= n*(n+1), dp=-2x, dq=-2,dr=0)
pqr_symmetry(p::Legendre) = true
pqr_weight(p::Legendre, n, x, dπx) = 2/(1-x^2)/dπx^2
gauss_nodes_weights(P::Type{<:Legendre}, n)  = glaser_liu_rokhlin_gauss_nodes(Polynomials.basis(P,n))
has_fast_gauss_nodes_weights(::Type{<: Legendre}) = true

(ch::Legendre{T})(x::S) where {T,S} = orthogonal_polyval(ch, x)


# from https://mathworld.wolfram.com/LegendrePolynomial.html
# we have xⁱ = \sum_l λ(i,l) L_l
# with
# l = n - 2k ; k = (n-l)/2
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

# A connection α(n,k) returns (k, α(k,k)), (k+1, α(k+1, k)), ..., (n, α(n,k))
function Base.iterate(o::Connection{P, Q}, state=nothing) where {P <: Legendre, Q <: Polynomials.StandardBasisPolynomial}

    k, n = o.k, o.n

    if state == nothing
        i = k
        i > n && return nothing
    elseif state + 2 > n # terminate
        return nothing
    else
        i = state
        i += 2
    end

    return(i, _legendre_lambda(i, k)), i
end


@inline function legendre_A(k, p, q)
    λs = legendre_lambda.((k, p-k, q-k, p+q-k))
    (2p + 2q - 4k +1)/(2p + 2q - 2k + 1) * λs[1] * λs[2] * λs[3] /  λs[4]
end
    
@inline legendre_lambda(n) = n <= 1 ? 1/1 : prod((2*i-1)/i for i in 2:n)
@inline legendre_phi(n) = n <= 0 ? 1 : (2n+1)/(n+1)
const _λ=legendre_lambda

    

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
            ls = legendre_lambda.((m,n,m+n))
            Ak = ls[1]*ls[2]/ls[3]
            k = 0
            for d in m+n:-2:0
                out[d+1] += ambn * Ak
                λ = (2k+1)*(m-k)*(n-k)*(2m+2n-4k-3)*(2m+2n-2k+1)/((k+1)*(m+n-k)*(2m-1-2k)*(2n-1-2k)*(2m+2n-4k+1))
                Ak *= λ
                k += 1
            end
        end
    end
    Legendre(out, p1.var)
end


# in https://www.cambridge.org/core/services/aop-cambridge-core/content/view/S2040618500035590
# we get linearization formula for Pn*Pm = ∑ Ak P_{n+m-2k};
# this allows use of linearization_product. However, it is slower than the direct approach above
# why does this allocate so much?
function Base.iterate(o::Linearization{<:Legendre}, state =  nothing)#

    l,m,n = o.l, o.m, o.n
    
    if state == nothing
        # k = 0
    
        l > m + n && return nothing
        
        p = min(l, n)
        q = l - p
        #val = (2p + 2q - 4k +1)/(2p + 2q - 2k + 1) * _λ(k) * _λ(p-k) * _λ(q-k) /  _λ(p+q-k)
        val = legendre_A(0, p, q)

        return (p,q,val), (0,p,q,val)

    else
        
        k, p,q,val =  state

        if p == 0 || q == m || p <= k 
            # bump k
            k += 1
            l + 2k > n+m && return  nothing
            p = min(l+k, n) # not l+2k, as then q < k and val=0
            q = l + 2k - p # p+q = l + 2k
            val = legendre_A(k, p, q)

            return (p,q,val), (k,p,q,val)

        else

            p -= 1
            q += 1

            λ = (p-k+1) / (q-k) * (2*(q-k)-1) / (2*(p-k)+1)
            val *= λ

            return (p,q,val), (k,p,q,val)
        end
    end
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


## utils


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

