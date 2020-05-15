## Hermite
@register0 Hermite
export Hermite



"""
    Hermite

"""
Hermite


basis_symbol(::Type{<:Hermite}) = "H"
Polynomials.domain(::Type{<:Hermite}) = Polynomials.Interval(-Inf, Inf)
weight_function(::Type{<: Hermite})  = x -> exp(-x^2)
generating_function(::Type{<:Hermite}) = (t, x)  -> exp(2*t*x - t^2)
function classical_hypergeometric(::Type{<:Hermite}, n, x)
    as = iseven(n) ? (-n ÷ 2, -(n-1)/2) : (-n/2, -(n-1)÷2)
    bs = ()
    (2x)^n * pFq(as, bs, -1/x^2)
end




abcde(::Type{<:Hermite})  = NamedTuple{(:a,:b,:c,:d,:e)}((1,0,0,-2,0))

function kn(::Type{<:Hermite}, n::Int)
    2^n
end
function k1k0(P::Type{<:Hermite}, k)
    val = 2*one(eltype(P))
    val
end
function k1k_1(P::Type{<:Hermite}, k)
    @assert k > 0
    #k  == 1  && return 2*one(eltype(P)) #???
    val = 4*one(eltype(P))
    return val
end
norm2(::Type{<:Hermite}, n) = sqrt(pi) * 2^n * gamma(n+1)

## Overrides
# Use override here, as we get  0/0 in  default  defn
Bn(P::Type{<:Hermite}, n::Int) = zero(eltype(P))
Cn(P::Type{<:Hermite}, n::Int) = 2*n*one(eltype(P))

# This is the issue
b̂n(::Type{<:Hermite}, n::Int) where {M} = error("Don't call me")#zero(S)
ĉn(::Type{<:Hermite}, n::Int) where {M} = error("Don't call me")#zero(S)

## https://arxiv.org/pdf/1901.01648.pdf. Connection formula (14)
##  x^n  = n! sum(H_{n-2j}/ (2^j(n-2j)!j!) j = 0:floor(n/2))

Base.convert(P::Type{<:Hermite}, q::Polynomial) = connection(P,q)
function Base.iterate(o::Connection{P, Q}, state=nothing) where
    {P<:Hermite,
     Q<:Polynomials.StandardBasisPolynomial}

    k, n = o.k, o.n

    if state == nothing
        i = k
        j = 0
        i > n && return nothing
        val = __hermite_lambda(i,k)
    elseif state[1] + 2 > n # terminate
        return nothing
    else
        j,val1 = state
        #val1 *= (i+1)*(i+2)/4/(j+1/4)
        j += 1
        i = k + 2j
        val = __hermite_lambda(i,k)
    end

    return(i, val), (j, val)
end

function __hermite_lambda(n,k)
    val = gamma(1+n)/2^n
    val /= gamma(1+k)
    val /= gamma(1 + (n-k)/2)
    val
end


# function Base.convert(P::Type{<:Hermite}, q::Polynomial)
#     d = degree(q)
#     R = eltype(one(eltype(1))/1)
#     ps = zeros(R, max(0,d)+1)
#     for i in 0:max(0,d)
#         ps[i+1] = sum(q[jj] * _hermite_lambda(jj, j-1) for (j, jj) in enumerate(i:2:d))
#     end
#     Hermite(ps, q.var)
# end

# # compute n!/(2^n * (n-2j)! * j!)
# function _hermite_lambda(n,j)
#     tot = 1/1
#     for i in 1:j
#         tot  /=  2
#         tot *= n
#         n -= 1
#     end
#     for i in 1:j
#         tot /= i
#         tot *= n
#         n -= 1
#     end
#     tot
# end

function Polynomials.derivative(p::P, order::Integer = 1) where {T, P <: Hermite{T}}

    order < 0 && throw(ArgumentError("Order of derivative must be non-negative"))
    order == 0 && return p
    hasnan(p) && return Hermite(T[NaN], p.var)
    order > length(p) && return zero(P, p.var)

    d = degree(p)
    qs = zeros(T, d+1-order)
    for i in order:d # Hn' =  2n Hn-1
        qs[i+1-order] = prod(2*(1 + i - j) for j in 1:order)  * p[i]
    end

    q = Hermite(qs, p.var)


end



function Polynomials.integrate(p::P, C::Number=0) where {T,P<:Hermite{T}}
    # int H_n = 1/(n+1) H_{n+1}
    R = eltype(one(T)/1)
    d = degree(p)
    qs = zeros(R, d+2)
    q = ⟒(P)(qs, p.var)

    for i in 0:d
        q[i+1] = p[i]/(2(i+1))
    end

    q = q - q(0) + R(C)

    return q
end


##
## --------------------------------------------------
##
@register0 dChebyshevHermite
export dChebyshevHermite



## ChebyshevHermite
@register0 ChebyshevHermite
export ChebyshevHermite
abcde(::Type{<:dChebyshevHermite})  = NamedTuple{(:a,:b,:c,:d,:e)}((0,0,1,-1,0))

kn(P::Type{<:dChebyshevHermite}, n::Int) =  n+1
k1k0(P::Type{<:dChebyshevHermite}, n)  = (n+2)/(n+1)
k1k_1(P::Type{<:dChebyshevHermite}, n) =  (n+2)/(n)



"""
    ChebyshevHermite

"""
ChebyshevHermite

basis_symbol(::Type{<:ChebyshevHermite}) = "Hₑ"
Polynomials.domain(::Type{<:ChebyshevHermite}) = Polynomials.Interval(-Inf, Inf)
weight_function(::Type{ChebyshevHermite{T}}) where {T} = x -> exp(-x^2/2)
generating_function(::Type{<:ChebyshevHermite}) = (t, x)  -> exp(t*x - t^2/2)
function classical_hypergeometric(::Type{<:ChebyshevHermite}, n, x)
    as = iseven(n) ? (-n ÷ 2, -(n-1)/2) : (-n/2, -(n-1)÷2)
    bs = ()
    (x)^n * pFq(as, bs, -1/x^2)
end

# https://arxiv.org/pdf/1901.01648.pdf eqn 17
abcde(::Type{<:ChebyshevHermite})  = NamedTuple{(:a,:b,:c,:d,:e)}((0,0,1,-1,0))

kn(P::Type{<:ChebyshevHermite}, n::Int) =  one(eltype(P))
k1k0(P::Type{<:ChebyshevHermite}, k)  = one(eltype(P))
k1k_1(P::Type{<:ChebyshevHermite}, k) =  one(eltype(P))



##
## --------------------------------------------------
##

##  Multiplication
## TODO: work  on  types
AbstractHermite = Union{Hermite, ChebyshevHermite}
⊗(p::Hermite, q::Hermite) = linearization_product(p, q)
⊗(p::ChebyshevHermite, q::ChebyshevHermite) = linearization_product(p, q)

function Base.iterate(o::Linearization{P,R}, state =  nothing) where {P <: AbstractHermite,R}

    l, m, n = o.l, o.m, o.n
    l  > m + n && return nothing
    
    if state == nothing

        #k =  0
        j = 0
        p = min(l, n)
        q = l - p
        val = linearization_α(P, R, j, l, p, q) 

    else
        # we work with l + 2j as there is a parity consideration
        j, p, q, val  = state
        s = l + j # l +2j + p + q = l +2j + l = 2l + 2j, so s=l+j
        if p == 0 ||  q == m || q  > s
            j += 1
            l + 2j > n+m && return nothing #  l + 2j  is too  big now
            p = min(l + j, n)  # p <= s
            q = l + 2j - p
            q > s + 1 && return nothing
            val = linearization_α(P, R, j, l, p, q)
        else
            p -= 1
            q += 1

            λ = q/(p+1)*(s-(q-1))/(s-p)            
            val *= λ
        end

    end
                
    return (p,q,val), (j, p, q, val)
end


#  https://arxiv.org/pdf/1901.01648.pdf equation (74)
function linearization_α(P::Type{<:AbstractHermite}, R, j, l, n, m)
    s = l + j # pass in j to avoid s = divrem(l + m + n,  2) call
    val = one(R)
    val *= gamma(1+m)*gamma(1+n)/gamma(1+s-m)/gamma(1+s-n)/gamma(1+s-l)
    val *= linearization_λ(P, l, m,n)
end
linearization_λ(::Type{<:Hermite}, l, m, n) = 2.0^((m+n-l)/2)
linearization_λ(::Type{<:ChebyshevHermite}, l, m, n) = 1    
