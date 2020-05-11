## Hermite
@register0 Hermite
export Hermite



"""
    Hermite

"""
Hermite


basis_symbol(::Type{<:Hermite}) = "H"
abcde(::Type{<:Hermite})  = NamedTuple{(:a,:b,:c,:d,:e)}((1,0,0,-2,0))

function kn(::Type{<:Hermite}, n::Int)
    2^n
end
function k1k0(::Type{<:Hermite}, k, ::Type{S}=Float64) where {S}
    val = 2*one(S)
    val
end
function k1k_1(P::Type{<:Hermite}, k, ::Type{S}=Float64) where {S}
    @assert k > 0
    k  == 1  && return 2*one(S) #???
    val = 4*one(S)
    return val
end

## Overrides
# Use override here, as we get  0/0 in  default  defn
Bn(::Type{Hermite{T,N}}, n::Int, ::Type{S}) where {T,N,S} = zero(S)
Cn(::Type{Hermite{T,N}}, n::Int, ::Type{S}) where {T,N,S} = 2*n*one(S)


b̂n(::Type{<:Hermite}, n::Val{M}, ::Type{S}) where {M,S} = zero(S)
ĉn(::Type{<:Hermite}, n::Val{M}, ::Type{S}) where {M,S} = zero(S)

## https://arxiv.org/pdf/1901.01648.pdf
##  x^n  = n! sum(H_{n-2j}/ (2^j(n-2j)!j!) j = 0:floor(n/2))
Base.convert(P::Type{<:Polynomial}, q::Hermite) = q(variable(P))
function Base.convert(P::Type{<:Hermite}, q::Polynomial)
    d = degree(q)
    R = eltype(one(eltype(p))/1)
    ps = zeros(R, d+1)
    for i in 0:d
        ps[i+1] = sum(q[jj] * _hermite_lambda(jj, j-1) for (j, jj) in enumerate(i:2:d))
    end
    Hermite(qs, q.var)
end

# compute n!/(2^n * (n-2j)! * j!)
function _hermite_lambda(n,j)
    tot = 1/1
    for i in 1:j
        tot  /=  2
        tot *= n
        n -= 1
    end
    for i in 1:j
        tot /= i
        tot *= n
        n -= 1
    end
    tot
end



##
## --------------------------------------------------
##

## ChebyshevHermite
@register0 ChebyshevHermite
export ChebyshevHermite



"""
    ChebyshevHermite

"""
ChebyshevHermite

basis_symbol(::Type{<:ChebyshevHermite}) = "Hₑ"
# https://arxiv.org/pdf/1901.01648.pdf eqn 17
abcde(::Type{<:ChebyshevHermite})  = NamedTuple{(:a,:b,:c,:d,:e)}((0,0,1,-1,0))

kn(::Type{<:ChebyshevHermite}, n::Int) = 1
k1k0(::Type{<:ChebyshevHermite}, k, ::Type{S}=Float64) where {S} = one(S)
k1k_1(P::Type{<:ChebyshevHermite}, k, ::Type{S}=Float64) where {S} =  one(S)



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
