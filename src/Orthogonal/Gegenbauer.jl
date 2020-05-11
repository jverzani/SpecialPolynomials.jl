## Gegenbauer Polynomials
@register1 Gegenbauer
export  Gegenbauer

"""
   Gegenbauer{α, T <: Number}

The Gegenbauer polynomials have weight function
`(1-x^2)^(α-1/2)` over the domain `[-1,1]`. The parameter `α` is
specified in the constructor. These are also called the ultra-spherical polynomials.
The Legendre polynomials are the specialization  `Gegenbauer{1/2}`.

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> p =  Gegenbauer{1/2}([1,2,3])
Gegenbauer(1⋅C^(0.5)_0(x) + 2⋅C^(0.5)_1(x) + 3⋅C^(0.5)_2(x))

julia> convert(Polynomial, p)
Polynomials.Polynomial(-0.5 + 2.0*x + 4.5*x^2)
```

"""
Gegenbaeur

basis_symbol(::Type{<:Gegenbauer{α}}) where {α} = "C^($α)"
Polynomials.domain(::Type{<:Gegenbauer{α}}) where {α} = Polynomials.Interval(-1,1)

abcde(::Type{<:Gegenbauer{α}})  where {α} = NamedTuple{(:a,:b,:c,:d,:e)}((-1,0,1,-(2α+1),0))

function kn(::Type{<:Gegenbauer{α}}, n::Int) where {α}
    SpecialPolynomials.Pochhammer_factorial(α, n) * 2^n
end

function k1k0(::Type{<:Gegenbauer{α}}, k, ::Type{S}=Float64) where {S,α}
    iszero(k) && return 2*one(S)*α
    val = one(S)
    val *= 2(α+k)
    val /= k + 1
    val
end
function k1k_1(P::Type{<:Gegenbauer{α}}, k, ::Type{S}=Float64) where {S,α}
    @assert k > 0
    val = one(S)
    if k == 1
        val *= 2α*(α+1)
    else
        val *= 4(α+k)*(α+k-1)
        val /= (k+1)*k
    end
    return val
end


# overrides
Bn(::Type{<:Gegenbauer{α}}, ::Val{0}, ::Type{S}) where  {α,S}  =  zero(S)
b̂n(::Type{<:Gegenbauer{α}}, ::Val{0}, ::Type{S}) where {α,S} = one(S) *  4/α/(α+2) 
ĉn(::Type{<:Gegenbauer{α}}, ::Val{0}, ::Type{S}) where {α,S} = one(S) *  4/α/(α+2) 


## multiplication
⊗(p::Gegenbauer{α}, q::Gegenbauer{α}) where {α} = linearization_product(p, q)
function Base.iterate(o::Linearization{P,R}, state =  nothing) where {α, P <: Gegenbauer{α}, R}
    
    l, m, n = o.l, o.m, o.n
    l  > m + n && return nothing

    if state == nothing

        j = 0
        p = min(l, n)
        q = l - p
        q > m &&  return  nothing

        val = linearization_α(P, R,  j, l, p, q) 

    else

        # we work with l + 2j as there is a parity consideration
        j, p, q, val  = state

        if p == 0 || p <= j || q == m || q < j-1
            j += 1
            p = max(j, min(l+2j, n))
            q = l + 2j - p  # p+q = l + 2j, so l > min(p,q) and l <=  p+q is guaranteed
            while q < j     # need p >= j and q >= j
                p-=1
                q+=1
            end
            q > m && return nothing
            val = linearization_α(P, R, j, l, p, q)
        else
            p -= 1
            q += 1
            ## could  do λ * val here.
            val = linearization_α(P, R, j, l, p, q)
        end

    end

    return (p,q,val), (j, p, q, val)
end


#  https://arxiv.org/pdf/1901.01648.pdf equation (74)
function linearization_α(::Type{P}, R, k, u,n,m) where {α, P <: Gegenbauer{α}}

    val =  one(R)
    val *= (n+m-2k+α) / (n+m-k+α)
    val *= gamma(1 + n+m-2k)
    val /= gamma(1+k)*gamma(1+n-k)*gamma(1+m-k)
    val *= Pochhammer(α, k) * Pochhammer(α, n-k) * Pochhammer(α, m-k) * Pochhammer(2α, n+m-k)
    val /= Pochhammer(α, n+m-k) * Pochhammer(2α, n+m-2k)

    val
    
end

