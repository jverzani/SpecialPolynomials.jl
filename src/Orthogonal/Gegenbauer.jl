"""
   Gegenbauer{α, T <: Number}

The Gegenbauer family of orthogonal polynomials has weight function
`(1-x^2)^(α-1/2)` over the domain `[-1,1]`. The parameter `α` is
specified in the constructor. These are also called the ultra-spherical polynomials.

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> p =  Gegenbauer{1/2}([1,2,3])
Gegenbauer(1⋅C^(0.5)_0(x) + 2⋅C^(0.5)_1(x) + 3⋅C^(0.5)_2(x))

julia> convert(Polynomial, p)
Polynomial(-0.5 + 2.0*x + 4.5*x^2)
```

"""
struct Gegenbauer{α, T <: Number} <: OrthogonalPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function Gegenbauer{α, T}(coeffs::AbstractVector{T}, var::Symbol=:x) where {α, T <: Number}
        length(coeffs) == 0 && return new{α, T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{α, T}(coeffs[1:last], var)
    end
end

export Gegenbauer

Polynomials.@register1 Gegenbauer
constructorof(P::Type{<:Gegenbauer{α}})  where {α} = Gegenbauer{α}

Polynomials.domain(::Type{<:Gegenbauer{α}}) where {α} = Polynomials.Interval(-1, 1, α < 1/2, α < 1/2)
basis_symbol(::Type{<:Gegenbauer{α}}) where {α} = "C^($(round(α,digits=2)))"

weight_function(::Type{<:Gegenbauer{α}}) where {α} = x -> (1-x^2)^(α-1/2)
generating_function(::Type{<:Gegenbauer{α}}) where {α} = (t,x) -> begin
    1/(1-2x*t +t^2)^α
end

# 2^n(α)_n /  n!
function leading_coefficient(::Type{<:Gegenbauer{α}}, n) where {α}
    a = α
    nn = n
    tot = one(n)/1
    for i in 0:n-1
        tot *= a * 2 / nn
        a += 1
        nn -= 1
    end
    tot
end


An(::Type{<:Gegenbauer{α}}, n) where {α} = 2(n+α)/(n+1)
Bn(::Type{<:Gegenbauer{α}}, n) where {α} = 0/1
Cn(::Type{<:Gegenbauer{α}}, n) where {α} = - (n - 1 + 2α)/(n+1)

# compute <Jn, Jn>
function norm2(::Type{<:Gegenbauer{α}}, n) where{α}
    pi * 2^(1-2α) * gamma(n + 2α) / (gamma(n+1) * (n+α) * gamma(α)^2)
end

classical_σ(::Type{<:Gegenbauer})  = x -> 1 - x^2
classical_τ(::Type{<:Gegenbauer{α}}) where {α} = x -> -(2α+1)*x
function classical_hypergeometric(::Type{<:Gegenbauer{α}}, n, x) where {α}
    (α > -1/2 && !iszero(α)) || throw(ArgumentError("α > -1/2 and α ≠ 0 is necessary"))

    as = (-n, 2α+n)
    bs = (α + 1/2,)
    return Pochhammer_factorial(2α,n) * pFq(as, bs, (1-x)/2)
          
end
    

(ch::Gegenbauer)(x::S) where {T,S} = orthogonal_polyval(ch, x)

#Base.convert(::Type{<:Gegenbauer{α}}, q::Gegenbauer{α,T}) where {α, T} = q
function Gegenbauer{α,T}(q::Gegenbauer{β,S})  where {α, β, T, S}
    β == α && return Gegenbauer{α, T}(T.(coeffs(q)), q.var)
    connection(Gegenbauer{α,T}, q)
end
Gegenbauer{α}(q::Gegenbauer{β,S})  where {α, β, S} = Gegenbauer{α,S}(q)

# connection between α and β
function Base.iterate(o::Connection{P, Q}, state=nothing) where
    {β, P <: Gegenbauer{β},
     α, Q <: Gegenbauer{α}}
    
    n,k = o.n, o.k
    if state  == nothing
        j = k
        val = connection_α(P,Q,j,k)
    else
        j, val = state
        j += 2  # likely j += 2 if parity condition
        #λ = connection_α(P,Q,j,k)/connection_α(P,Q,j-2,k)
        #val *= λ  # XXX can do ratio?
        val = connection_α(P,Q,j,k)        
    end
    j > n && return nothing
    (j,val), (j, val)
end

# p29 of thesis
# k = n - 2m; m = (n-k)/2
# replace n-2m with m, (so m -> as written
function connection_α(::Type{<:Gegenbauer{β}},
                      ::Type{<:Gegenbauer{α}},
                      n,k) where {α, β}
    α == β && return (n==k ? one(α) : zero(α))
    
    # definitely replace with ratio
    val = gamma(β)/gamma(α)/gamma(α-β)
    val *= (k + β)
    val *= gamma((n-k)/2 + α - β) * gamma((n+k)/2 + α)
    val /= gamma(1+(n-k)/2)
    val /= gamma((n+k)/2 + β + 1)

    val
                              
end


## Inversion
function Base.iterate(o::Connection{P, Q}, state=nothing) where
    {α, P <: Gegenbauer{α},
     Q <: Polynomials.StandardBasisPolynomial}
    
    n,k = o.n, o.k
    if state  == nothing
        j = k 
        j > n  && return nothing
        val = connection_α(P,Q,j,k)
    else
        j, val = state
        j += 2  #
        j > n && return nothing

        m  = (j-k)/2
        λ = j*(j-1)/4/(m)/(α+k+m)
        val *= λ
        #val = connection_α(P,Q,j,k)
    end
    
    j > n && return nothing
    (j,val), (j, val)
end

# p29 of thesis
function connection_α(::Type{<:Gegenbauer{α}},
                      ::Type{<:Polynomials.StandardBasisPolynomial},
                      n,l) where {α}

    m = (n-l)÷2
    val = gamma(1+n)/2^n
    val *= (l + α)
    val /= gamma(1 + m) * Pochhammer(α, l+1+m)

    val
                              
end


## Multiplication

function Base.iterate(o::Linearization{P}, state =  nothing) where {α, P <: Gegenbauer{α}}

    l, m, n = o.l, o.m, o.n
    l  > m + n && return nothing

    if state == nothing

        j = 0
        p = min(l, n)
        q = l - p
        q > m &&  return  nothing

        val = linearization_α(P, j, l, p, q) 

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
            val = linearization_α(P, j, l, p, q)
        else
            p -= 1
            q += 1
            ## could  do λ * val here.
            val = linearization_α(P, j, l, p, q)
        end

    end

    return (p,q,val), (j, p, q, val)
end


#  https://arxiv.org/pdf/1901.01648.pdf equation (74)
function linearization_α(::Type{P}, k, u,n,m) where {α, P <: Gegenbauer{α}}


    val = (n+m-2k+α) / (n+m-k+α)
    val *= gamma(1 + n+m-2k)
    val /= gamma(1+k)*gamma(1+n-k)*gamma(1+m-k)
    val *= Pochhammer(α, k) * Pochhammer(α, n-k) * Pochhammer(α, m-k) * Pochhammer(2α, n+m-k)
    val /= Pochhammer(α, n+m-k) * Pochhammer(2α, n+m-2k)

    val
    
end




# XXX could tidy up by doing higher orders at once
# XXX this returns answer in a different basis.
function Polynomials.derivative(p::Gegenbauer{α, T}, order::Int=1) where {α,T}
    order == 0 && return p
    order < 0 && throw(ArgumentError("order must be ≥ 0"))

    xs = coeffs(p)
    qs = 2 * α * xs[2:end]
    q = Gegenbauer{α+1}(qs, p.var)

    if order > 1
        return derivative(q, order-1)
    else
        return q
    end

end


# integral: https://mathworld.wolfram.com/GegenbauerPolynomial.html
# \int Pn = 1/(2(n+α)) * [P_{n+1} - P_{n-1}]
function Polynomials.integrate(p::Gegenbauer{α, T}, C::Number=0) where {α,T}

    R = eltype(one(T)/1)
    d = degree(p)
    qs = zeros(R, d+2)
    qs[1] = -p[1]/(2*(1+α))

    for i in 1:(d-1)
        qs[i+1] = p[i-1]/(2(i-1+α)) - p[i+1]/(2(i+1+α))
    end
    qs[d+1] = p[d-1] / (2(d-1+α))
    qs[d+2] = p[d] / (2(d+α))

    q = Gegenbauer{α}(qs, p.var)
    q = q - q(0) + R(C)

    return q
end
