"""
    Jacobi{α,  β, T}

Implements the [Jacobi](https://en.wikipedia.org/wiki/Jacobi_polynomials) polynomials. These have weight function `w(x) = (1-x)^α ⋅ (1+x)^β` over the domain `[-1,1]`. Many orthogonal polynomial types are special cases. The parameters are  specified to the constructors:

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> p = Jacobi{-1/2, -1/2}([0,0,1])
Jacobi(1⋅J^(α, β)_2(x))

julia> convert(Polynomial, p)
Polynomials.Polynomial(-0.375 + 0.75*x^2)

julia> monic(p::Polynomial) = p/p[end];

julia> (monic ∘ convert).(Polynomial, (p, basis(Chebyshev, 2))) # scaled version of each other
(Polynomials.Polynomial(-0.5 + 1.0*x^2), Polynomials.Polynomial(-0.5 + 1.0*x^2))
```

"""
struct Jacobi{α,  β, T <: Number} <: OrthogonalPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function Jacobi{α, β, T}(coeffs::AbstractVector{T}, var::Symbol=:x) where {α, β, T <: Number}
        α <= -1 && throw(ArgumentError("α > -1 is necessary"))
        β <= -1 && throw(ArgumentError("β > -1 is necessary"))
        length(coeffs) == 0 && return new{α, β, T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{α, β, T}(coeffs[1:last], var)
    end
end

export Jacobi

Polynomials.@register2 Jacobi
constructorof(P::Type{<:Jacobi{α,β}})  where {α,β} = Jacobi{α, β}
  
basis_symbol(::Type{<: Jacobi{α, β}}) where {α, β} = "J^(α, β)"

Polynomials.domain(::Type{<:Jacobi{α, β}}) where {α, β} = Polynomials.Interval(-1, 1, β >= 0, α >= 0)
weight_function(::Type{<:Jacobi{α, β}}) where {α, β} = x -> (1-x)^α *  (1+x)^β
generating_function(::Type{<:Jacobi{α, β}}) where {α, β} = (t,x) -> begin
    R = sqrt(1 - 2x*t+t^2)
    2^(α + β) * 1/R  * (1 - t + R)^(-α) * (1 + t + R)^(-β)
end

# Pochhammer(α+β+n+1)_n / (2^n ⋅ n!)
function leading_coefficient(::Type{<:Jacobi{α, β}}, n) where {α,β}
    a = α+β+n+1
    nn = n
    tot = one(n)/1
    for i in 0:n-1
        tot *= a/(2*nn)
        a += 1
        nn -= 1
    end
    tot
end

An(::Type{<:Jacobi{α, β}}, n) where {α, β} = iszero(n) ? (α+β+2)/2  : (2n+α+β+1) * (2n+α+β+2) / (2(n+1)*(n+α+β+1))
Bn(::Type{<:Jacobi{α, β}}, n) where {α, β} = iszero(n) ? (α-β)/2    : (α^2-β^2)  * (2n+α+β+1) / (2(n+1)*(n+α+β+1)*(2n+α+β))
Cn(::Type{<:Jacobi{α, β}}, n) where {α, β} = iszero(n) ? α*β*(α+β+2)/((α+β)*(α+β+1)) : -(n+α)*(n+β)*(2n+α+β+2)/((n+1)*(n+α+β+1)*(2n+α+β))

(ch::Jacobi)(x::S) where {T,S} = orthogonal_polyval(ch, x)

# compute <Jn, Jn>
function norm2(::Type{<:Jacobi{α, β}}, n) where{α, β}
    α > -1 && β > -1 || throw(ArgumentError("α, β > -1 is necessary"))
    2^(α+β+1)/(2n+α+β+1) * (Γ(n+α+1) *  Γ(n + β +1))/(Γ(n+α+β+1)*Γ(n+1))
end

classical_σ(::Type{<:Jacobi})  = x -> 1 - x^2
classical_τ(::Type{<:Jacobi{α, β}}) where {α, β} = x -> (β-α) - (α+β+2)*x
function classical_hypergeometric(::Type{<:Jacobi{α, β}}, n, x) where {α,β}

    (α ≤ -1 || β ≤ -1) && throw(ArgumentError("α and β must be > -1"))
    
    as = (-n, n+α+β+1)
    bs = (α+1,)
    
    Pochhammer_factorial(α+1,n) * pFq(as, bs, (1 - x)/2)
end
    
    
# Conversion from J{α, β} to J{γ,δ}
function Base.iterate(o::Connection{P, Q}, state=nothing) where
    {γ, δ, P <: Jacobi{γ,δ},
     α, β, Q <: Jacobi{α,β}}
    
    n,k = o.n, o.k
    if state  == nothing
        j = k
        val = connection_α(P,Q,j,k)        
    else
        j, val = state
        j += 1  # likely j += 2 if parity condition
        val = connection_α(P,Q,j,k)        # XXX can do ratio
    end
    j > n && return nothing

    (j,val), (j, val)
end

# p29 of thesis
# gets pretty inaccurate quickly
# might be best to write in terms of ratios, then the
# firsrt three terms could be handled better.
function connection_α(::Type{<:Jacobi{γ, δ}},
                      ::Type{<:Jacobi{α, β}},
                      n,m) where {α, β, γ, δ}

    val = one(α) * Pochhammer_factorial(m+α+1,n-m)
    val *= Pochhammer(n+α+β+1,m)
    val /= Pochhammer(m+γ+δ+1,m)
    as = (m-n,n+m+α+β+1,m+γ+1)
    bs = (m+α+1,2m+γ+δ+2)    
    val *= pFq(as, bs, 1)

    
    val
    
end





# a_nx^n + ... + a_0  x^0 = b_n (1-x)^n + ... + b_0 (1-x)^0
# bs = A \ as, where A is from pascal's triangle
@memoize function _pascal(n)
    A = zeros(n+1,n+1)
    A[1,1] = 1
    for j = 2:n+1
        A[1,j]= 1
        A[j,j] = 1
        for i in 2:j-1
            A[i,j] = A[i,j-1] + A[i-1,j-1]
    end
    end
    for i = 2:2:n+1
        A[i,:] .*= -1
    end
    A
end

    
    

# Inversion formula for standard basis -> Jacobi
# inversion formula used is (1-x)^n = 2^n∑ α(n,k) J_k
# we use pascals triangle to express q in the bases (1-x)^n
function Base.convert(::Type{P}, q::Q) where
    {α, β, P <:Jacobi{α, β},
     T, Q <: Polynomials.StandardBasisPolynomial{T}}

    n = degree(q)
    n < 1 && return P(q[0], q.var)
    qq = Q(_pascal(n)\coeffs(q))

    n = degree(q)
    R = eltype(one(T)/1)
    cs = zeros(R, n+1)
    for k in 0:n
        cs[1+k] = sum(qq[i] * val for (i,val) in Connection{P,Q}(n, k))
    end
    
    ⟒(P)(cs, q.var)
end

function Base.iterate(o::Connection{P, Q}, state=nothing) where
    {α, β, P <: Jacobi{α,β},
     Q <: Polynomials.StandardBasisPolynomial}

    n,k = o.n, o.k
    
    if state  == nothing
        j = k
        val = connection_α(P,Q,j,k)
        
    else
        j, val = state
        j += 1  
        λ= one(α) * 2*(-j)/(-j+k) * (α+j) / (α + β + j + k + 1)
        val *= λ
        #val = connection_α(P,Q,j,k)        
    end
    j > n && return nothing
    
    (j,val), (j, val)
end

# https://arxiv.org/pdf/math/9908148.pdf equation (13)
function connection_α(::Type{<:Jacobi{α,β}},
                      ::Type{<:Polynomials.StandardBasisPolynomial},
                      n,k) where {α,β}

    tot = one(α) * one(β)
    tot *= 2^n
    tot *= Pochhammer(α+β+1,k) * (α+β+2k+1)
    tot *= Pochhammer(-n,k) * Pochhammer(α+k+1,n-k)
    tot /= prod(α+β+i for i in 1:n+k+1) # keeps rational Γ(α+β+1)/Γ(α+β+n+k+2)
    
    return tot
end



# https://arxiv.org/pdf/math/9703217.pdf corollary 1
function Polynomials.integrate(p::P, C::S) where
    {α,β,T, P<:Jacobi{α,β,T},
     S <: Number}
    
    R = promote_type(eltype(one(α) * one(β) *one(T) / 1), S)
    if hasnan(p) || isnan(C)
        return ⟒(P)([NaN])
    end
   d = degree(p)
    if d == 0
        return ⟒(P)([C, p[0]])
    end
    
    as = zeros(R, d + 2)
    
    @inbounds for n in 0:d
        pn = p[n]
        
        as[1 + n + 1] += pn * 2 * checked_div(n+α+β+1, (2n+α+β+1)*(2n+α+β+2) )
        as[1 + n]     += pn * 2 * checked_div(α-β, (2n+α+β)*(2n+α+β+2))
        if  n > 0
                as[1 + n - 1] += pn * (-2) *checked_div((n+α)*(n+β), (n+α+β)*(2n+α+β)*(2n+α+β+1))
        end
    end
    
    ∫p = ⟒(P)(as,  p.var)
    ∫p[0] = R(C) - ∫p(0)

    return  ∫p
    
end
