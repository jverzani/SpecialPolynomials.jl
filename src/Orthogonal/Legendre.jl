## Legendre Polynomials = Gegenbauer{1/2}
@register0 Legendre AbstractCCOP0
export Legendre


"""
    Legendre{T}

Implements the [Legendre](https://en.wikipedia.org/wiki/Legendre_polynomials) polynomials. These have weight function `w(x) = 1` over the domain `[-1,1]`.

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> p = Legendre([1,2,3])
Legendre(1⋅P_0(x) + 2⋅P_1(x) + 3⋅P_2(x))

julia> convert(Polynomial, p)
Polynomials.Polynomial(-1//2 + 2//1*x + 9//2*x^2)

julia> p2m, p2m1 = basis.(Legendre, (8,9)) # evaluation P_{2m+k}(-1) =  (-1)^k
(Legendre(1⋅P_8(x)), Legendre(1⋅P_9(x)))

julia> p2m(-1) == 1
true

julia> p2m1(-1) == -1
true

julia> n = 5  # verify  Rodrigues' formula 
5

julia> x = Polynomial(:x)
Polynomials.Polynomial(x)

julia> derivative((x^2-1)^n, n) - 2^n *  factorial(n) * basis(Legendre, n)
Polynomials.Polynomial(0//1)

julia> p4, p5  =  basis.(Legendre, (4,5)) # verify  orthogonality  of  P₄,P₅
(Legendre(1⋅P_4(x)), Legendre(1⋅P_5(x)))

julia> SpecialPolynomials.innerproduct(Legendre, p4,  p5)
0.0
```
"""
Legendre

abcde(::Type{<:Legendre})  = NamedTuple{(:a,:b,:c,:d,:e)}((-1,0,1,-2,0))

basis_symbol(::Type{<:Legendre})  = "P"
Polynomials.domain(::Type{<:Legendre}) = Polynomials.Interval(-1, 1)


kn(P::Type{<:Legendre}, n::Int)  = kn(Gegenbauer{1/2, eltype(P)}, n)
k1k0(P::Type{<:Legendre}, n::Int) = k1k0(Gegenbauer{1/2, eltype(P)}, n)
k1k_1(P::Type{<:Legendre}, n::Int) = k1k_1(Gegenbauer{1/2, eltype(P)}, n)

norm2(::Type{<:Legendre}, n) = 2/(2n+1)
weight_function(::Type{<:Legendre})  = x -> one(x)
generating_function(::Type{<:Legendre}) = (t, x)  -> 1/sqrt(1 - 2x*t +t^2)

# gauss nodes
function gauss_nodes_weights(P::Type{<:Legendre}, n)
    xs,  ws =  glaser_liu_rokhlin_gauss_nodes(basis(MonicLegendre,n))
    λ = kn(P,n)^2
    xs, ws/λ
end


# overrides
B̃n(P::Type{<:Legendre}, ::Val{0}) = B̃n(Gegenbauer{1/2, eltype(P)}, Val(0))
C̃n(P::Type{<:Legendre}, ::Val{0}) = zero(eltype(P))
b̂̃n(P::Type{<:Legendre}, ::Val{0}) = b̂̃n(Gegenbauer{1/2, eltype(P)}, Val(0))
ĉ̃n(P::Type{<:Legendre}, ::Val{0}) = ĉ̃n(Gegenbauer{1/2, eltype(P)}, Val(0))

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

##
## --------------------------------------------------
##

# Monic
@register0 MonicLegendre AbstractCCOP0
export MonicLegendre
ϟ(::Type{<:MonicLegendre}) = Legendre
ϟ(::Type{<:MonicLegendre{T}}) where {T} = Legendre{T}
@register_monic(MonicLegendre)

## fast  gauss nodes
pqr_symmetry(::Type{<:MonicLegendre}) = true
pqr_weight(P::Type{<:MonicLegendre}, n, x, dπx) = (one(eltype(P)) * 2)/(1-x^2)/dπx^2

##
## --------------------------------------------------
##

@register0 ShiftedLegendre AbstractCCOP0
"""
    ShiftedLegendre

Type for the shifted Legendre polynomials: `Pˢᵢ(x) =  Pᵢ(2x-1)` for `x ∈ [0,1]`.
"""
ShiftedLegendre

ϟ(::Type{<:ShiftedLegendre})=Legendre
ϟ(::Type{<:ShiftedLegendre{T}}) where {T} = Legendre{T}
export ShiftedLegendre
@register_shifted(ShiftedLegendre, 2, -1)

B̃n(P::Type{<:ShiftedLegendre}, ::Val{0}) = -one(eltype(P))/2

##
## --------------------------------------------------
##

@register0 MonicShiftedLegendre AbstractCCOP0
export MonicShiftedLegendre
ϟ(::Type{<:MonicShiftedLegendre}) = ShiftedLegendre
ϟ(::Type{<:MonicShiftedLegendre{T}}) where {T} = ShiftedLegendre{T}
@register_monic(MonicShiftedLegendre)


