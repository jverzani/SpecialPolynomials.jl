## Bessel
@registerN Bessel AbstractCCOP1 α

"""
    Bessel{α}

Implements the [Bessel](https://dlmf.nist.gov/18.34) polynomials, introduced by [Krall and Frink](https://www.ams.org/journals/tran/1949-065-01/S0002-9947-1949-0028473-1/S0002-9947-1949-0028473-1.pdf) (with `b=2`). The  case `a=2` corresponds to the
[Bessel](https://en.wikipedia.org/wiki/Bessel_polynomials) polynomials of Wikipedia. The Bessel  polynomials are not orthogonal over  a domain of the real  line, rather over an arbitray curve in the complex plane enclosing the  origin.  The weight  function is `ρ(x)=(2πi)^(-1)∑Γ(α)/Γ(α+n-1)(-β/x)^n`,   where `β=2`.

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> 𝐐 = Rational{Int}
Rational{Int64}

julia> x = variable(Polynomial{𝐐})
Polynomials.Polynomial(x)

julia> [basis(Bessel{3//2, 𝐐}, i)(x) for i in 0:5]
6-element Vector{Polynomial{Rational{Int64}, :x}}:
 Polynomials.Polynomial(1//1)
 Polynomials.Polynomial(1//1 + 3//4*x)
 Polynomials.Polynomial(1//1 + 5//2*x + 35//16*x^2)
 Polynomials.Polynomial(1//1 + 21//4*x + 189//16*x^2 + 693//64*x^3)
 Polynomials.Polynomial(1//1 + 9//1*x + 297//8*x^2 + 1287//16*x^3 + 19305//256*x^4)
 Polynomials.Polynomial(1//1 + 55//4*x + 715//8*x^2 + 10725//32*x^3 + 182325//256*x^4 + 692835//1024*x^5)
```

"""
Bessel
export Bessel
basis_symbol(::Type{<:Bessel{α}}) where {α} = "Cᵅ"
Polynomials.domain(::Type{<:Bessel}) = Polynomials.Interval(0, Inf)

abcde(::Type{<:Bessel{α}}) where {α} = NamedTuple{(:a, :b, :c, :d, :e)}((1, 0, 0, α, 2))

function k1k0(::Type{P}, k::Int) where {α,P<:Bessel{α}}
    k < 0 && return zero(eltype(P)) * one(α) / 1
    iszero(k) && return (one(eltype(P)) * α) / 2

    val = one(eltype(P)) * one(α)
    val *= (2k + α) * (2k + α - 1)
    val /= (k + α - 1) * 2
    val
end

norm2(::Type{<:Bessel{α}}, n) where {α} =
    -1^(n + α - 1) * Γ(1 + n) * 2^(α - 1) / (Γ(n + α - 2) * (2n + α - 1))

# From https://www.ams.org/journals/tran/1949-065-01/S0002-9947-1949-0028473-1/S0002-9947-1949-0028473-1.pdf we use
# kn = (n + α - 1)_n / 2^n not (n+α+1)_n/2^n
# Koepf suggests kn = (n + α + 1)_n/2^n
function leading_term(P::Type{<:Bessel{α}}, n::Int) where {α}
    one(eltype(P)) / 2^n * Pochhammer(n + α - 1, n)
end

weight_function(::Type{<:Bessel{α}}) where {α} =
    z -> (2π * im)^(-1) * z^(α - 2) * exp(-2 / z)
generating_function(::Type{<:Bessel}) = (t, x) -> error("XXX")
function classical_hypergeometric(::Type{<:Bessel{α}}, n, x) where {α}
    as = (-n, n + α - 1)
    bs = ()
    pFq(as, bs, -x / 2) #/ kn
end

## Overrides XXX fails wih 1 and 2
function B̃n(P::Type{<:Bessel{α}}, n::Int) where {α}
    val = one(eltype(P))
    (iszero(n) && α == 2) && return val
    val *= 2 * (α - 2)
    val /= (2 * n + α) * (2 * n + α - 2)
    val
end

function C̃n(P::Type{<:Bessel{α}}, n::Int) where {α}
    val = one(eltype(P))*one(α)
    n == 1 && return -(val * 4) / (α^2 * (α + 1))

    val *= -4 * n * (n + α - 2)
    val /= (2 * n + α - 3) * (2 * n + α - 2)^2 * (2 * n + α - 1)
    val
end

#Bn(::Type{<:Bessel{1}}, n::Int, ::Type{S}) where {S} = error("α=1 is not correct")
#B̃n(P::Type{<:Bessel{2}}, ::Val{0}) where {α} =  one(eltype(P))  # 0  otherwise
#iszero(N) #?  one(eltype(P)) : zero(eltype(P)))
#C̃n(P::Type{<:Bessel{α}}, ::Val{1}) where {α} =  -(one(eltype(P))*4)/(α^2*(α + 1))

b̂̃n(::Type{<:Bessel{2}}, n::Int) = (one(eltype(P)) * 2) / (n * (2n + 2))
b̂̃n(::Type{<:Bessel{2}}, ::Val{0}) = one(eltype(P)) * Inf
