## ShiftedJacobi Polynomials
struct ShiftedJacobiBasis{α,β} <: AbstractCCOPBasis end

"""
    ShiftedJacobi{α,  β, T}

Implements the [ShiftedJacobi](https://arxiv.org/abs/2004.09801) orthogonal polynomials. These have weight function `w(x) = (1-x)^α ⋅ x^β` over the domain `[0,1]`. The parameters are  specified to the constructor.

## Example

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> b = ShiftedJacobi{1,1}([0,0,1])
ShiftedJacobi{1,1}(1⋅Rᵅᵝ₂(x))

julia> convert(Polynomial, b)
Polynomial(3.0 - 15.0*x + 15.0*x^2)
```

Note: These are `Jᵢᵅᵝ(2x-1)`, `J=basis(Jacobi{α,β},i)`.
"""
ShiftedJacobi = MutableDensePolynomial{ShiftedJacobiBasis{α,β}} where {α,β}
export ShiftedJacobi

Polynomials._typealias(::Type{P}) where {α, β, P<:ShiftedJacobi{α, β}} = "ShiftedJacobi{$α,$β}"

basis_symbol(::Type{<:AbstractUnivariatePolynomial{ShiftedJacobiBasis{α,β}}}) where {α,β} = "Rᵅᵝ"
Polynomials.domain(::Type{<:ShiftedJacobiBasis{α,β}}) where {α,β} =
    Polynomials.Interval{β >= 0 ? Open : Closed,α >= 0 ? Open : Closed}(0, 1)

# use ζ₀(n)Rₙ + ζ₁(n)Rₙ₊₁ + ζ₂(n)Rₙ₊₂ = 0
function ζ₀(P::Type{<:ShiftedJacobiBasis{α,β}}, n::Int) where {α,β}
    σ = α + β + 1
    -2(n+α+1)*(n+β+1)*(2n+σ+3)
end
function ζ₁ₓ(P::Type{<:ShiftedJacobiBasis{α,β}}, n::Int) where {α,β}
    σ = α + β + 1
    (2n+σ+2) * (2n+σ+1) * (2n+σ+3) * 2
end
function ζ₁c(P::Type{<:ShiftedJacobiBasis{α,β}}, n::Int) where {α,β}
    σ = α + β + 1
    (2n+σ+2) * ((2n+σ+1) * (2n+σ+3) * (-1) + α^2 - β^2)
end
function ζ₂(P::Type{<:ShiftedJacobiBasis{α,β}}, n::Int) where {α,β}
    σ = α + β + 1
    -2 * (n+2) * (n+σ+1) * (2n+σ+1)
end



# Pₙ₊₁ = (Aₙ*x + Bₙ)Pₙ - CₙPₙ₋₁
function An(P::Type{<:ShiftedJacobiBasis{α,β}}, n::Int) where {α,β}
    n == 0 && return α + β + 2
    -ζ₁ₓ(P,n-1)/ζ₂(P,n-1)
end
function Bn(P::Type{<:ShiftedJacobiBasis{α,β}}, n::Int) where {α,β}
    if iszero(n)
        num,den = (α^2 - β^2 - (α + β)*(α + β + 2)), (2*(α + β))
        iszero(num) && return -1 - β
        return num/den
    end
    -ζ₁c(P,n-1)/ζ₂(P,n-1)
end
function Cn(P::Type{<:ShiftedJacobiBasis{α,β}}, n::Int) where {α,β}
    if n == 0
        num, den = α*β*(α + β + 2), ((α + β)*(α + β + 1))
        if iszero(num)
            if iszero(α)
                d = (β+2)^2
                return 2/d2 + 5β/(2d) + β^2/(2d)
            elseif iszero(β)
                d = (α + 2)^2
                return 2/d2 + 5α/(2d) + α^2/(2d)
            end
        end
        return num/den
    end
    ζ₀(P,n-1)/ζ₂(P,n-1)
end

#kn =  1/2^n ⋅ choose(2n+α+β, n) = (2n+α+β)_n / (2^b  n!)
function leading_term(P::Type{<:ShiftedJacobiBasis{α,β}}, n::Int) where {α,β}
    a = α + β + n + 1
    nn = n
    tot = one(eltype(P)) / 1
    for i in 0:(n - 1)
        tot *= a / (2 * nn)
        a += 1
        nn -= 1
    end
    tot
end

weight_function(::Type{<:ShiftedJacobiBasis{α,β}}) where {α,β} = x -> (1 - x)^α * x^β
#generating_function(::Type{<:ShiftedJacobiBasis{α,β}}) where {α,β}

function classical_hypergeometric(::Type{<:ShiftedJacobiBasis{α,β}}, n, x) where {α,β}
    (α ≤ -1 || β ≤ -1) && throw(ArgumentError("α and β must be > -1"))

    Rₙ =  Pochhammer(α+1,n)/factorial(n)
    Rₙ *= sum(Pochhammer(-n,k)*Pochhammer(n + α + β + 1,k)/factorial(k) / Pochhammer(α+1,k) * (1-x)^k for k in 0:n; init=zero(x))
    Rₙ
end
