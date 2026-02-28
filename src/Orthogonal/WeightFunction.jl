## General  weight function
## Wheeler algorithm, https://authors.library.caltech.edu/86868/1/1.4822929.pdf
## algorithm not considered, just text

abstract type AbstractWeightFunction{T,X} <: AbstractContinuousOrthogonalPolynomial{T,X} end
# parent types for `@registerX` constructor
abstract type WeightFunction{T,X} <: AbstractWeightFunction{T,X} end
abstract type WeightFunction1{α,T,X} <: AbstractWeightFunction{T,X} end
export WeightFunction, WeightFunction1

"""
    WeightFunction{T}

A type for orthogonal polynomials relative to some weight
function. The Wheeler or modified Chebyshev algorithm
([Gautschi](https://www.cs.purdue.edu/homes/wxg/Madrid.pdf),
[Press and Teukolsky](https://doi.org/10.1063/1.4822929))
is used to generate the three-term recurrence relation.

If the second order differential equation, `σ⋅p'' + τ⋅p' + λ⋅p=-` is
known, using that to define the polynomial type would be preferred, as
then several additional properties follow for free.

The key computation is the modified moment, `∫πⱼ dw` where `πⱼ` is the
`j`th basis vector for an associated *monic* system, `P`, and `w` is the
weight function.  These values are registered through the
`@register_weight_function(Type, P, w)` macro, as
illustrated in the examples.


#  Example.

Toy example with `ChebyshevU` being derived using the  `Chebyshev` system.

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> const SP = SpecialPolynomials
SpecialPolynomials

julia> SP.@register0 Toy SP.WeightFunction   # register a  Toy  example

julia> SP.@register_weight_function Toy MonicChebyshev SP.weight_function(ChebyshevU)

julia> [SP.Cn.(Toy, 1:5) SP.Cn.(MonicChebyshevU, 1:5)]
5×2 Matrix{Float64}:
 0.25  0.25
 0.25  0.25
 0.25  0.25
 0.25  0.25
 0.25  0.25
```

Elliptic orthogonal polynomials on  `[-1,1]`. Demo 2 of Gautschi.

```jldoctest WeightFunction
julia> using Polynomials, SpecialPolynomials

julia> const SP=SpecialPolynomials
SpecialPolynomials

julia> N,  ω² = 40, 0.999
(40, 0.999)

julia> w(t) = ((1-ω²*t^2)*(1-t^2))^(-1/2)
w (generic function with 1 method)

julia> SP.@register0 WF SP.WeightFunction

julia> SP.@register_weight_function WF MonicChebyshev w

julia> αs, βs = -SP.Bn.(WF, 0:5), SP.Cn.(WF, 0:5);

julia> [αs βs]
6×2 Matrix{Float64}:
 -1.87309e-15  9.68226
 -7.11555e-18  0.793782
 -1.76472e-15  0.119868
 -2.89401e-15  0.22704
 -4.11827e-15  0.241061
 -5.47762e-15  0.245428
```


The main computation involved in this is the modified moment, `νⱼ =
∫πⱼ dw`, computed with `QuadGK.quadgk`. For some examples, this
computation can be completed directly and the `modified_moment` method
may be overloaded for the type. This example is from Press and
Teukolsky, where the modified moments are given through the function
`v(j)` defined below.

```jldoctest
julia> using Polynomials, SpecialPolynomials, SpecialFunctions

julia> const SP=SpecialPolynomials
SpecialPolynomials

julia> w(t) = -log(t)
w (generic function with 1 method)

julia> SP.@register0 WF1 SP.WeightFunction   #  register  type WF1 as a weight function

julia> SP.@register_weight_function WF1  MonicShiftedLegendre w

julia> ν(j) = iszero(j) ? 1 : (-1)^j * gamma(j+1)^2 / (j*(j+1)*gamma(2j+1)) # help  out
ν (generic function with 1 method)

julia> SP.modified_moment(::Type{WF1},  j::Int) = ν(j)

julia> αs, βs = -SP.Bn.(WF1, 0:5), SP.Cn.(WF1, 0:5);

julia> [αs βs]
6×2 Matrix{Float64}:
 -0.75      1.0
 -0.381148  0.423611
 -0.504058  0.172646
 -0.516236  0.203166
 -0.517168  0.222419
 -0.513854  0.239667
```

"""
WeightFunction

## We need w, π
## w a function; π a *monic* system (ismonic(P) == true)
weight_function(::Type{<:AbstractWeightFunction}) =
    throw(ArgumentError("No default method"))
ϟ(P::Type{<:AbstractWeightFunction}) = throw(ArgumentError("No default method"))

Polynomials.domain(W::Type{<:AbstractWeightFunction}) = (domain ∘ ϟ)(W)
Polynomials.basis_symbol(::Type{<:AbstractWeightFunction}) = "W"

## An,Bn,Cn for the monic orthogonal polynomials for the given weight function
## We  have An=1;  α(n) =  -Bn(n); β(n) = Cn(n)
An(::Type{W}, n) where {W<:AbstractWeightFunction} = one(eltype(W))

function Bn(::Type{W}, k::Int) where {W<:AbstractWeightFunction}
    P = ϟ(W)

    if k == 0
        ν₀, ν₁ = modified_moment(W, 0), modified_moment(W, 1)
        out = Bn(P, 0) + ν₁ / ν₀
    else
        out =
            Bn(P, k) - sigma(W, k - 1, k) / sigma(W, k - 1, k - 1) +
            sigma(W, k, k + 1) / sigma(W, k, k)
    end
    -out
end

function Cn(::Type{W}, k::Int) where {W<:AbstractWeightFunction}
    iszero(k) && return innerproduct(W, one, one)
    (sigma(W, k, k) / sigma(W, k - 1, k - 1))
end

@memoize function sigma(W::Type{<:AbstractWeightFunction}, k, l)
    P = ϟ(W)

    T = eltype(W)
    k == -1 && return zero(T) / 1
    k == 0 && return modified_moment(W, l) # ν_l = ∫π_l dw
    k > l && return zero(T) / 1

    # otherwise get through recursion (this uses Gautschi's formula, P&T have apparent error
    sigma(W, k - 1, l + 1) - (Bn(W, k - 1) - Bn(P, l)) * sigma(W, k - 1, l) -
    (Cn(W, k - 1)) * sigma(W, k - 2, l) + Cn(P, l) * sigma(W, k - 1, l - 1)
end

## Compute the modified moment
## v_j =  ∫ π_j ω(x) dx
## (The moment is ∫ x^j ω(x) dx)
## This can  be overridden when an exact calculation is known, e.g.:
function modified_moment(::Type{W}, j::Int) where {W<:AbstractWeightFunction}
    P = ϟ(W)
    @assert ismonic(P)
    modified_moment(basis(P, j), W)
end

@memoize function modified_moment(πⱼ::AbstractOrthogonalPolynomial, W)
    P = ϟ(W)
    ω = weight_function(W)
    πⱼdω = x -> πⱼ(x) * ω(x)
    return ∫(πⱼdω, P, atol=sqrt(eps(float(eltype(W)))))
end

# integrate fdw over region specified by domain(P)
function ∫(fdw, P; kwargs...)
    dom = domain(P)
    a, b = first(dom), last(dom)
    # revert #25; have issue with ∞ away from a,b, but near
    # quadgk would typically handle this, but doesn't get the chance.
    if first(bounds_types(dom)) == Open
        a += eps(float(one(a)))
    end
    if last(bounds_types(dom)) == Open
        b -= eps(float(one(b)))
    end
    out, err = quadgk(fdw, a, b; kwargs...)
    out
end
