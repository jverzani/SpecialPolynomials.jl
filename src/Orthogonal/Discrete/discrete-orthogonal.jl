##################################################

const ∑ = sum

function innerproduct(P::Type{<:AbstractDiscreteOrthogonalPolynomial}, f, g)
    dom = domain(P)
    fn = x -> f(x) * g(x) * weight_function(P)(x)
    a, b = first(dom), last(dom)
    if !isinf(a) && !isinf(b)
        return ∑(fn(x) for x in first(dom):last(dom))
    else
        ## what to do if infinite
    end
end

##
## --------------------------------------------------
##

abstract type AbstractDiscreteWeightFunction{T,X} <:
              AbstractDiscreteOrthogonalPolynomial{T,X} end
abstract type DiscreteWeightFunction{T,X} <: AbstractDiscreteWeightFunction{T,X} end
export DiscreteWeightFunction

"""
    DiscreteWeightFunction

For a discrete measure, `dλ = ∑ wᵢ δ(x - xᵢ)`, specified through two
vectors, `xs` and `ws`, a collection of monic orthogonal polynomials is
produced through Darboux's formula for `α_n` and `β_n` using the
3-term recurrence defined by `π_{n+1} = (x-α_n)⋅π_n - β_n⋅π_{n-1}` (`An=1`, `Bn=-α_n`, `Cn=β_n`)
and the discrete Stieltjes method [Guatschi §3.1](https://www.cs.purdue.edu/homes/wxg/Madrid.pdf).

# Example

Discrete Chebyshev by its weight function (uniform  on 0,1,…,N-1)

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> const SP = SpecialPolynomials;

julia> N = 9
9

julia> xs, ws = collect(0:N-1), ones(N);   # w(x) = ∑ wⱼ⋅δ(x-xⱼ)

julia> SP.@register0 DWF DiscreteWeightFunction

julia> SP.@register_discrete_weight_function(DWF, xs, ws)

julia> [SP.Bn.(DWF, 0:N-1) SP.Cn.(DWF, 0:N-1)]
9×2 Matrix{Float64}:
 -4.0  9.0
 -4.0  6.66667
 -4.0  5.13333
 -4.0  4.62857
 -4.0  4.12698
 -4.0  3.53535
 -4.0  2.83217
 -4.0  2.01026
 -4.0  1.06667

julia> i,j = 3,4; ## check  that ∫pᵢpⱼdw  = 0    for i,j=3,4

julia> sum(basis(DWF,i)(x) *  basis(DWF,j)(x) * w for  (x,w) in zip(xs, ws))
5.684341886080802e-14

julia> ## Gogin, Hirvensalo (https://doi.org/10.1007/s10958-017-3410-8) characterization
       D(k,N,x) =  sum((-1)^l * binomial(k+l,k) * binomial(N-l,k-l) *  SP.generalized_binomial(x,l) for l in 0:k)
D (generic function with 1 method)

julia> x = variable()
Polynomial(x)

julia> ps,qs = [D(k,N-1,x)  for  k in 0:N-1], [basis(DWF, k)(x) for k  in 0:N-1];
ERROR: MethodError: no method matching eval_cop(::Type{DWF{Float64, :x}}, ::Vector{Float64}, ::Polynomial{Int64, :x})

Closest candidates are:
  eval_cop(!Matched::Type{<:SpecialPolynomials.AbstractCOP{T, X}}, ::Any, ::S) where {T, X, S}
   @ SpecialPolynomials ~/julia/SpecialPolynomials/src/Orthogonal/cop.jl:150

Stacktrace:
 [1] evalpoly(x::Polynomial{Int64, :x}, p::DWF{Float64, :x})
   @ Main ~/julia/SpecialPolynomials/src/Orthogonal/abstract.jl:47
 [2] polynomial_composition(p::DWF{Float64, :x}, q::Polynomial{Int64, :x})
   @ Polynomials ~/julia/Polynomials/src/common.jl:1017
 [3] (::DWF{Float64, :x})(x::Polynomial{Int64, :x})
   @ Main ~/julia/Polynomials/src/abstract.jl:124
 [4] (::var"#6#8")(k::Int64)
   @ Main ./none:0
 [5] iterate(::Base.Generator{Base.Iterators.Enumerate{Vector{Any}}, RecipesBase.var"#sub#4"})
   @ Base ./generator.jl:47 [inlined]
 [6] collect(itr::Base.Generator{UnitRange{Int64}, var"#6#8"})
   @ Base ./array.jl:834
 [7] top-level scope
   @ none:1

julia> all(qs .* [p[end] for p  in ps] .≈ ps)
ERROR: UndefVarError: `qs` not defined
Stacktrace:
 [1] top-level scope
   @ none:1
```

"""
DiscreteWeightFunction

basis_symbol(::Type{<:AbstractDiscreteWeightFunction}) = "W"

xs_ws(::Type{<:AbstractDiscreteWeightFunction}) = throw(ArgumentError("No default method"))

# (Gautschi](https://www.cs.purdue.edu/homes/wxg/Madrid.pdf), section 3.1
# compute α_n = <tπ_n,π_n>/<π_n,π_n>, β_n = <π_n,π_n>/<π_{n-1},π_{n-1}>
# wherer <p,q> = ∑_1^N w_k p(k) q(k)
function discrete_stieltjes(W::Type{<:AbstractDiscreteWeightFunction})
    xs, ws = xs_ws(W)

    N = length(xs)
    n = N

    # k = 0 case
    πk_1, πk = zeros(N), ones(N)

    βk = β0 = norm_k = sum(ws) / 1 # <π_0, π_0> = ∑ w_k 1 ⋅ 1
    αk = α0 = ∑(ws[k] * xs[k] for k in eachindex(ws)) / norm_k
    αs = [α0]
    βs = [β0]

    for _ in 1:(n - 1)
        πk1 = (xs .- αk) .* πk - βk * πk_1 # use just computed αk, βk to find π_{k+1} = (x-αk)⋅π_k - βk⋅π_{k-1}
        norm_k1 = ∑(ws[k] * πk1[k]^2 for k in eachindex(ws)) # <π_{k+1}, π_{k+1}>

        # Darboux
        # α_{k+1} = <x ⋅ π_{k+1}, π_{k+1}> / <π_{k+1}, π_{k+1}>,
        # β_{k+1} = <π_{k+1}, π_{k+1}> / <π_k, π_k>,
        αk1 = ∑(ws[k] * xs[k] * πk1[k]^2 for k in eachindex(ws)) / norm_k1
        βk1 = norm_k1 / norm_k

        push!(αs, αk1)
        push!(βs, βk1)

        αk, βk = αk1, βk1
        πk_1, πk = πk, πk1
        norm_k = norm_k1
    end

    (-αs, βs)
end

An(::Type{W}, n) where {W<:AbstractDiscreteWeightFunction} = one(eltype(W))

Bn(::Type{W}, k::Int) where {W<:AbstractDiscreteWeightFunction} =
    discrete_stieltjes(W)[1][k + 1]
Cn(::Type{W}, k::Int) where {W<:AbstractDiscreteWeightFunction} =
    discrete_stieltjes(W)[2][k + 1]

Polynomials.domain(::Type{<:DiscreteWeightFunction}) = Polynomials.Interval(-Inf, Inf)

innerproduct(W::Type{<:AbstractDiscreteWeightFunction}, f, g) =
    ∑(wk * f(xk) * g(xk) for (xk, wk) in zip(xs_ws(W)...))
