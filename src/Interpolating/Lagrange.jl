"""
    Lagrange(xs, [ws], coeffs, [var])

Lagrange interpolation of points `(xᵢ, fᵢ)` for `i ∈ 0..n`.

* `xs`, `coeffs`: the interpolating coordinates.
* `ws`: weights used in the barycentric representation. (From `SpecialPolynomials.lagrange_barycentric_weights` or `SpecialPolynomials.lagrange_barycentric_nodes_weights`.)
* var: the polynomial indeterminate

# Extended help

The Lagrange interpolation of points `(xᵢ, fᵢ)` for `i ∈ 0:n` is the polynomial `p(x) = ∑ᵢ lⱼ(x) fⱼ`.

The basis vectors `lⱼ(x)` are `1` on `xⱼ` and `0` on `xᵢ` when `i ≠ j`. That is, `lⱼ(x) = Π_{i ≠ j}(x-xᵢ)/Π_{i ≠j}(xⱼ-xᵢ)`. These can be rewritten in terms of weights, `wⱼ`, depending on the `xᵢ` only, yielding `lⱼ = l(x) wⱼ/(x - xⱼ)` with `l(x) = Π(x-xᵢ)`. Going further, yields the barycentric formula:

`p(x) = (∑ wⱼ / (x - xⱼ) ⋅ fⱼ) /  (∑ wⱼ / (x - xⱼ) )`.

This representation has several properties, as detailed in Berrut and Trefethen [Barycentric Lagrange Interpolation](https://doi.org/10.1137/S0036144502417715).

## Examples

```jldoctest Lagrange
julia> using Polynomials, SpecialPolynomials

julia> p =  Lagrange([1,2,3], [1,2,3])
Lagrange(1⋅ℓ_0(x) + 2⋅ℓ_1(x) + 3⋅ℓ_2(x))

julia> p.([1,2,3]) # the coefficients
3-element Vector{Int64}:
 1
 2
 3

julia> convert(Polynomial,  p)
Polynomials.Polynomial(1.0*x)
```

The instances hold the nodes and weights, which are necessary for
representation, so the type alone can not be used for functions such
as `variable` or `convert(Lagrange, ...)`. For the former we can use
an instance, for the latter we can use `fit`:

```jldoctest Lagrange
julia> p =  Lagrange([1,2,3], [1,2,3])
Lagrange(1⋅ℓ_0(x) + 2⋅ℓ_1(x) + 3⋅ℓ_2(x))

julia> variable(p)
Lagrange(1⋅ℓ_0(x) + 2⋅ℓ_1(x) + 3⋅ℓ_2(x))

julia> q = Polynomial([0,0,1])
Polynomials.Polynomial(x^2)

julia> qq = fit(Lagrange, p.xs, p.ws, q)
Lagrange(1⋅ℓ_0(x) + 4⋅ℓ_1(x) + 9⋅ℓ_2(x))

julia> convert(Polynomial, qq)
Polynomials.Polynomial(1.0*x^2)
```

For a given set of nodes,
`SpecialPolynomials.lagrange_barycentric_weights` can compute the
weights.  For all but modest values of `n`, interpolating polynomials
suffer from the Runge phenomenon unless the nodes are well
chosen. (They should have asymptotic density of `1/√(1-x^2)` over
`[-1,1]`.) For `P=Chebyshvev` and `P=ChebyshevU`, the function
`SpecialPolynomials.lagrange_barycentric_nodes_weights(P, n)` will
return a good choice of `n+1` points over `[-1,1]` along with
precomputed weights.

```jldoctest Lagrange
julia> xs, ws = SpecialPolynomials.lagrange_barycentric_nodes_weights(Chebyshev, 64);


julia> f(x) = exp(-x)*sinpi(x)
f (generic function with 1 method)

julia> p = fit(Lagrange, xs, f.(xs));


julia> degree(p)
64

julia> maximum(abs.(f(x) - p(x) for x in range(-1, 1, length=20))) <= 1e-14
true
```

!!! note
    The above example is more directly done through
    `fit(Chebyshev, f, 64)`, though the resulting polynomial will
    reference a different basis.

"""
struct Lagrange{N,S<:Number,R<:Number,T<:Number,X} <: AbstractInterpolatingPolynomial{T,X}
    xs::Vector{S}
    ws::Vector{R}
    coeffs::Vector{T}
    var::Symbol
    function Lagrange(
        xs::Vector{S},
        ws::Vector{R},
        coeffs::Vector{T},
        var::Polynomials.SymbolLike=:x,
    ) where {S,R,T}
        xs = unique(xs)
        N = length(xs)
        N == length(coeffs) ||
            throw(ArgumentError("the unique xs and the coeffs must have the same length"))
        new{N,S,R,T,Symbol(var)}(xs, ws, coeffs)
    end

    # no weights
    function Lagrange(xs::Vector{S},
                      coeffs::Vector{T},
                      var::Polynomials.SymbolLike=:x) where {S,T}
        xs = unique(xs)
        N = length(xs)
        ws = lagrange_barycentric_weights(xs)
        R = eltype(ws)
        length(coeffs) == length(xs) ||
            throw(ArgumentError("the unique xs and the coeffs must have the same length"))
        new{N,S,R,T,Symbol(var)}(xs, ws, coeffs)
    end
end

export Lagrange

basis_symbol(::Type{<:Lagrange}) = "ℓ"

## Boilerplate code reproduced here, as there are three type parameters
Base.convert(::Type{P}, p::P) where {P<:Lagrange} = p
Base.convert(::Type{Lagrange{N,S,R,T}}, p::Lagrange) where {N,S,R,T} =
    Lagrange{N,S,R,T}(p.xs, p.ws, coeffs(p), p.var)
Base.promote_rule(
    ::Type{Lagrange{N,S,R,T}},
    ::Type{Lagrange{N,S,R,Q}},
) where {N,S,R,T,Q<:Number} = Lagrange{N,S,R,promote_type(T, Q)}
Base.promote_rule(::Type{Lagrange{N,S,R,T}}, ::Type{Q}) where {N,S,R,T,Q<:Number} =
    Lagrange{N,S,R,promote_type(T, Q)}

Polynomials.domain(::Type{<:Lagrange}) = Polynomials.Interval(-Inf, Inf)

function Polynomials.variable(
    p::Lagrange{S,R,T},
    var::Polynomials.SymbolLike=:x,
) where {S,R,T}
    _fit(Lagrange, p.xs, p.ws, p.xs, var)
end

function Polynomials.one(p::Lagrange)
    _fit(
        Lagrange,
        p.xs,
        p.ws,
        ones(eltype(p.xs), length(p.xs)),
        Polynomials.indeterminate(p),
    )
end

Polynomials.zero(p::Lagrange{N,S,R,T}) where {N,S,R,T} = 0 * p

##  Evaluation
##  From https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
function Polynomials.evalpoly(x, p::Lagrange)
    ws, xs, cs = p.ws, p.xs, coeffs(p)
    a = b = zero(x / 1)
    for (xⱼ, wⱼ, cⱼ) in zip(xs, ws, cs)
        Δ = (x - xⱼ)
        iszero(Δ) && return cⱼ # xj ∈ xs, so sum c_j l_i(xj) = sum c_j δ_{ij} = c_j
        l = wⱼ / Δ
        a += l * cⱼ
        b += l
    end
    a / b
end
(p::Lagrange)(x::Number) = evalpoly(x, p)

## convert to poly
## avoids division in barycentric form, which isn't available for Polynomial  arguments, such as `x-xs[j]`
## This has numeric issues, unlike the barycentric form used for evaluation
function Base.convert(Q::Type{<:Polynomial}, p::Lagrange{N,S,R,T}) where {N,S,R,T}
    q = zero(Q) / 1
    x = variable(q)
    xs, ws, cs = p.xs, p.ws, coeffs(p)

    if length(xs) == 1
        return cs[1] + 0 * x #Polynomial(cs[1]*ones(T,1))
    else
        for (i, (wᵢ, cᵢ)) in enumerate(zip(ws, cs))
            lᵢ = prod(x - xⱼ for (j,xⱼ) in enumerate(xs) if j != i)
            q += lᵢ * wᵢ * cᵢ
        end
    end
    q
end

## ----- nodes and weights -----

## return w⁻¹s = {1/w_j, 0 <= j <= n} as 1 based ws[j+1] = w_j
"""
    lagrange_barycentric_weights(xs)

Return `[1/w_i for 0 <i <= n]` where `wi = ∏_{j≠i} (xi - xj)`. If `xs` is a range, the weights have a closed form formula, though for large `n` interpolation will be unstable.
"""
function lagrange_barycentric_weights(xs)
    n = length(xs)
    ws = ones(eltype(xs), n)
    for j in 2:n
        for k in 1:j
            ws[k] = (xs[k] - xs[j]) * ws[k]
        end
        ws[j] = prod(xs[j] - xs[i] for i in 1:(j - 1))
    end
    1 ./ ws
end

# https://people.maths.ox.ac.uk/trefethen/barycentric.pdf (5.1)
function lagrange_barycentric_weights(
    xs::Union{AbstractRange{T},AbstractUnitRange{T}},
) where {T}
    n = length(xs) - 1
    ws = zeros(float(T), n + 1)
    o = one(T)
    for j in 0:n
        ws[1 + j] = o * binomial(n, j)
        o *= -1
    end
    ws
end

##  add a new point
## check for  xk ∈ xs?
"""
    update_lagrange_barycentric_weights(ws, xs, xn1)

Adds `xn1` to the nodes `xs` and updates the weights `ws` accordingly.
"""
function update_lagrange_barycentric_weights(ws, xs, xn1)
    ys = [wj / (xj - xn1) for (wj, xj) in zip(ws, xs)]
    wn1 = one(eltype(ws))
    for xj in xs
        wn1 *= (xn1 - xj)
    end
    push!(ys, 1 / wn1)
    ys
end

# combine two non-equal pairs of (xs, ws), (us, vs)
function merge_nodes_weights(
    p1::Lagrange{N,S,R,T,X},
    p2::Lagrange{M,S1,R1,T1,Y},
) where {N,S,R,T,X,M,S1,R1,T1,Y}

    same_nodes(p1, p2) && return (p1.xs, p1.ws)

    p, q = N >= M ? (p1, p2) : (p2, p1)
    S2 = typeof(one(S) * one(S1) / 1)
    xs = convert(Vector{S2}, copy(p.xs))
    ws = copy(p.ws)

    # new_xs ## XXX this is where the work is
    new_xs = _new_nodes(xs, q.xs)
    for y in new_xs
        ws = update_lagrange_barycentric_weights(ws, xs, y)
        push!(xs, y)
    end

    xs, ws
end

"""
    lagrange_barycentric_nodes_weights(::Type{<:SpecialPolynomial}, n::Int)

Return a collection of `n+1` nodes and weights the given family of
polynomials. Defined for `ChebyshevT` and `ChebyshevU` familes, as there
are `O(n)` algorithms available. Otherwise, `lagrange_barycentric_weights`
is needed.


"""
lagrange_barycentric_nodes_weights(::Type{<:AbstractSpecialPolynomial}, n::Int) =
    throw(MethodError())


# do we have the same nodes (which makes easier to combine)
same_nodes(p1::Lagrange, p2::Lagrange) = false
same_nodes(p1::Lagrange{N,S}, p2::Lagrange{N,S}) where {N,S} =
    p1.xs == p2.xs

## ----- arithmetic operations

# call op on coeffs when ps all match (assumed)
function _lagrange_op(f, ps...)
    p = first(ps)
    xs, ws, X = p.xs, p.ws, Polynomials.indeterminate(p)
    Lagrange(xs, ws, f([p.coeffs for p ∈ ps]...), X)
end

function Base.:+(
    p1::Lagrange{N,S,R,T,X},
    p2::Lagrange{M,S1,R1,T1,Y},
) where {N,S,R,T,X,M,S1,R1,T1,Y}

    assert_same_variable(p1, p2) || throw(ArgumentError("`p` and `q` have different indeterminate"))

    ## same nodes/weights
    same_nodes(p1, p2) && return _lagrange_op(+, p1, p2)

    ## o/w we have to merge, evaluate, fit
    xs, ws = merge_nodes_weights(p1, p2)
    ys = p1.(xs) + p2.(xs)
    return Lagrange(xs, ws, ys, X)
end

function Base.:*(
    p1::Lagrange{N,S,R,T,X},
    p2::Lagrange{M,S1,R1,T1,Y},
) where {N,S,R,T,X,M,S1,R1,T1,Y}

    assert_same_variable(p1, p2) || throw(ArgumentError("`p` and `q` have different indeterminate"))

    # same nodes/weigths
    same_nodes(p1, p2) && return _lagrange_op(*, p1, p2)

    xs, ws = merge_nodes_weights(p1, p2)
    ys = p1.(xs) .* p2.(xs)
    return Lagrange(xs, ws, ys, X)

end

Base.:-(p::P) where {P<:Lagrange} = _lagrange_op(-, p)

# Scalar ops
function Base.:+(p::Lagrange, c::Number)
    xs, ws, X = p.xs, p.ws,  Polynomials.indeterminate(p)
    ys = c * one.(ws)
    q = Lagrange(xs, ws, ys, X)
    p + q
end

Base.:*(c::Number, p::Lagrange) = p * c
function Base.:*(p::P, c::Number) where {P<:Lagrange}
    xs, ws, X = p.xs, p.ws, Polynomials.indeterminate(p)
    ys = c * coeffs(p)
    Lagrange(xs, ws, ys, X)
end


## derivative

## The derivative can be computed exactly by
## s(x) = ∑ wⱼ(x-xᵢ)/(x-xⱼ)
## l′ⱼ s(x) + lⱼ s'(x) = wⱼ((x-xᵢ)/(x-xⱼ))' for some chosen xᵢ
## However, he we *fit* the derivative using a polynomial to *evaluate* it, as
## l′ⱼ(xᵢ) = wⱼ/wᵢ / (xᵢ - xⱼ) are known.
## Using p(x) = ∑ lⱼ fⱼ we get p'(xᵢ) = ∑ l′ⱼ(xᵢ)fⱼ so we have (xᵢ, p′(xᵢ)) values easily computed.


# 0-based enumeration p
l′(p::Lagrange, j, i) = (p.ws[1+j] / p.ws[1+i]) / (p.xs[1+i] - p.xs[1+j])
l′(p::Lagrange, j) = - sum(l′(p, j′, j) for j′ ∈ eachindex(p) if j′ ≠ j)

"""
    derivative(p::Lagrange)

In [C. Schneider and W. Werner](https://doi.org/10.2307/2008095) is a formula for the ``k``th derivative of polynomial in barycentric form for x not in {xᵢ}:

r^(k)(x)/k! = ∑(wᵢ/(x-xᵢ)⋅r[x,x,...x,xᵢ]) / ∑wᵢ/(x - xᵢ)

We don't use that here, rather, we use the fact that the derivative evaluated at xᵢ can be evaluated at the basis polynomials using the formulas to evaluate `lⱼ′(xᵢ)`, where `p=∑lⱼ fⱼ`. , so we can evaluate the derivative at the same nodes as the polynomial and fit the polynomial that way. Up to numeric issues, this should be the same polynomial.


"""
function Polynomials.derivative(p::Lagrange)

    p′s = zero.(coeffs(p))

    for i ∈  eachindex(p)
        # p′(xᵢ) = ∑ l′ⱼ(xᵢ) fⱼ
        for (j, fⱼ) ∈ pairs(p)
            l′ⱼᵢ = i == j ? l′(p, j) : l′(p, j, i)
            p′s[1+i] += l′ⱼᵢ * fⱼ
        end
    end

    xs, ws, X = p.xs, p.ws, Polynomials.indeterminate(p)
    Lagrange(xs, ws, p′s, X)
end

## ----- fitting -----

## short cut if weights are known.
_fit(P::Type{Lagrange}, xs, ws, ys, var=:x) = Lagrange(xs, ws, ys, var)

function Polynomials.fit(
    P::Type{<:Lagrange},
    xs::AbstractVector{S},
    ys::AbstractVector{T};
    var=:x,
) where {S,T}
    ws = lagrange_barycentric_weights(xs)
    _fit(P, xs, ws, ys, var)
end

function Polynomials.fit(P::Type{<:Lagrange},
                         xs::AbstractVector{S},
                         f; var=:x) where {S}
    ws = lagrange_barycentric_weights(xs)
    _fit(P, xs, ws, f.(xs), var)
end

# fit with weights specified
function Polynomials.fit(P::Type{<:Lagrange},
                         xs::AbstractVector{S}, ws::AbstractVector,
                         f; var=:x) where {S}
    _fit(P, xs, ws, f.(xs), var)
end
