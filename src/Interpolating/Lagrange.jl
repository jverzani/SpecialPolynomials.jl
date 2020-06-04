"""
    Lagrange{N, S, R, T}

Represent a polynomial in Lagrange form using nodes `xs`, weights `ws`, and coefficients `coeffs`.
The Lagrange form does polynomial interpolation between `xs` and `ys` through `p(x) = Σ_{0..n} ℓ_i(x) y_i`,
where if `ℓ(x) = prod(x-x_i)`, `w_i = 1/prod_{j≠i}(x_i - x_j)`, then `ℓ_i(x) = ℓ(x) w_i/(x-x_i)`. The `ℓ_i`
satisfy `ℓ_i(x_j) = δ_{ij}`, so the coefficients are just the `ys`.

```jldoctest Lagrange
julia> using Polynomials, SpecialPolynomials

julia> p =  Lagrange([1,2,3], [1,2,3])
Lagrange(1⋅ℓ^2_0(x) + 2⋅ℓ^2_1(x) + 3⋅ℓ^2_2(x))

julia> p.([1,2,3]) # the coefficients
3-element Array{Int64,1}:
 1
 2
 3

julia> convert(Polynomial,  p)
Polynomial(1.0*x)
```

The instances hold the nodes and weights, which are necessary for
representation, so the type alone can not be used for functions such
as `variable` or `convert(Lagrange, ...)`. For the former we can  use an instance, for the latter we can use `fit`:

```jldoctest Lagrange
julia> p =  Lagrange([1,2,3], [1,2,3])
Lagrange(1⋅ℓ^2_0(x) + 2⋅ℓ^2_1(x) + 3⋅ℓ^2_2(x))

julia> variable(p)
Lagrange(1⋅ℓ^2_0(x) + 2⋅ℓ^2_1(x) + 3⋅ℓ^2_2(x))

julia> q = Polynomial([0,0,1])
Polynomial(x^2)

julia> qq = fit(Lagrange, p.xs, q)
Lagrange(1⋅ℓ^2_0(x) + 4⋅ℓ^2_1(x) + 9⋅ℓ^2_2(x))

julia> convert(Polynomial, qq)
Polynomial(1.0*x^2)
```

Interpolating polynomials suffer from the Runge phenomenon unless the
nodes are well chosen. For `P=Chebyshvev` and `P=ChebyshevU`, the
function `SpecialPolynomials.lagrange_barycentric_nodes_weights(P, n)`
will return a good choice of `n+1` points over `[-1,1]` along with
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
    The above example  is  more directly  done through `fit(Chebyshev, f, 64)`, though  the resulting
    polynomial will reference a different  basis.

"""
struct Lagrange{N, S<:Number, R <: Number, T <: Number} <: AbstractInterpolatingPolynomial{T}
    xs::Vector{S}
    ws::Vector{R}
    coeffs::Vector{T}
    var::Symbol
    function  Lagrange(xs::Vector{S}, ws::Vector{R}, coeffs::Vector{T}, var::Symbol=:x) where {S,R,T}
        xs = unique(xs)
        N = length(xs)
        N == length(coeffs) || throw(ArgumentError("the unique xs and the coeffs must have the same length"))
        new{N, S,R,T}(xs, ws, coeffs, var)
    end
    function  Lagrange(xs::Vector{S}, coeffs::Vector{T}, var::Symbol=:x) where {S,T}
        xs = unique(xs)
        N = length(xs)
        ws = lagrange_barycentric_weights(xs)
        R = eltype(ws)
        length(coeffs) == length(xs) || throw(ArgumentError("the unique xs and the coeffs must have the same length"))
        new{N, S, R, T}(xs, ws, coeffs, :x)
    end
end

export Lagrange

basis_symbol(::Type{Lagrange{N,S,R,T}}) where {N,S,R,T} = "ℓ^$(N-1)"

## Boilerplate code reproduced here, as there are three type parameters
Base.convert(::Type{P}, p::P) where {P <: Lagrange} =  p
Base.convert(::Type{Lagrange{N,S,R,T}}, p::Lagrange) where {N,S,R,T} = Lagrange{N,S,R,T}(p.xs, p.ws, coeffs(p), p.var)
Base.promote_rule(::Type{Lagrange{N, S, R, T}}, ::Type{Lagrange{N, S, R, Q}}) where {N, S, R, T, Q <: Number} = Lagrange{N, S, R, promote_type(T,Q)}
Base.promote_rule(::Type{Lagrange{N, S, R, T}}, ::Type{Q}) where {N, S, R, T, Q <: Number} = Lagrange{N, S, R, promote_type(T,Q)}






Polynomials.domain(::Type{<:Lagrange}) = Polynomials.Interval(-Inf, Inf)
function Polynomials.variable(p::Lagrange{S, R, T}, var::Polynomials.SymbolLike=:x) where {S, R, T}
    _fit(Lagrange, p.xs, p.ws, p.xs, var)
end
function Polynomials.one(p::Lagrange)
    _fit(Lagrange, p.xs, p.ws, ones(eltype(p.xs),length(p.xs)), p.var)
end
Polynomials.zero(p::Lagrange{N, S, R, T}) where {N,S,R,T} = 0*p

##  Evaluation
##  From https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
function (p::Lagrange)(x::Number)
    ws, xs, cs = p.ws, p.xs, coeffs(p)
    a = b = zero(x/1)
    for j in eachindex(xs)
        Δ = (x - xs[j])
        iszero(Δ) && return cs[j] # xj ∈ xs, so sum c_j l_i(xj) = sum c_j δ_{ij} = c_j
        l = ws[j]/Δ
        a += l * cs[j]
        b += l
   end
   a/b
end

## convert to poly
## avoids division in barycentric form, which isn't available for Polynomial  arguments, such as `x-xs[j]`
## This has numeric issues, unlike the barycentric form used for evaluation
function  Base.convert(Q::Type{<:Polynomial},  p::Lagrange{N, S,R,T})  where  {N, S, R, T}
    q = zero(Q)/1
    x =  variable(q)
    xs = p.xs
    cs = coeffs(p)
    ws = p.ws

    if length(xs) == 1
       return cs[1] + 0*x #Polynomial(cs[1]*ones(T,1))
    else
        l = fromroots(Polynomial, xs)
        for i in  eachindex(xs)
            li = prod(x - xs[j] for j in eachindex(xs) if j != i);
            q += li * ws[i] * cs[i]
        end
    end
    q
end

## return w⁻¹s = {1/w_j, 0 <= j <= n} as 1 based ws[j+1] = w_j
"""
   lagrange_barycentric_weights(xs)

return [1/w_i for 0 <i <= n] where wi = ∏_{j≠i} (xi - xj)
"""
function lagrange_barycentric_weights(xs)
    n =  length(xs)
    ws =  ones(eltype(xs), n)
    for j in 2:n
        for k in 1:j
            ws[k] =  (xs[k] - xs[j]) * ws[k]
        end
        ws[j] = prod(xs[j]-xs[i] for i in 1:j-1)
    end
    1 ./ ws
end

# https://people.maths.ox.ac.uk/trefethen/barycentric.pdf (5.1)
function lagrange_barycentric_weights(xs::Union{AbstractRange{T},AbstractUnitRange{T}}) where {T}
    n = length(xs)-1
    ws = zeros(float(T), n+1)
    o = one(T)
    for j = 0:n
        ws[1+j] = o * binomial(n, j)
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
    ys = [wj/(xj - xn1) for (wj,xj) in zip(ws, xs)]
    wn1 = one(eltype(ws))
    for xj in xs
        wn1 *= (xn1 - xj)
    end
    push!(ys, 1/wn1)
    ys
end

"""
    lagrange_barycentric_nodes_weights(::Type{<:SpecialPolynomial}, n::Int)

Return a collection of `n+1` nodes and weights the given family of
polynomials. Defined for `ChebyshevT` and `ChebyshevU` familes, as there
are `O(n)` algorithms available. Otherwise, `lagrange_barycentric_weights`
is needed.


"""
lagrange_barycentric_nodes_weights(::Type{<:AbstractSpecialPolynomial}, n::Int) = throw(MethodError())

## can add easily if the nodes are shared
## otherwise, this needs to interpolate
function Base.:+(p1::Lagrange{N,S,R,T}, p2::Lagrange{M,S1,R1,T1}) where {N,S,R,T,M,S1,R1,T1}
    p1.var == p2.var || throw(ArgumentError("p1 and p2 must share the same variable"))

    # can just add if using the same set of nodes
    if N == M && p1.xs == p2.xs
        cs =  coeffs(p1) +  coeffs(p2)
        return Lagrange(p1.xs, p1.ws, cs, p1.var)
    else
        p,q = N >= M ? (p1, p2) : (p2, p1)
        xs, ws, cs = p.xs, p.ws, copy(coeffs(p))
        for i in eachindex(xs)
            cs[i] += q(xs[i])
        end
        return _fit(Lagrange, xs,  ws, cs, p1.var)
    end
end

## XXX is there a better way?
## This refits using  updated weights,
## but this can be wildly inaccurate
function Base.:*(p1::Lagrange{N,S,R,T}, p2::Lagrange{M,S1,R1,T1}) where {N,S,R,T,M,S1,R1,T1}

    p1.var == p2.var || throw(ArgumentError("p1 and p2 must share the same variable"))

    p, q =  N >= M ? (p1, p2)  : (p2, p1)
    S2 = typeof(one(S) * one(S1)/1)
    xs = convert(Vector{S2}, copy(p.xs))
    ws = copy(p.ws)

    # new_xs ## XXX this is where the work is
    new_xs = _new_nodes(xs, q.xs)
    for y in new_xs
        ws = update_lagrange_barycentric_weights(ws, xs, y)
        push!(xs, y)
    end
    ys = p.(xs) .* q.(xs)
    _fit(Lagrange, xs,  ws, ys, p1.var)

end


# handle scalar case
function Base.:+(p::Lagrange, c::Number)
    q = Lagrange(p.xs[1:1], [c], p.var)
    p + q
end

Base.:*(c::Number, p::Lagrange) = p*c
function Base.:*(p::P, c::Number) where {P <: Lagrange}
    Lagrange(p.xs, c*coeffs(p), p.var)
end

Base.:-(p::P) where {P<:Lagrange} = (-1)*p



## short cut if weights are known.
_fit(P::Type{Lagrange}, xs, ws, ys, var=:x) = Lagrange(xs, ws, ys, var)

function Polynomials.fit(P::Type{<:Lagrange},
                         xs::AbstractVector{S},
                         ys::AbstractVector{T};
                         var = :x) where {S, T}
    ws = lagrange_barycentric_weights(xs)
    _fit(P, xs, ws, ys, var)
end

function Polynomials.fit(P::Type{<:Lagrange}, xs::AbstractVector{S}, f; var=:x) where {S}
     ws = lagrange_barycentric_weights(xs)
    _fit(P, xs, ws, f.(xs), var)
end
