"""
    DiscreteWeightFunction

For a discrete measure `dλ = ∑ wᵢ δ(x - xᵢ)` specified through two
vectors `xs` and `ws` a family of monic orthogonal polynomials is
produced through Darboux's formula for `α_n` and `β_n` using the
3-term recurrence defined by `π_{n+1} = (x-α_n)⋅π_n - β_n⋅π_{n-1}`
and the discrete Stieltjes method.

# Example

Discrete Chebyshev by its weight function

```jldoctest
julia> using SpecialPolynomials
julia> SP = SpecialPolynomials
julia> N = 9
julia> xs = collect(0:N-1)
julia> ws = ones(N)
julia> p = DiscreteWeightFunction(xs, ws, [0,0,1], :x)
julia> [SP.beta.(p, 0:N-1) SP.beta.(DiscreteChebyshev{N,Float64}, 0:N-1)]
```


"""
struct DiscreteWeightFunction{S1, S2, R1, R2, T <: Number} <: AbstractWeightFunction{T}
    xs::Vector{S1}
    ws::Vector{S2}
    αs::Vector{R1}
    βs::Vector{R2}
    coeffs::Vector{T}
    var::Symbol
    # no inner constructor here
end

export DiscreteWeightFunction

function DiscreteWeightFunction(xs, ws, coeffs, var=:x)
    N = length(xs)
    N == length(ws) || throw(ArgumentError("Nodes and weights are vectors of the same length"))
    αs, βs = discrete_stieltjes(xs, ws, N)
    DiscreteWeightFunction(xs, ws, αs, βs, coeffs, var)
end

(ch::DiscreteWeightFunction)(x::S) where {S} = orthogonal_polyval(ch, x)
Polynomials.domain(::Type{<:DiscreteWeightFunction}) = Polynomials.Interval(-Inf, Inf)

alpha(p::DiscreteWeightFunction, n) = n < 0 ? 0 : p.αs[n+1]
beta(p::DiscreteWeightFunction, n) = p.βs[n+1]

innerproduct(p::DiscreteWeightFunction, f, g) = ∑(wk * f(xk) * g(xk) for (xk, wk) in zip(p.xs, p.ws))

# (Gautschi](https://www.cs.purdue.edu/homes/wxg/Madrid.pdf), section 3.1
# compute α_n = <tπ_n,π_n>/<π_n,π_n>, β_n = <π_n,π_n>/<π_{n-1},π_{n-1}>
# wherer <p,q> = ∑_1^N w_k p(k) q(k)
function discrete_stieltjes(xs, ws, n)
    N  =  length(xs)
    n > N &&  throw(ArgumentError("n ≤ N is required"))

    # k = 0 case
    πk_1, πk = zeros(N), ones(N)
    
    βk = β0 = norm_k = sum(ws) # <π_0, π_0> = ∑ w_k 1 ⋅ 1
    αk = α0 = ∑(ws[k] * xs[k] for k in eachindex(ws)) / norm_k
    αs = [α0]
    βs = [β0]

    for _ ∈ 1:(n-1)

        πk1 = (xs .- αk) .* πk - βk * πk_1 # use just computed αk, βk to find π_{k+1} = (x-αk)⋅π_k - βk⋅π_{k-1}
        norm_k1 =  ∑(ws[k] * πk1[k]^2  for k in eachindex(ws)) # <π_{k+1}, π_{k+1}>
        
        # Darboux
        # α_{k+1} = <x ⋅ π_{k+1}, π_{k+1}> / <π_{k+1}, π_{k+1}>,
        # β_{k+1} = <π_{k+1}, π_{k+1}> / <π_k, π_k>,
        αk1 = ∑(ws[k] * xs[k] * πk1[k]^2 for k in eachindex(ws)) / norm_k1
        βk1 = norm_k1/norm_k
        
        push!(αs, αk1)
        push!(βs, βk1)

        αk, βk = αk1, βk1
        πk_1, πk = πk, πk1
        norm_k = norm_k1

    end

    (αs, βs)
end
       
        
##################################################

const ∑ = sum

function innerproduct(P::Type{<:DiscreteOrthogonalPolynomial}, f, g)
    dom = domain(P)
    fn = x -> f(x) * g(x) * weight_function(P)(x)
    a, b = first(dom), last(dom)
    if !isinf(a) && !isinf(b)
        return ∑(fn(x)  for x in first(dom):last(dom))
    else
        ## what to do if infinite
    end
end
