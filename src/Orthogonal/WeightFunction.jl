## General  weight function
## Wheeler algorithm, https://authors.library.caltech.edu/86868/1/1.4822929.pdf
## algorithm not considered, just text

abstract type AbstractWeightFunction{T} <: AbstractOrthogonalPolynomial{T} end

"""
    WeightFunction{P,F,T}

A type for orthogonal polynomials relative to some weight
function. The Wheeler or modified Chebyshev algorithm
([Gautschi](https://www.cs.purdue.edu/homes/wxg/Madrid.pdf),
[Press and Teukolsky](https://doi.org/10.1063/1.4822929))
is used
to generate the three-term recurrence relation.  This uses a known
family of orthogonal polynomials over the domain of the weight
function.  These are specified through `pis`.

#  Example. 

Elliptic orthogonal polynomials on  `[-1,1]`. Demo 2 of Gautschi.

```jldoctest WeightFunction
julia> using Polynomials, SpecialPolynomials

julia> N,  ω² = 40, 0.999
(40, 0.999)

julia> w(t) = ((1-ω²*t^2)*(1-t^2))^(-1/2)
w (generic function with 1 method)

julia> πs = Chebyshev{Float64}
Chebyshev{Float64}

julia> p = WeightFunction(πs, w)
WeightFunction(1⋅e_1(x))

julia> αs = SpecialPolynomials.alpha.(p, 0:5); # α_0, α_1, ..., α_5

julia> βs = SpecialPolynomials.beta.(p, 0:5);  # β_0, β_1, ..., β_5

julia> round.([αs βs], digits=8)
6×2 Array{Float64,2}:
  0.0  9.68226
 -0.0  0.793782
 -0.0  0.119868
 -0.0  0.22704
 -0.0  0.241061
 -0.0  0.245429
```


The main computation involved in this is the modified moment `ν_l = ∫
π_l dw`. If πs has a fast gauss nodes available, then those are used,
otherwise `QuadGK.quadgk` is. For  some examples, this  can be overridden. This one is from  Press and Teukolsky, where the modified moments are given through  the  function  `v(j)` defined below.

```jldoctest
julia> using Polynomials, SpecialPolynomials, SpecialFunctions

julia> πs = ShiftedLegendre
ShiftedLegendre

julia> w(t) = -log(t)
w (generic function with 1 method)

julia> p = WeightFunction(πs, w)
WeightFunction(1⋅e_1(x))

julia> v(j) = iszero(j) ? 1 : (-1)^j*gamma(j+1)^2/(j*(j+1)*gamma(2j+1))
v (generic function with 1 method)

julia> WW = typeof(w)
typeof(w)

julia> SpecialPolynomials.modified_moment(w::W, j) where {W <: WW} = v(j)

julia> αs = SpecialPolynomials.alpha.(p, 0:5);

julia> βs = SpecialPolynomials.beta.(p, 0:5);

julia> [αs βs]
6×2 Array{Float64,2}:
 0.25      1.0
 0.464286  0.0486111
 0.485482  0.0586848
 0.492103  0.0607286
 0.495028  0.061482
 0.49658   0.0618408
```

"""
struct WeightFunction{P, F, T <: Number} <: AbstractWeightFunction{T}
    pis::P
    W::F
    coeffs::Vector{T}
    var::Symbol
    # no inner constructor here
end

WeightFunction(pis::Type{<:AbstractOrthogonalPolynomial}, W, coeffs::Vector{T}=[0,1]; var::Symbol=:x) where {T} = WeightFunction(pis, W, coeffs, var)


export WeightFunction

## conversion functions ....

basis_symbol(::Type{<:AbstractWeightFunction}) = "W"

Polynomials.domain(::Type{<:WeightFunction}) = throw(ArgumentError("Must pass instance, not type"))
Polynomials.domain(w::WeightFunction) = domain(w.pis)
weight_function(::Type{<:WeightFunction}) = throw(ArgumentError("Must pass instance, not type"))
weight_function(w::WeightFunction) = w.W


(ch::WeightFunction)(x::S) where {S} = orthogonal_polyval(ch, x)


## An,Bn,Cn for the orthogonal polynomials for the given weight function
An(w::AbstractWeightFunction, n) = 1.0
Bn(w::AbstractWeightFunction, k) = - alpha(w, k)
Cn(w::AbstractWeightFunction, k) = - beta(w, k)

## methods to work with instances, not types
P0(p::P, x) where {P <: AbstractWeightFunction} = one(x)
P1(p::P, x) where {P <: AbstractWeightFunction} = An(p,0)*x + Bn(p,0)

function jacobi_matrix(p::P, n) where {P <: AbstractWeightFunction}
    LinearAlgebra.SymTridiagonal([alpha(p,i) for i in 0:n-1], [sqrt(beta(p,i)) for i in 1:n-1])
end

function gauss_nodes_weights(p::P, n) where {P <: AbstractWeightFunction}
    J = jacobi_matrix(p, n)
    eig = eigen(J, extrema(p)...)
    # Is this necessary?
    nm  = 1 # diag(eig.vectors * eig.vectors')
    wts =  beta(p,0) * eig.vectors[1,:] / nm
    (eig.values,  wts)
end

function Base.extrema(p::P) where {P <: WeightFunction}
    dom = domain(p)
    first(dom), last(dom)
end


## Try to speed up  quadgk by not specializing on F
mutable struct Wrapper
    F
end
(F::Wrapper)(x) = F.F(x)

# integrate fdw over region specified by domain(P)
function ∫(fdw, P; kwargs...)
    dom = domain(P)
    a, b = first(dom), last(dom)
    if !first(Polynomials.inclusivity(dom))
        a += eps(float(one(a)))
    end
    if !last(Polynomials.inclusivity(dom))
        b -= eps(float(one(b)))
    end

    F = Wrapper(fdw)

    quadgk(F, a, b; kwargs...)[1]
end
    
function innerproduct(p::WeightFunction, f, g)
    P = p.pis
    fgdw = x -> f(x) * g(x) * weight_function(p)(x)
    ∫(fgdw, P)
end



## Compute the modified moment
## v_j =  ∫ π_j W(x) dx
## (The moment is ∫ x^j W(x) dx)
## This can  be overridden when an exact calculuation is known, e.g.:
## q0 = WeightFunction(P,f,[1],:x)
## WW = typeof(q0)
## SpecialPolynomials.modified_moment(w::W, j) where {W <: WW} = v(j)
#@memoize
function modified_moment(w::WeightFunction,  j; n = 10_000, kwargs...)
    P = w.pis
    val = Val(has_fast_gauss_nodes_weights(P))
    _modified_moment(val, w, P, j; n=n, kwargs...)
end

# brute force, not  adaptive quadrature
# 10_000 might be too big/too small
# XXX this could use some theory behind it.
function  _modified_moment(use_gauss::Val{true}, w, P, j; n=10_000, kwargs...)
    xs,  ws = gauss_nodes_weights(P, n)
    wf = weight_function(P)
    π_jdw = x -> monic_orthogonal_basis_polyval(P, j, x)/wf(x)*w.W(x)
    dot(π_jdw.(xs), ws)
end

function  _modified_moment(use_gauss::Val{false}, w, P, j; n=0, kwargs...)
    π_jdw = x -> monic_orthogonal_basis_polyval(P, j, x)*w.W(x)
    return ∫(π_jdw, P; kwargs...)
end



## If the moments, v_k = ∫_a^b x^k W dx are known, then Golub and
## Welsh
## (https://web.stanford.edu/class/cme335/S0025-5718-69-99647-1.pdf)
## give an algorithm tot find An, Bn, Cn.
## However, the moments arer numerically unstable.
## The Wheeler approach, detailed in Press and Tuloksky (https://authors.library.caltech.edu/86868/1/1.4822929.pdf)
## uses a different formula. In [Gautschi](https://www.cs.purdue.edu/homes/wxg/Madrid.pdf), p10, this is called
## the mmodified Chebyshev algorithm
#
## σ(k,l) = <p_k, π_l>, where π is the reference set of Orthogonal polynomials, and p_k are the
## orthogonal polynomials for w being generated
@memoize function sigma(w::W, k, l) where {W <: WeightFunction}
    P = w.pis

    T = eltype(w)
    k == -1 && return zero(T)/1
    k == 0  && return modified_moment(w, l) # ν_l = ∫π_l dw
    k > l   && return zero(T)/1

    # otherwise get through recursion (this uses Gautschi's formula, P&T have apparent error
    sigma(w, k-1, l+1) - (alpha(w, k-1) - alpha(P, l)) * sigma(w, k-1, l) - (beta(w, k-1))*sigma(w, k-2, l) + beta(P,l)*sigma(w, k-1, l-1)
    
end

## Compute the α, β functions for the 3-point recurrence for the orthogonal polynomials relative to the weight function
@memoize function alpha(w::WeightFunction, k)
    P = w.pis
    if k == 0
        ν₀, ν₁ =  modified_moment(w, 0), modified_moment(w, 1)
        alpha(P, 0) + ν₁ / ν₀
    else
        alpha(P, k) - sigma(w, k-1, k)/sigma(w, k-1, k-1)  + sigma(w, k, k+1)/sigma(w, k, k) 
    end
end

@memoize function beta(w::WeightFunction, k)
    iszero(k) && return innerproduct(w, one, one)
    (sigma(w, k, k)/sigma(w, k-1, k-1))
end
