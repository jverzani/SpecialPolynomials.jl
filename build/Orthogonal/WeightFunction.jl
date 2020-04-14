## General  weight function
## Wheeler algorithm, https://authors.library.caltech.edu/86868/1/1.4822929.pdf
## algorithm not considered, just text

abstract type AbstractWeightFunction{T} <: AbstractOrthogonalPolynomial{T} end

"""
    WeightFunction{P,F,T}

A type for orthogonal polynomials relative to some weight
function. The Wheeler or modified Chebyshev algorithm
([Gautschi](https://www.cs.purdue.edu/homes/wxg/Madrid.pdf)) is used
to generate the three-term recurrence relation.  This uses a known
family of orthogonal polynomials over the domain of the weight
function.  These are specified through `pis`.

#  Example. 

Elliptic orthogonal polynomials on  `[-1,1]`. Demo 2 of Gautschi.

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> N,  ω² = 40, 0.999
(40, 0.999)

julia> w(t) = ((1-ω²*t^2)*(1-t^2))^(-1/2)
w (generic function with 1 method)

julia> πs = ChebyshevTT{Float64}
ChebyshevTT{Float64}

julia> p = WeightFunction(πs, w, [0,0,1], :x)
WeightFunction(1⋅e_2(x))

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

The integrals in the mixed moment, `σ_{kl} = ∫ p_k⋅π_l dw`,
computation are done through recursion and are slow to compute.


"""
struct WeightFunction{P, F, T <: Number} <: AbstractWeightFunction{T}
    pis::P
    W::F
    coeffs::Vector{T}
    var::Symbol
    # no inner constructor here
end

WeightFunction(pis::Type{<:AbstractOrthogonalPolynomial}, W, coeffs::Vector{T}) where {T} = WeightFunction(pis, W, coeffs, :x)

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

function innerproduct(p::WeightFunction, f, g)
    dom = domain(p)
    fn = x -> f(x) * g(x) * weight_function(p)(x)
    ∫(fn, first(dom)+eps(), last(dom)-eps())
end


## Compute the modified moment
## v_j =  ∫ π_j W(x) dx
## (The moment is ∫ x^j W(x) dx)
## This can  be overridden when an exact calculuation is known, e.g.:
## q0 = WeightFunction(P,f,[1],:x)
## WW = typeof(q0)
## SpecialPolynomials.modified_moment(w::W, j) where {W <: WW} = v(j)
@memoize function modified_moment(w::WeightFunction,  j)

    π_j = convert(Polynomial, Polynomials.basis(w.pis, j))
    π_j /= π_j[end]

    return innerproduct(w, π_j, one)
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

    if k == -1
        return 0
    elseif  iszero(k)
        return modified_moment(w, l)
    elseif iszero(l)

        π_l = 1

        # get basis poly p_l
        zs = zeros(Int, k+1); zs[end]=1
        p_k = WeightFunction(w.pis, w.W, zs, w.var)

        fn = x -> p_k(x) * π_l * weight_function(w.pis)(x)
        return innerproduct(w.pis, p_k, π_l)
    end
    # l+1 or l-1?
    sigma(w, k-1, l+1) - (alpha(w, k-1) - alpha(P, l)) * sigma(w, k-1, l) - (beta(w, k-1))*sigma(w, k-2, l) + beta(P,l)*sigma(w, k-1, l-1)
end

## Compute the α, β functions for the 3-point recurrence for the orthogonal polynomials relative to the weight function
@memoize function alpha(w::WeightFunction, k)
    P = w.pis
    if k == 0
        v0, v1 =  modified_moment(w, 0), modified_moment(w, 1)
        alpha(P, 0) + v1 / v0
    else
        alpha(P, k) - sigma(w, k-1, k)/sigma(w, k-1, k-1) + sigma(w, k, k+1)/sigma(w, k, k)
    end
end

@memoize function beta(w::WeightFunction, k)
    iszero(k) && return innerproduct(w, one, one)
    (sigma(w, k, k)/sigma(w, k-1, k-1))
end
