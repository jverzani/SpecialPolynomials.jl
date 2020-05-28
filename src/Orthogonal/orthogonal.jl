## Abstract  types for  orthogonal  polynomials

## Has An(P), Bn(P), Cn(P)
abstract type AbstractOrthogonalPolynomial{T} <: AbstractSpecialPolynomial{T} end
abstract type AbstractContinuousOrthogonalPolynomial{T} <: AbstractOrthogonalPolynomial{T} end
abstract type AbstractDiscreteOrthogonalPolynomial{T} <: AbstractOrthogonalPolynomial{T} end

abstract type AbstractCOP{T,N} <: AbstractOrthogonalPolynomial{T} end

"""
    AbstractCCOP{T,N}

Following [Koepf  and Schmersau](https://arxiv.org/pdf/math/9703217.pdf), a family `y(x)=p_n(x)=k_x⋅x^n +  ...`  
for  `n  ∈  {0, 1,…}, k_n ≠ 0` of polynomials is a family of classic *continuous* orthogonal polynomials if it  is  a
solution of a differential equation

(a⋅x²+b⋅x+c) ⋅ y'' + (d⋅x + e) ⋅ y' + λᵢ⋅ y = 0.

A family is characterized by the 5 coefficients: a,b,c,d,e.
Let σ = (a⋅x²+b⋅x+c), τ = (d⋅x + e).

From these  5  coefficients several structural  equations are represented. For example
the three-point recusion.

P₍ᵢ₊₁) = (Aᵢ⋅x + Bᵢ) * Pᵢ - Cᵢ *  P₍ᵢ₋₁₎

Here `Aᵢ,Bᵢ,Cᵢ`  can be represented in formulas involving just  `a,b,c,d,e` and `i`.

Rearraging   gives the structural equation

x⋅p_n   = [an, bn, cn]    ⋅ [p_{n+1}, p_n, p_{n-1}]     #  Eqn (7)


The other structural equations are (equation  references are for Koepf  and  Schmerrsaus):

σ⋅p'_n  = [αn, βn, γn]    ⋅  [p_{n+1}, p_n, p_{n-1}]    # Eqn (9), n ≥ 1
p_n    = [ân, b̂n, ĉn]    ⋅  [p'_{n+1}, p'_n, p'_{n-1}] # Eqn (19)
x⋅p'_n  = [αᴵn, βᴵn, γᴵn] ⋅  [p'_{n+1}, p'_n, p'_{n-1}] # Eqn  (14) with  α^*, β^*,  γ^* 

Using (7), Clenshaw polynomial evaluation using the three  point recursion is defined.

Using (19), expressions for derivatives are found.

Using  (19),  expressions   for  integration are  found  (p7).

Using Thms 2,4, and 5, connection coefficients,  C(n,m) satisfying 
P_n(x) =  ∑  C(n,m)  Q_m(x) (n ≥ 0, 0 ≤  m ≤ n) are  found. These 
allow  fallback  definitions for `convert(Polynomial,p)`,  `convert(P, p::Polynomial)`,
`convert(P{α…}, p::P(β…))` and through composition  `p*q`

Subtypes of `AbstractCCOP` are  created through  the  `@registerN` macros, where `N` is the number  of  parameters used to describe the family.

If non-monic versions are desired, then the  leading  term can be  specified through   `kn()`.  

This is  sufficient for many cases, though the general  equations may need specializations when algebraic cancellation is required. 

## Example

For this example, the value of `Bn` at `0` needs help:

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> const SP=SpecialPolynomials
SpecialPolynomials

julia> SP.@register0 MonicLegendre SP.AbstractCCOP0

julia> SP.abcde(::Type{<:MonicLegendre})  = (-1,0,1,-2,0)

julia> SP.Bn(P::Type{<:MonicLegendre}, ::Val{0}) =  0

julia> 𝐐  =  Rational{Int}
Rational{Int64}

julia> x = variable(Polynomial{𝐐})
Polynomial(x)

julia> [basis(MonicLegendre{𝐐}, i)(x) for i  in 0:5]
6-element Array{Polynomial{Rational{Int64}},1}:
 Polynomial(1//1)
 Polynomial(x)
 Polynomial(-1//3 + x^2)
 Polynomial(-3//5*x + x^3)
 Polynomial(3//35 - 6//7*x^2 + x^4)
 Polynomial(5//21*x - 10//9*x^3 + x^5)
```

[Koekoek and Swarttouw](https://arxiv.org/pdf/math/9602214.pdf)
present an encyclopedia of formula characterizing families of
orthogonal polynomials.


"""
abstract type AbstractCCOP{T,N} <: AbstractCOP{T,N} end


"""
     AbstractCDOP{T,N}

Following [Koepf  and Schmersau](https://arxiv.org/pdf/math/9703217.pdf), a family `y(x)=p_n(x)=k_x⋅x^n +  ...`  
for  `n  ∈  {0, 1,…}, k_n ≠ 0` of polynomials is a family of classic *discrete* orthogonal polynomials if it  is  a
solution of a differential equation

(a⋅x²+b⋅x+c) ⋅ Δ∇y + (d⋅x + e) ⋅ ∇' + λᵢ⋅ y = 0,

where  `Δy(x) = y(x+1) - y(x)` and `∇y(x) = y(x) - y(x-1)`.

A family is characterized by the 5 coefficients: a,b,c,d,e.
Let σ = (a⋅x²+b⋅x+c), τ = (d⋅x + e).

As in the classical coninuous orthogonal polynomial case
[`AbstractCCOP`](@ref), from these 5 values the cofficients in the
there-point recursion, and other structural equations can be
represented. These allow polynomial multiplication, integration,
differentiation, conversion, etc. to be defined generically.

[Koekoek and Swarttouw](https://arxiv.org/pdf/math/9602214.pdf)
present an encyclopedia of formula characterizing families of
orthogonal polynomials. 

For example, on p29 they give  formula for Hahn polynomials through:

`n(n+α+β+1)y(x) = B(x)y(x+1) -[B(x)+D(x)]y(x) + D(x)y(x-1)`,  with  explicit values  for  `B` and `D`. Reexpressing gives:
`BΔy(x) - D∇y(x) -λ y(x)  = 0`. From the rexpressed Eqn (4) for Koepf & Schemersau we have the identification:
`σ+τ =  B; σ=D`,  so  `τ=B-D`. From this `a,b,c,d,e` can be  gleaned.


"""
abstract type AbstractCDOP{T,N} <: AbstractCOP{T,N} end

##
## --------------------------------------------------
##


Polynomials.variable(P::Type{<:AbstractOrthogonalPolynomial},  var::Polynomials.SymbolLike=:x) =
    (basis(P,1,var) - Bn(P,0)) / An(P,0)

## Evaluation

# from type, cs, x
function clenshaw_eval(P::Type{<:AbstractOrthogonalPolynomial{T}}, cs, x::S) where {T,S}

    N = length(cs)
    N == 0 && return zero(T)*zero(S)
    N == 1 && return cs[1] * one(S)

    Δ0 = cs[end - 1]
    Δ1 = cs[end]
    @inbounds for i in N-1:-1:2
        Δ0, Δ1 = cs[i - 1] - Δ1 * Cn(P, i-1), Δ0 + Δ1 * muladd(x, An(P,i-1),Bn(P,i-1))
    end

    return Δ0 + Δ1 * muladd(x, An(P,0),  Bn(P,0))
end

# from instance, x
function clenshaw_eval(p::AbstractOrthogonalPolynomial, x::S) where {S}

    T, cs = eltype(p), coeffs(p)
    N = length(cs)
    N == 0 && return zero(T)*zero(S)
    N == 1 && return cs[1] * one(S)

    Δ0 = cs[end - 1]
    Δ1 = cs[end]
    @inbounds for i in N-1:-1:2
        @show cs[i-1]
        Δ0, Δ1 = cs[i - 1] - Δ1 * Cn(p, i-1), Δ0 + Δ1 * muladd(x, An(p,i-1),Bn(p,i-1))
    end

    return Δ0 + Δ1 * muladd(x, An(p,0),  Bn(p,0))
end


## Connection/Linearization

#  Structs for connection, linearization (connection.jl)
struct Connection{P,Q}
    n::Int
    k::Int
end

struct Linearization{P,V}
    l::Int
    n::Int
    m::Int
end


## Collect facts

"""
    weight_function(p)
    weight_function(::Type{P})

For an orthogonal polynomial type, a function `w` with `∫ B_n(t) B_m(t) w(t) dt = 0` when n and m are not equal.

"""
weight_function(::Type{P}) where {P <: AbstractOrthogonalPolynomial} = throw(MethodError("Not implemented"))
weight_function(::P) where {P <: AbstractOrthogonalPolynomial} = weight_function(P)

"""
    generating_function(p)
    generating_function(::Type{P})

The generating function is a function defined by: `(t,x) -> sum(t^n Pn(x) for n in 0:oo)`.
"""
generating_function(::Type{P}) where {P <: AbstractOrthogonalPolynomial} = throw(MethodError("Not implemented"))
generating_function(::P) where {P <: AbstractOrthogonalPolynomial} = generating_function(P)

ismonic(::Type{P}) where {P <: AbstractOrthogonalPolynomial} = false

# cf. https://en.wikipedia.org/wiki/Orthogonal_polynomials#Recurrence_relation
# Orthogonal polynomials have a three-point recursion formula
# parameterized here through:
# P_{n+1} = (An*x + Bn) * P_n + Cn P_{n-1}
# also commonly  written as
# An⋅x⋅P_n = -P_{n+1} + Bn⋅P_n + C_n P_{n-1}
"""
    An(::Type{P},n)
    An(p::P, n)


Orthogonal polynomials defined by a weight function satisfy a three point recursion formula of the form:

`P_{n+1} = (A_n x + B_n) P_{n} - C_n P_{n-1}`

If the polynomials are monic, this is usually parameterized as:

`π_{n+1} = (x - α̃_n) π_n - β̃_n π_{n-1}`

These functions are used through recursion when evaluating the polynomials, converting to `Polynomial` format, for constructing the Vandermonde matrix, for construction the Jacobi matrix, and elsewhere.
"""
An(::Type{P}, n) where {P <: AbstractOrthogonalPolynomial} = throw(MethodError())

"""
    Bn(::Type{P},n)
    Bn(p::P, n)

cf. [`An()`](@ref)
"""
Bn(::Type{P}, n) where {P <: AbstractOrthogonalPolynomial} = throw(MethodError())


"""
    Cn(::Type{P},n)
    Cn(p::P, n)

cf. [`An()`](@ref)
"""
Cn(::Type{P}, n) where {P <: AbstractOrthogonalPolynomial} = throw(MethodError())
An(p::P, n) where {P <: AbstractOrthogonalPolynomial} = An(P,n)
Bn(p::P, n) where {P <: AbstractOrthogonalPolynomial} = Bn(P,n)
Cn(p::P, n) where {P <: AbstractOrthogonalPolynomial} = Cn(P,n)


## For monic polynomials, we have
## π_{n+1} = (x - α̃(n)) π_{n} - β̃(n)π_{n-1}
"""
    π̃αn(::Type{P}, n)

cf. [`An()`](@ref)
"""
π̃αn(P::Type{<:AbstractOrthogonalPolynomial}, n)        = - Bn(P,n) / An(P,n)

"""
    β̃n(::Type{P}, n)

cf. [`An()`](@ref)
"""
function π̃βn(P::Type{<:AbstractOrthogonalPolynomial}, n)
    iszero(n) &&  return innerproduct(P, one, one)
    Cn(P,n)/An(P,n)/An(P, n-1)
end


"""
    leading_term(::Type{P},n)

Return leading term of `basis(P,n)` in the  standard basis. By default this is generated through the three-point recursion.
"""
function leading_term(::Type{P}, n::Int) where {P <: AbstractOrthogonalPolynomial}
    n < 0 && throw(ArgumentError("n must be a non-negative integer"))
    n == 0 && return one(eltype(P))
    prod(An(P,i) for i in n-1:-1:0)
end

"""
    monic(p::AbstractOrthogonalPolynomial)

Return `p` as a monic polynomial *when* represented in the standard basis. Retursn the zero polynomial if the degree of `p` is `-1`. 
"""
function monic(p::P) where {P <: AbstractOrthogonalPolynomial}
    n = degree(p)
    n == -1 && return ⟒(P)(0/one(eltype(p)))
    p / (p[end]*leading_term(P, n))
end


# return the domain as a tuple, not an interval object
function Base.extrema(::Type{P}) where {P <: AbstractOrthogonalPolynomial}
    dom = domain(P)
    first(dom), last(dom)
end
Base.extrema(p::P) where {P <: AbstractOrthogonalPolynomial} = extrema(P)



##
## Vandermonde matrix can be generated through the 3-point recursion formula
##
function Polynomials.vander(p::Type{P}, x::AbstractVector{T}, n::Integer) where
    {P <:  AbstractOrthogonalPolynomial,
     T <: Number}
    
    A = Matrix{T}(undef, length(x), n + 1)
    A[:, 1] .= P0(P, one(T))

    if n > 0
        A[:, 2] .= P1(P, x)
        @inbounds for i in 1:n-1
            A[:, i+2] .= A[:, i+1] .* (An(p, i)*x .+ Bn(p,i)) .+ (Cn(p,i) * A[:, i])
        end
    end
    
    return A
    
end

# Jacobi matrix
# roots(basis(P,n)) = eigvals(jacobi_matrix(P,n)), but this is more stable
# https://en.wikipedia.org/wiki/Gaussian_quadrature#The_Golub-Welsch_algorithm
"""
    jacobi_matrix(::Type{P}, n)
    jacobi_matrix(p::P, n)

The Jacobi Matrix is a symmetric tri-diagonal matrix. The diagonal entries are the `alpha_i` values, the off diagonal entries,
the square root of the  `beta_i` values. This matrix has the properties that

* the eigenvalues are the roots of the corresponding basis vector. As these roots are important in quadrature, and finding eigenvalues of qa symmetric tri-diagonal matrix yields less error than finding the eigenvalues of the companion matrix, this can be used for higher degree basis polynomials.
* the normalized eigenvectors have initial term proportional to the weights in a quadrature formula

See the [FastGaussQuadrature](https://github.com/JuliaApproximation/FastGaussQuadrature.jl) package for faster implementations.

"""
function jacobi_matrix(::Type{P}, n) where {P <: AbstractOrthogonalPolynomial}
    LinearAlgebra.SymTridiagonal([π̃αn(P,i) for i in 0:n-1], [sqrt(π̃βn(P,i)) for i in 1:n-1])
end
jacobi_matrix(p::P, n) where {P <: AbstractOrthogonalPolynomial} = jacobi_matrix(P,n)

##  Compute weights and nodes for quadrature
"""
    gauss_nodes_weights(::Type{P}, n)
    gauss_nodes_weights(p::P, n)

Returns a tuple of nodes and weights for Gauss quadrature for the given orthogonal type.

For some types, a method from  A. Glaser, X. Liu, and V. Rokhlin. "A fast algorithm for the calculation of the roots of special functions." SIAM J. Sci. Comput., 29 (2007), 1420-1438. is used. 

For others the Jacobi matrix, J_n, for which the Golub-Welsch] algorithm The nodes  are computed from the eigenvalues of J_n, the weights a scaling of the first component of the normalized eigen vectors (β_0 * [v[1] for v in vs])

!!! note
    See the [FastGaussQuadrature](https://github.com/JuliaApproximation/FastGaussQuadrature.jl) package for faster, vastly more engineered implementations.

"""
function gauss_nodes_weights(p::Type{P}, n) where {P <: AbstractOrthogonalPolynomial}
    J = jacobi_matrix(P, n)
    eig = eigen(J, extrema(P)...)
    # Is this necessary?
    nm  = 1 #diag(eig.vectors * eig.vectors')
    wts =  π̃βn(P,0) * (eig.vectors[1,:] ./ nm).^2
    eig.values,  wts
end
gauss_nodes_weights(p::P, n) where {P <: AbstractOrthogonalPolynomial} = gauss_nodes_weights(P, n)

# Trait to indicate if computation of nodes and weights is is O(n) or O(n^2)
has_fast_gauss_nodes_weights(p::P)  where {P <: AbstractOrthogonalPolynomial} = has_fast_gauss_nodes_weights(P)
has_fast_gauss_nodes_weights(::Type{P})  where {P <: AbstractOrthogonalPolynomial} = false


##
## --------------------------------------------------
##

## <f,g> = ∫ f⋅g⋅w dx
function innerproduct(P::Type{<:Union{AbstractCCOP,AbstractContinuousOrthogonalPolynomial}}, f, g)
    dom = domain(P)
    a, b = first(dom), last(dom)
    if !first(Polynomials.inclusivity(dom))
        a += eps(float(one(a)))
    end
    if !last(Polynomials.inclusivity(dom))
        b -= eps(float(one(b)))
    end
    fn = x -> f(x) * g(x) * weight_function(P)(x)
    
    return _quadgk(fn, a, b)
end

innerproduct(p::P, f, g) where {P <: AbstractContinuousOrthogonalPolynomial} =  innerproduct(P, f, g)

## Compute <p_i, p_i> = \| p \|^2; allows export, but work is in norm2
Base.abs2(::Type{P}, n) where {P <: AbstractOrthogonalPolynomial} = norm2(P, n)
Base.abs2(p::P, n) where {P <: AbstractOrthogonalPolynomial} = norm2(p, n)

## Compute <p_i, p_i> = \| p \|^2
## Slow default; generally should  be directly expressed for each type
function norm2(::Type{P}, n) where {P <: AbstractOrthogonalPolynomial}
    p = basis(P,n)
    innerproduct(P, p, p)
end
norm2(p::P, n) where {P <: AbstractOrthogonalPolynomial} = norm2(P, n)


"""
    fit([Val(S)], P::Type{<:AbstractOrthogonalPolynomial}, f, n::Int; var=:x)

    Find an approximating polynomial of degree `n` or less for a function `f`, that is returns `p(x) = ∑ᵢⁿ cᵢ Pᵢ(x)` for some coefficients `cᵢ`.

Defaults to an interpolating polynomial. To specify others, use one of `Val(:interpolating)`, `Val(:lsq)` (least squares), or `Val(:series)` (trunated series expansion) as the first argument. See [`SpecialPolynomials.cks`](@ref) for some more detail.

"""
Polynomials.fit(P::Type{<:AbstractOrthogonalPolynomial}, f, n::Int; var=:x) =
    fit(Val(:interpolating), P, f, n, var=var)

"""
    fit(val::Val{:interpolating}, P::Type{<:AbstractOrthogonalPolynomial}, f, deg::Int; var=:x)

Fit `f` with an interpolating polynomial of degree `n` orless using nodes
`x0,x1, ..., xn`, the  zeros of `P_{n+1} = basis(P, n+1)`. and `p(xᵢ)=f(xᵢ)`.
"""
Polynomials.fit(val::Val{:interpolating},
                P::Type{<:AbstractOrthogonalPolynomial}, f, n::Int;
                var=:x) =
                    P(cks(val, P,f,n), var)

"""
    fit(Val(:lsq), P::Type{<:AbstractOrthogonalPolynomial}, f, n::Int)

Fit `f` with `p(x)=∑ d_i P_i(x)` where `p` had degree `n` or less using least squares.
"""
Polynomials.fit(val::Val{:lsq},
                P::Type{<:AbstractOrthogonalPolynomial}, f, n::Int;
                var=:x) =
                    P(cks(val, P,f,n), var)

"""
    fit(Val(:series), P::Type{<:AbstractOrthogonalPolynomial}, f, n::Int)

If `f(x)` is written as an infinite sum `∑ c_kP_k(x)`, this returns a truncated sum `∑_0^n c̃_k P_k(x)` where `n` ischosen algorithmically and `c̃_k`is chosen using an efficient manner, not necessarily through the orthogonality condition, `<f,p_k>/<p_k,p_k>`.

"""
Polynomials.fit(val::Val{:series},
                P::Type{<:AbstractOrthogonalPolynomial}, f;
                var=:x,kwargs...) = 
                    P(cks(val, P,f; kwargs...), var)


"""
    cks(::Val{:interpolating}, ::Type{P}, f, n::Int)

Fit `f` with the interpolating polynomial using roots of P_{n+1}.

Let xs,ws be the gauss nodes and weights of P (xs are the zeros of P_{n+1}).

Then if we interpolate `f` at the `xs` to get `p`, then
`p(xᵢ) = f(xᵢ)` and 

`∑ p_n(xᵢ) p_m(xᵢ) wᵢ = K_m δ_{nm}` and

`p(x) = ∑ᵢ cᵢ Pᵢ(x)`.

Using this:
`∑ᵢ f(xᵢ)P_k(xᵢ) wᵢ =` 
`∑ᵢ (∑_j c_j P_j(xᵢ)) P_k(xᵢ) wᵢ =`
`∑_j c_j (K_k δ_{ik}) = c_k K_k`, So 
`c_k = (1/K_k) ∑ᵢ f(xᵢ)P_k(xᵢ) wᵢ`
"""
function cks(::Val{:interpolating}, ::Type{P}, f, n::Int) where {P  <:  AbstractOrthogonalPolynomial}
    
    xs, ws = gauss_nodes_weights(P, n)
    return [sum(f(xⱼ) * basis(P, k)(xⱼ) * wⱼ for (xⱼ,wⱼ) in zip(xs,ws)) / norm2(P,k) for k in 0:n]
end


"""
    cks(::Val{:lsq}, ::Type{P}, f, n::Int)

Fit `f` with a polynomial `∑ᵢⁿ cᵢ Pᵢ` chosen so `<f-p,f-p>_w` is as small as possible.
 Using the normal equations, the coefficients are found to be `c_k = <f,P_k>_w / <P_k,P_k>_w.` 
For some types an approximation to the inner product, `<f,P_k>_w` may be used.


ref: http://www.math.niu.edu/~dattab/MATH435.2013/APPROXIMATION
"""
function cks(::Val{:lsq}, ::Type{P}, f, n::Int) where {P  <:  AbstractOrthogonalPolynomial}
## return ck =  <f,P_k>/<P_k,P_k>, k =  0...n
    [innerproduct(P, f, basis(P,k)) /norm2(P,k) for k in 0:n]
end

"""
    cks(::Val{:series}, ::Type{P}, f, n::Int)

If `f(x)` is written as an infinite sum `∑ c_kP_k(x)`, then this 
tries to identify an `n` for which the series expansion is a good approximation and returns the coefficients. 

"""
function cks(::Val{:series}, ::Type{P}, f, n::Int) where {P  <:  AbstractOrthogonalPolynomial}
    throw(MethodError())
end

