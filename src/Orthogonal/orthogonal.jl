## Abstract  types for  orthogonal  polynomials

abstract type  AbstractOrthogonalBasis <: AbstractSpecialPolynomialBasis end

## Has An(P), Bn(P), Cn(P)
## These are now Basis; have to adjust for non basis polys (weightfunction?)
#abstract type AbstractOrthogonalPolynomial{B,T,X} <: AbstractSpecialPolynomial{B,T,X} end
#abstract type AbstractContinuousOrthogonalPolynomial{B,T,X} <:
#              AbstractOrthogonalPolynomial{B,T,X} end
#abstract type AbstractDiscreteOrthogonalPolynomial{B,T,X} <: AbstractOrthogonalPolynomial{B,T,X} end

"""
    AbstractOrthogonalPolynomial{T,X}

This is an alias for polys with an Orthogonal Basis (`AbstractOrthogonalBasis`) specified.

These polynomials have  several properties, including an accompanying inner product satsifying  `⟨yᵢ, yⱼ⟩ = cᵢδᵢⱼ`.

In addition to methods inherited from the underlying `AbstractPolynomial`  type, orthogonal polynomial  types may have methods   `weight_function`, `generating_function`, `leading_term`, `norm2`, `jacobi_matrix`, and `gauss_nodes_weights`,  though none are  exported.


Subtypes of `AbstractCOPBasis <: AbstractOrthogonalBasis` utilize the fact that the basis  polynomials  satisfy

`(ax² + bx + c)yᵢ''(x) + (dx+e)yᵢ'(x) + λᵢyᵢ(x) = 0` (or a discrete analogue)

where the structural relations are functions of `a,b,c,d,e`. These allow default definitions for polynomial evaluation,   addition, multiplication, differentiation, integration, and  conversion to and from  the `Polynomial` type (the `FallingFactorial` type in the discrete  c case),

A key structural relation is the three-term recursion,  `yᵢ₊₁ =  (Aᵢx +  Bᵢ)yᵢ -  Cᵢyᵢ₋₁`. For systems  specfied by  a  weight function, the  values of `Aᵢ`, `Bᵢ`, and `Cᵢ` can  be  generated, yielding formulas for polynomial evaluation, addition, and conversion to the `Polynomial`  type throughe evaluation.
"""
const AbstractOrthogonalPolynomial  = AbstractUnivariatePolynomial{<:AbstractOrthogonalBasis,T,X} where {T,X}


##
## --------------------------------------------------
##
function Polynomials.fromroots(
    P::Type{<:AbstractOrthogonalPolynomial},
    roots::AbstractVector;
    var::Polynomials.SymbolLike=:x,
)
    x = variable(P, var)
    p = prod(x - r for r in roots)
end

Base.convert(P::Type{<:AbstractOrthogonalPolynomial}, c::Number) = c * one(P)

function Base.one(::Type{P})  where {P<:AbstractOrthogonalPolynomial}
    B = Polynomials.basistype(P)
    basis(P, 0) / k0(B)
end

function Polynomials.variable(::Type{P})  where {P<:AbstractOrthogonalPolynomial}
    B = Polynomials.basistype(P)
    (basis(P, 1) / k0(B) - Bn(B, 0)) / An(B, 0)
end

##
## --------------------------------------------------
##

function Polynomials.derivative(p::P) where {P<:AbstractOrthogonalPolynomial}
    T = eltype(one(P))
    q = convert(Polynomial{T}, p)
    convert(⟒(P), derivative(q))
end

function Polynomials.integrate(p::P) where {P<:AbstractOrthogonalPolynomial}
    q = convert(Polynomial, p)
    integrate(q)
end

# XXX How to use ngcd here, as there are round off errors leading to nonsenxe
function Base.divrem(num::P, den::P) where {P<:AbstractOrthogonalPolynomial}
    p1 = convert(Polynomial, num)
    p2 = convert(Polynomial, den)
    q, r = divrem(promote(p1, p2)...)
    convert.(P, (q, r))
end

function Base.gcd(p1::P, p2::Q;
                  atol::Real=zero(real(T)),
                  rtol::Real=Base.rtoldefault(real(T)),
                  method=:numerical, # lossy conversion makes this desirable
                  kwargs...
                  ) where {T, X, P<:AbstractOrthogonalPolynomial{T,X},
                           S,    Q<:AbstractOrthogonalPolynomial{S,X}}
    u, v = convert.(Polynomial, (p1, p2))
    gcd(u, v; atol=atol, rtol=rtol, method=method, kwargs...)
end

function Polynomials.companion(p::P) where {P<:AbstractOrthogonalPolynomial}
    companion(convert(Polynomial, p))
end

function Polynomials.vander(
    ::Type{P},
    x::AbstractVector{T},
    n::Integer,
) where {B<:AbstractOrthogonalBasis, P<:AbstractUnivariatePolynomial{B},T}
    N = length(x) - 1
    R = typeof(one(T) / one(T))
    V = zeros(R, N + 1, n + 1)

    for j in 0:n
        bj = Polynomials.basis(P, j)
        for i in 0:N
            V[i + 1, j + 1] = bj(x[i + 1])
        end
    end

    V
end

## Evaluation

# from type, cs, x
function clenshaw_eval(::Type{B}, cs, x::S) where {B<:AbstractOrthogonalBasis,S}
    T, N = eltype(cs), length(cs)
    p₀ = k0(B)
    R = promote_type(typeof(p₀), promote_type(promote_type(T, S), typeof(An(B, 0))))
    N == 0 && return zero(R)
    N == 1 && return (cs[1] * p₀) * one(R)

    Δ0::R = cs[end - 1]
    Δ1::R = cs[end]
    @inbounds for i in (N - 1):-1:2
        Aᵢ,Bᵢ,Cᵢ = ABCₙ(B, i-1)
        Δ0, Δ1 =
            cs[i - 1] - Δ1 * Cᵢ, Δ0 + Δ1 * muladd(x, Aᵢ, Bᵢ)
    end
    A₀,B₀,_ = ABCₙ(B, 0)
    p₁ = muladd(x, A₀, B₀) * p₀
    return Δ0 * p₀ + Δ1 * p₁
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

## Collect facts  about  the orthogonal polynomial system

"""
    weight_function(p)
    weight_function(::Type{P})

For an orthogonal polynomial type, a function `w` with `∫ B_n(t) B_m(t) w(t) dt = 0` when `n` and `m` are not equal.

"""
weight_function(::Type{B}) where {B<:AbstractOrthogonalBasis}# = throw(ErrorException("Not implemented"))
weight_function(::P) where {B,P<:AbstractOrthogonalPolynomial{B}} = weight_function(B)

"""
    generating_function(p)
    generating_function(::Type{P})

The generating function is a function defined by: `(t,x) -> sum(t^n Pn(x) for n in 0:oo)`.
"""
generating_function(::Type{B}) where {B<:AbstractOrthogonalBasis} #  = throw(ArgumentError("Not implemented"))
generating_function(::P) where {B,P<:AbstractOrthogonalPolynomial{B}} = generating_function(B)

"""
    leading_term(::Type{P},n)

Return leading term of `basis(P,n)` in the  standard basis. By default this is generated through the three-point recursion.
"""
function leading_term(::Type{B}, n::Int) where {B<:AbstractOrthogonalBasis}
    n < 0 && throw(ArgumentError("n must be a non-negative integer"))
    n == 0 && return 1
    prod(An(B, i) for i in (n - 1):-1:0)
end

# is P a monic polynomial system?
ismonic(::Type{B}) where {B<:AbstractOrthogonalBasis} = false
ismonic(::Type{P}) where {B<:AbstractOrthogonalBasis, P<:AbstractUnivariatePolynomial{B}} = ismonic(B)
isorthonormal(::Type{B}) where {B<:AbstractOrthogonalBasis} = false
isorthonormal(::Type{P}) where {B<:AbstractOrthogonalBasis, P<:AbstractUnivariatePolynomial{B}} = isorthonormal(B)

# cf. https://en.wikipedia.org/wiki/Orthogonal_polynomials#Recurrence_relation
# Orthogonal polynomials have a three-point recursion formula
# parameterized here through:
# P_{n+1} = (An*x + Bn) * P_n + Cn P_{n-1}
# also commonly  written as
# An⋅x⋅P_n = -P_{n+1} + Bn⋅P_n + C_n P_{n-1}
"""
    An(::Type{P},n)


Orthogonal polynomials defined by a weight function satisfy a three point recursion formula of the form:

`P_{n+1} = (A_n x + B_n) P_{n} - C_n P_{n-1}`

If the polynomials are monic, this is usually parameterized as:

`π_{n+1} = (x - α̃_n) π_n - β̃_n π_{n-1}`

These functions are used through recursion when evaluating the polynomials, converting to `Polynomial` format, for constructing the Vandermonde matrix, for construction the Jacobi matrix, and elsewhere.
"""
An(::Type{B}, n) where {B<:AbstractOrthogonalBasis} # = throw(ArgumentError("No default method"))

"""
    Bn(::Type{B},n)
    Bn(p::P, n)

cf. [`An()`](@ref)
"""
Bn(::Type{B}, n) where {B<:AbstractOrthogonalBasis} # = throw(ArgumentError("No default method"))

"""
    Cn(::Type{B},n)
    Cn(p::P, n)

cf. [`An()`](@ref)
"""
Cn(::Type{B}, n) where {B<:AbstractOrthogonalBasis} # = throw(ArgumentError("No default method"))

## For monic polynomials, we have
## π_{n+1} = (x - α̃(n)) π_{n} - β̃(n)π_{n-1}
"""
    π̃αn(::Type{P}, n)

cf. [`An()`](@ref)
"""
π̃αn(B::Type{<:AbstractOrthogonalBasis}, n) = -Bn(B, n) / An(B, n)

"""
    β̃n(::Type{B}, n)

cf. [`An()`](@ref)
"""
function π̃βn(B::Type{<:AbstractOrthogonalBasis}, n)
    iszero(n) && return innerproduct(B, one, one)
    Cn(B, n) / An(B, n) / An(B, n - 1)
end

"""
    monic(p::AbstractOrthogonalPolynomial)

Return `p` as a monic polynomial *when* represented in the standard basis. Returns the zero polynomial if the degree of `p` is `-1`.
"""
function monic(p::P) where {B<:AbstractOrthogonalBasis, P<:AbstractOrthogonalPolynomial{B}}
    n = degree(p)
    n == -1 && return ⟒(P)(0 / one(eltype(p)))
    p / (p[end] * leading_term(B, n))
end

##
## Work with compact  basis
##
"""
     Basis(P,n)
     Basis{P}(n)

The  command `basis(P,n, [var])` realizes the polynomial. `Basis(P,n)` does  not.  This  can  be  useful for  evaluation.
"""
struct Basis{Π}
    n::Int
end
Basis(B, n) = Basis{B}(n)
Basis(::Type{P}, n) where {B<:AbstractBasis,P<:AbstractUnivariatePolynomial{B}} = Basis(B, n)
(b::Basis{B})(x) where {B} = basis(MutableDensePolynomial{B},b.n)(x)
basistype(b::Basis{B}) where {B} = B
Base.show(io::IO, mimetype::MIME"text/plain", b::Basis{B}) where {B} =
    print(io, "$(B)($(b.n))")

function innerproduct(::B, f::Basis{B}, g::Basis{B}) where {B<:AbstractOrthogonalBasis}
    n, m = f.n, g.n
    if n == m
        return norm2(B, n)
    else
        return zero(B)
    end
end
innerproduct(::P, f::Basis{B}, g::Basis{B}) where {B <: AbstractOrthogonalBasis, P<:AbstractUnivariatePolynomial{B}} =
    innerproduct(B, f, g)
##
## Vandermonde matrix can be generated through the 3-point recursion formula
##
function Polynomials.vander(
    p::Type{P},
    x::AbstractVector{T},
    n::Integer,
) where {B<:AbstractOrthogonalBasis, P<:AbstractUnivariatePolynomial{B},T<:Number}
    A = Matrix{T}(undef, length(x), n + 1)

    # for i in 0:n
    #     A[:, i+1] = basis(P, i).(x)
    # end
    # return A

    A[:, 1] .= basis(B, 0)(one(T)) #P0(P, one(T))

    if n > 0
        A[:, 2] .= basis(B, 1).(x) #P1(P, x)
        @inbounds for i in 1:(n - 1)
            # `P_{n+1} = (A_n x + B_n) P_{n} - C_n P_{n-1}`
            n′ = i + 1
            A[:, n′ + 1] =
                (An(B, n′) * x .+ Bn(B, n′)) .* A[:, n′] .- (Cn(B, n′) * A[:, n′ - 1])
        end
    end

    return A
end

# Jacobi matrix
# roots(basis(P,n)) = eigvals(jacobi_matrix(P,n)), but this is more stable
# https://en.wikipedia.org/wiki/Gaussian_quadrature#The_Golub-Welsch_algorithm
"""
    jacobi_matrix(::Type{P}, n)

The Jacobi Matrix is a symmetric tri-diagonal matrix. The diagonal entries are the `alphaᵢ` values, the off diagonal entries,
the square root of the  `betaᵢ` values. This matrix has the properties that

* the eigenvalues are the roots of the corresponding basis vector. As these roots are important in quadrature, and finding eigenvalues of a symmetric tri-diagonal matrix yields less error than finding the eigenvalues of the companion matrix, this can be used for higher degree basis polynomials.
* the normalized eigenvectors have initial term proportional to the weights in a quadrature formula

"""
function jacobi_matrix(::Type{B}, n::Int) where {B<:AbstractOrthogonalBasis}
    a, b = [π̃αn(B, i) for i in 0:(n - 1)], [sqrt(π̃βn(B, i)) for i in 1:(n - 1)]
    LinearAlgebra.SymTridiagonal(promote(a,b)...)
end
jacobi_matrix(::Type{P}, n::Int) where {B, P<:AbstractUnivariatePolynomial{B}} = jacobi_matrix(B, n)
jacobi_matrix(p::P, n::Int) where {B, P<:AbstractUnivariatePolynomial{B}} = jacobi_matrix(B, n)

##  Compute weights and nodes for quadrature
"""
    gauss_nodes_weights(::Type{P}, n)

Returns a tuple of nodes and weights for Gauss quadrature for the given orthogonal type.

When loaded, the values are computed through  the `FastGaussQuadrature` package.

For some types, a method from  A. Glaser, X. Liu, and V. Rokhlin. "A fast algorithm for the calculation of the roots of special functions." SIAM J. Sci. Comput., 29 (2007), 1420-1438. is used.

For others the Jacobi matrix, J_n, for which the Golub-Welsch] algorithm The nodes  are computed from the eigenvalues of J_n, the weights a scaling of the first component of the normalized eigen vectors (β_0 * [v[1] for v in vs])

"""
function gauss_nodes_weights(::Type{B}, n) where {B<:AbstractOrthogonalBasis}#, P<:AbstractUnOrthogonalPolynomial{B}}
    J = jacobi_matrix(B, n)
    eig = eigen(J, extrema(domain(B))...)
    # Is this necessary?
    nm  = 1 #diag(eig.vectors * eig.vectors')
    wts = π̃βn(B, 0) * (eig.vectors[1, :] ./ nm) .^ 2
    eig.values, wts
end
gauss_nodes_weights(b::Basis{B}) where {B} = gauss_nodes_weights(B, b.n)
gauss_nodes_weights(p::P,n) where {B<:AbstractOrthogonalBasis, P<:AbstractUnivariatePolynomial{B}} = gauss_nodes_weights(B, n)




##
## --------------------------------------------------
##

"""
      innerproduct(::Type{P}, f, g)

Compute  `<f,g> = ∫ f⋅g⋅w dx` where  `w` is the weight function of the  type  `P`  and the integral is  taken  over  the domain of the type `P`.
"""
function innerproduct(
    B::Type{<:AbstractOrthogonalBasis},
    f,
    g;
    atol=sqrt(eps()),
)
    dom = domain(B)
    a, b = first(dom), last(dom)
    if first(bounds_types(dom)) == Open
        a += eps(float(one(a)))
    end
    if last(bounds_types(dom)) == Open
        b -= eps(float(one(b)))
    end
    fn = x -> f(x) * g(x) * weight_function(B)(x)

    return quadgk(fn, a, b; atol=atol)[1]
end

## Compute <p_i, p_i> = \| p \|^2; allows export, but work is in norm2
Base.abs2(::Type{P}, n) where {B<:AbstractOrthogonalBasis, P<:AbstractUnivariatePolynomial{B}} = norm2(B, n)

## Compute <p_i, p_i> = \| p \|^2
## Slow default; generally should  be directly expressed for each type
function norm2(::Type{B}, n) where {B<:AbstractOrthogonalBasis}
    p = basis(MutableDensePolynomial{B}, n)
    innerproduct(B, p, p)
end

##
## --------------------------------------------------
##

"""
    fit([Val(S)], P::Type{<:AbstractOrthogonalPolynomial}, f, n::Int; var=:x)

Find an approximating polynomial of degree `n` or less for a function `f`, that is returns `p(x) = ∑ᵢⁿ cᵢ Pᵢ(x)` for some coefficients `cᵢ`.


Defaults to an interpolating polynomial. To specify others, use one of `Val(:interpolating)`, `Val(:lsq)` (least squares), or `Val(:series)` (trunated series expansion) as the first argument. See [`SpecialPolynomials.cks`](@ref) for some more detail.

"""
Polynomials.fit(P::Type{<:AbstractOrthogonalPolynomial}, f, n::Int; var=:x) =
    fit(Val(:interpolating), P, f, n, var=var)

"""
    fit(val::Val{:interpolating}, P::Type{<:AbstractOrthogonalPolynomial}, f, deg::Int; var=:x)

Fit `f` with an interpolating polynomial of degree `n` or less using nodes
`x0,x1, ..., xn`, the  zeros of `P_{n+1} = basis(P, n+1)`. and `p(xᵢ)=f(xᵢ)`.
"""
Polynomials.fit(
    val::Val{:interpolating},
    P::Type{<:AbstractOrthogonalPolynomial},
    f,
    n::Int;
    var=:x,
) = P(cks(val, P, f, n), var)

"""
    fit(Val(:lsq), P::Type{<:AbstractOrthogonalPolynomial}, f, n::Int)

Fit `f` with `p(x)=∑ d_i P_i(x)` where `p` had degree `n` or less using least squares.
"""
Polynomials.fit(
    val::Val{:lsq},
    P::Type{<:AbstractOrthogonalPolynomial},
    f,
    n::Int;
    var=:x,
) = P(cks(val, P, f, n), var)

"""
    fit(Val(:series), P::Type{<:AbstractOrthogonalPolynomial}, f, n::Int)

If `f(x)` is written as an infinite sum `∑ c_kP_k(x)`, this returns a truncated sum `∑_0^n c̃_k P_k(x)` where `n` is chosen algorithmically and `c̃_k` is chosen using an efficient manner, not necessarily through the orthogonality condition, `<f,p_k>/<p_k,p_k>`.

"""
Polynomials.fit(
    val::Val{:series},
    P::Type{<:AbstractOrthogonalPolynomial},
    f;
    var=:x
) = P(cks(val, P, f), var)

"""
    cks(::Val{:interpolating}, ::Type{P}, f, n::Int)

Fit `f` with the interpolating polynomial using roots of `P_{n+1}`.

Let `xs`, `ws` be the gauss nodes and weights of `P` (`xs` are the zeros of `P_{n+1}`).

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
function cks(
    ::Val{:interpolating},
    ::Type{P},
    f,
    n::Int,
) where {B<:AbstractOrthogonalBasis, P<:AbstractUnivariatePolynomial{B}}
    xs, ws = gauss_nodes_weights(B, n)
    return [
        sum(f(xⱼ) * basis(P, k)(xⱼ) * wⱼ for (xⱼ, wⱼ) in zip(xs, ws)) / norm2(B, k) for
        k in 0:n
    ]
end

"""
    cks(::Val{:lsq}, ::Type{P}, f, n::Int)

Fit `f` with a polynomial `∑ᵢⁿ cᵢ Pᵢ` chosen so `<f-p,f-p>_w` is as small as possible.
 Using the normal equations, the coefficients are found to be `c_k = <f,P_k>_w / <P_k,P_k>_w.`
For some types an approximation to the inner product, `<f,P_k>_w` may be used.


ref: [http://www.math.niu.edu/~dattab/MATH435.2013/APPROXIMATION](http://www.math.niu.edu/~dattab/MATH435.2013/APPROXIMATION)
"""
function cks(::Val{:lsq}, ::Type{P}, f, n::Int) where {P<:AbstractOrthogonalPolynomial}
    ## return ck =  <f,P_k>/<P_k,P_k>, k =  0...n
    [innerproduct(P, f, basis(P, k)) / norm2(P, k) for k in 0:n]
end

"""
    cks(::Val{:series}, ::Type{P}, f, n::Int)

If `f(x)` is written as an infinite sum `∑ c_kP_k(x)`, then this
tries to identify an `n` for which the series expansion is a good approximation and returns the coefficients.

"""
cks(::Val{:series}, ::Type{P}, f, n::Int) where {P<:AbstractOrthogonalPolynomial} #   throw(ArgumentError("No default method"))
