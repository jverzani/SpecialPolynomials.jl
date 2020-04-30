## References
## Gautschi: https://www.cs.purdue.edu/homes/wxg/Madrid.pdf
##

"""
    AbstractOrthogonalPolynomial{T}

Abstract type for an orthogonal polynomial family.

An orthogonal polynomial family, `P_0`, `P_1`, ..., `P_n`, ... has the properties that:

* the degree  of `P_i` is `i`
* the set `P_0`, `P_1`, ..., `P_n` forms a basis for the polynomials of degree `n` or less
* the polynomials are orthogonal with respect to  a weight function over a domain. That is if `i ≠ j` then
`∫_a^b P_i ⋅ P_j ⋅ w dt = 0`

The fact that these form a basis, allows the representation a degree
`n` polynomial through a coefficient vector with `n+1` elements. That is
`p = a_0⋅P_0 + a_1⋅P_1 + ⋅⋅⋅ + a_n⋅P_n` can be represented by the vector 
`[a_0, a_1, ..., a_n]` in the same manner that `p = b_0 + b_1⋅x + ⋅⋅⋅ + b_n⋅x^n`
is represented by `[b_0, b_1, ..., b_n]`. The coefficients depend on the 
underlying family, which is implicit in the choice of defining type.

The polynomials have numerous properties. For purposes of definitions,
the three-point recursion property allows the implementation of much
of the functionality, including evaluation and conversion to the
standard basis (the `Polynomial` type).

"""
abstract type AbstractOrthogonalPolynomial{T} <: AbstractSpecialPolynomial{T} end
abstract type OrthogonalPolynomial{T} <: AbstractOrthogonalPolynomial{T} end
abstract type DiscreteOrthogonalPolynomial{T} <: AbstractOrthogonalPolynomial{T} end

basis_symbol(::Type{<:AbstractOrthogonalPolynomial}) = "P"
function Polynomials.showterm(io::IO, ::Type{P}, pj::T, var, j, first::Bool, mimetype) where {N, T, P <: OrthogonalPolynomial}
    iszero(pj) && return false
    !first &&  print(io, " ")
    print(io, Polynomials.hasneg(T)  && Polynomials.isneg(pj) ? "- " :  (!first ? "+ " : ""))
    print(io, "$(abs(pj))⋅$(basis_symbol(P))_$j($var)")
    return true
end


"""
    weight_function(p)
    weight_function(::Type{P})

For an orthogonal polynomial family, a function `w` with `integrate B_n(t) B_m(t) w(t) dt = 0` when n and m are not equal.

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

# cf. https://en.wikipedia.org/wiki/Orthogonal_polynomials#Recurrence_relation
# Orthogonal polynomials have a three-point recursion formula
# parameterized here through:
# P_{n+1} = (An*x + Bn) * P_n + Cn P_{n-1}
"""
    An(::Type{P},n)
    An(p::P, n)


Families of orthogonal polynomials defined by a weight function satisfy a three point recursion formula of the form:

`P_{n+1} = (A_n x + B_n) P_{n} + C_n P_{n-1}`

If the polynomials are monic, this is usually parameterized as:

`π_{n+1} = (x - alpha_n) π_n - beta_n π_{n-1}`

These functions are used through recursion when evaluating the polynomials, converting to `Polynomial` format, for constructing the Vandermonde matrix, for construction the Jacobi matrix, and elsewhere.
"""
An(::Type{P}, n) where {P <: AbstractOrthogonalPolynomial} = throw(MethodError("Not implemented"))

"""
    Bn(::Type{P},n)
    Bn(p::P, n)

cf. [`An()`](@ref)
"""
Bn(::Type{P}, n) where {P <: AbstractOrthogonalPolynomial} = throw(MethodError("Not implemented"))


"""
    Cn(::Type{P},n)
    Cn(p::P, n)

cf. [`An()`](@ref)
"""
Cn(::Type{P}, n) where {P <: AbstractOrthogonalPolynomial} = throw(MethodError("Not implemented"))
An(p::P, n) where {P <: AbstractOrthogonalPolynomial} = An(P,n)
Bn(p::P, n) where {P <: AbstractOrthogonalPolynomial} = Bn(P,n)
Cn(p::P, n) where {P <: AbstractOrthogonalPolynomial} = Cn(P,n)


## For monic polynomials, we have
## π_{n+1} = (x - α(n)) π_{n} - β(n)π_{n-1}
"""
    alpha(::Type{P}, n)
    alpha(p::P, n)

cf. [`An()`](@ref)
"""
alpha(P::Type{<:AbstractOrthogonalPolynomial}, n)        = - Bn(P,n) / An(P,n)
alpha(p::P, n) where {P <: AbstractOrthogonalPolynomial} = alpha(P, n)

"""
    beta(::Type{P}, n)
    beta(p::P, n)

cf. [`An()`](@ref)
"""
function beta(P::Type{<:AbstractOrthogonalPolynomial}, n)
    iszero(n) &&  return innerproduct(P, one, one)
    -Cn(P,n)/An(P,n)/An(P, n-1)
end
beta(p::P, n) where {P <: AbstractOrthogonalPolynomial} = beta(P, n)

# typically, P_0 = 1; P_{-1} = 0 and P_1 from recursion, but here we compute P1
P0(::Type{P}, x) where {P <: AbstractOrthogonalPolynomial} = one(x)
P0(p::P, x) where {P <: AbstractOrthogonalPolynomial} = P0(P,x)
P1(::Type{P}, x) where {P <: AbstractOrthogonalPolynomial} = An(P,0)*x .+ Bn(P,0)
P1(p::P, x) where {P <: AbstractOrthogonalPolynomial} = P1(P,x)
dP0(::Type{P},  x) where {P <: AbstractOrthogonalPolynomial} = zero(x) 

Polynomials.variable(::Type{P}, var::Polynomials.SymbolLike=:x) where {P <: AbstractOrthogonalPolynomial} = P([-Bn(P,0), 1]/An(P,0))
Polynomials.variable(p::P, var::Polynomials.SymbolLike=:x) where {P <: AbstractOrthogonalPolynomial} = variable(P)


function Base.extrema(::Type{P}) where {P <: AbstractOrthogonalPolynomial}
    dom = domain(P)
    first(dom), last(dom)
end
Base.extrema(p::P) where {P <: AbstractOrthogonalPolynomial} = extrema(P)


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
    LinearAlgebra.SymTridiagonal([alpha(P,i) for i in 0:n-1], [sqrt(beta(P,i)) for i in 1:n-1])
end
jacobi_matrix(p::P, n) where {P <: AbstractOrthogonalPolynomial} = jacobi_matrix(P,n)

##  Compute weights and nodes for quadrature
"""
    gauss_nodes_weights(::Type{P}, n)
    gauss_nodes_weights(p::P, n)

Returns a tuple of nodes and weights for Gauss quadrature for the given orthogonal family.

For some families, a method from  A. Glaser, X. Liu, and V. Rokhlin. "A fast algorithm for the calculation of the roots of special functions." SIAM J. Sci. Comput., 29 (2007), 1420-1438. is used. 

For others the Jacobi matrix, J_n, for which the Golub-Welsch] algorithm The nodes  are computed from the eigenvalues of J_n, the weights a scaling of the first component of the normalized eigen vectors (β_0 * [v[1] for v in vs])

!!! note
    See the [FastGaussQuadrature](https://github.com/JuliaApproximation/FastGaussQuadrature.jl) package for faster, vastly more engineered implementations.

"""
function gauss_nodes_weights(p::Type{P}, n) where {P <: AbstractOrthogonalPolynomial}
    J = jacobi_matrix(P, n)
    eig = eigen(J, extrema(P)...)
    # Is this necessary?
    nm  = 1 #diag(eig.vectors * eig.vectors')
    wts =  beta(P,0) * (eig.vectors[1,:] ./ nm).^2
    eig.values,  wts
end
gauss_nodes_weights(p::P, n) where {P <: AbstractOrthogonalPolynomial} = gauss_nodes_weights(P, n)

has_fast_gauss_nodes_weights(p::P)  where {P <: AbstractOrthogonalPolynomial} = has_fast_gauss_nodes_weights(P)
has_fast_gauss_nodes_weights(::Type{P})  where {P <: AbstractOrthogonalPolynomial} = false



## Evaluate an orthogonal polynomial using the three-point recursion representation and
## Clenshaw recurrence.
## For each type define (ch::Type)(x) = orthogonal_polyval(ch, x)
## as we can't call methods defined for an abstract type
function orthogonal_polyval(ch::P, x::S) where {P <: AbstractOrthogonalPolynomial, S}
    S <: Number && !(first(domain(ch)) <= x <= last(domain(ch))) && throw(ArgumentError("$x outside of domain"))
    T = eltype(ch)
    oS = one(x)
    R = typeof(one(T) * oS)
    length(ch) == 0 && return zero(R)
    length(ch) == 1 && return R(ch[0])


    c0 = ch[end - 1]
    c1 = ch[end]
    @inbounds for i in degree(ch):-1:2
        c0, c1 = ch[i - 2] + c1 * Cn(ch, i-1), c0 + c1 * muladd(x, An(ch,i-1),Bn(ch,i-1))
    end
    # Now have c0 P_0(x) + c1 P_1 = c0 + c1*x
    return c0 * P0(ch, x) + c1 * P1(ch, x)
end

function monic_orthogonal_basis_polyval(P::Type{<:OrthogonalPolynomial}, j,  x::S) where {S}

    j == 0 && return one(S)
    j == 1 && return x - alpha(P,0)

    oS = one(S)
    c0 = zero(S)
    c1 = oS
    β = α = alpha(P,0)
    p0 = oS
    p1 = (x - α)*p0
    @inbounds for i in j:-1:2
        α, β = alpha(P,i-1), beta(P, i-1)
        c0, c1 = c1 * (-β), c0 + c1 * (x - α)
    end
    # Now have c0 P_0(x) + c1 P_1 = c0 + c1*x
    return c0 * p0 + c1 * p1
end

## Evaluate an orthogonal polynomial and its  derivative using the three-point recursion representation
function orthogonal_polyval_derivative(p::P, x::S) where {P <: AbstractOrthogonalPolynomial, S}

    T = eltype(p)
    oS = one(x)
    R = eltype(one(T) * oS)
    
    length(p) == 0 && return (zero(R), zero(R))
    d = length(p) - 1

    # P0,  P_{-1} case
    pn,  pn_1 = P0(P,x), zero(x)
    dpn, dpn_1  = dP0(P,x), zero(x)

    ptot = p[0]*P0(P,x)
    dptot  = p[0]*dP0(P,x)
    
    
    ## we compute  forward here, finding pi, dp_i,  p_{i+1}, dp_{i+1}, ...
    for i in 0:d-1
        an, bn, cn=  An(p,i), Bn(p,i), Cn(p,i)
        dpn_1, dpn = dpn,  an*pn + (an*x + bn)*dpn + cn*dpn_1
        pn_1,  pn =  pn, (an*x + bn)*pn + cn*pn_1
        ptot += p[i+1]*pn
        dptot += p[i+1]*dpn
    end
    return (ptot, dptot)
end

## sn: (x,n) -> ? is used to scale the polyonmials
## We march forward here, not backwards as is done with Clenshaw recursion
# function orthogonal_polyval_derivative(p::P, x::S,
#                                        scale=(x,n)-> (one(x),one(x),zero(x),zero(x))
#                                        ) where {P <: AbstractOrthogonalPolynomial, S}
#     T = eltype(p)
#     oS = one(x)
#     R = eltype(one(T) * oS)
    
#     length(p) == 0 && return (zero(R), zero(R))
#     d = length(p) - 1

    
#     λ1, λ0, dλ1, dλ0 = scale(x, 0)
#     pn_1 = zero(R)  # π_{-1}
#     pn = one(R) * λ0 # π_0
#     dpn_1 = zero(R) # dπ_{-1}
#     dpn = dλ0

#     ptot, dptot =  p[0]*pn,  p[0]*dpn


#     ## we compute  forward here, finding pi, dp_i,  p_{i+1}, dp_{i+1}, ...
#     for i in 0:d-1
#         λ1, λ0, dλ1, dλ0 = scale(x, i)
#         an, bn, cn=  An(p,i), Bn(p,i), Cn(p,i)
#         dpn_1, dpn = dpn,  (an*λ1 +  (an*x+bn)*dλ1)*pn  + cn*dλ0*pn_1 + ((an*x+bn)*λ1)*dpn + cn*λ0*dpn_1
#         pn_1,  pn =  pn, (an*x + bn)*λ1*pn + cn*λ0*pn_1
#         ptot += p[i+1]*pn
#         dptot += p[i+1]*dpn
#     end

#     return (ptot, dptot)
# end







# conversion to Polynomial is simply evaluating
# the polynomial on `variable(Polynomial)`
function Base.convert(P::Type{<:Polynomial}, ch::AbstractOrthogonalPolynomial)
    if length(ch) == 1
        return P(ch.coeffs, ch.var)
    end
    T = eltype(ch)
    ##
    ## P_{n+1} = (An x + Bn)P_n + CnP_{n-1}
    ## so ... c0 P_{n-1} + c1 P_n
    ## c1 = c0 + (An x+ Bn)
    ## c0 = p[i-1] + Cn
    c0 = P(ch[end - 1], ch.var)
    c1 = P(ch[end], ch.var)
    x = variable(P)
    Q= typeof(ch)
    @inbounds for i in degree(ch):-1:2
        c0, c1 = ch[i - 2] + c1 * Cn(ch, i-1), c0 + c1 * (An(ch, i-1) * x + Bn(ch, i-1))
    end
    return c0 *  P0(ch,  x) + c1 * P1(ch, x)
end

# brute force this, by matching up leading terms and subtracting
function Base.convert(P::Type{J}, p::Polynomial) where {J <: AbstractOrthogonalPolynomial}
    d = degree(p)
    R = eltype(one(eltype(p))/1)
    pp =  convert(Polynomial{R}, p)
    qs = zeros(R, d+1)
    for i in d:-1:0
        Pn = convert(Polynomial{R}, basis(P,  i))
        qs[i+1] = lambda =  pp[i]/Pn[i]
        pp = pp - lambda*Pn
    end
    P(qs, p.var)
end



function Polynomials.vander(p::Type{P}, x::AbstractVector{T}, n::Integer) where {P <:  AbstractOrthogonalPolynomial, T <: Number}
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








# <f,g> = ∫ f⋅g⋅w dx
function innerproduct(P::Type{<: OrthogonalPolynomial}, f, g)
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

innerproduct(p::P, f, g) where {P <: AbstractOrthogonalPolynomial} =  innerproduct(P, f, g)

## Compute <p_i, p_i> = \| p \|^2; allows export, but work is in norm2
Base.abs2(::Type{P}, n) where {P <: AbstractOrthogonalPolynomial} = norm2(P, n)
Base.abs2(p::P, n) where {P <: AbstractOrthogonalPolynomial} = norm2(p, n)

## Compute <p_i, p_i> = \| p \|^2
## Slow default; generally should  be directly expressed for each family
function norm2(::Type{P}, n) where {P <: AbstractOrthogonalPolynomial}
    p = basis(P,n)
    innerproduct(P, p, p)
end
norm2(p::P, n) where {P <: AbstractOrthogonalPolynomial} = norm2(P, n)

# compute ∫x^k dw
moment(P::Type{<:AbstractOrthogonalPolynomial}, k::Int) = innerproduct(p, x -> x^k, one)
moment(p::AbstractOrthogonalPolynomial, k::Int) = innerproduct(p, x->x^k, one)

function dot(p::P, q::P) where {P <: AbstractOrthogonalPolynomial}
    # use orthogonality  here

    n, m = degree(p), degree(q)
    tot = zero(eltype(P))
    for i in 0:minimum(n,m)
        tot += p[i] * q[i] * norm2(P, i)
    end

    tot
end
dot(::Type{P}, f, i::Int) where {P <: AbstractOrthogonalPolynomial} = innerproduct(P, f, Polynomials.basis(P, i))
dot(::Type{P}, f, p::P) where {P <: AbstractOrthogonalPolynomial} = innerproduct(P, f, p)

## return ck =  <f,P_k>/<P_k,P_k>, k =  0...n
## innerproduct is likely  slow here in default
## Some types (ChebyshevT) have a faster alternative
function cks(::Type{P}, f, n::Int) where {P  <:  AbstractOrthogonalPolynomial}

    return [innerproduct(P, f, basis(P, k))/norm2(P,k) for k in 0:n]
end



## for an orthogonal family, fit f with a0 p_0 + ... + a_d P_d
## where a_i = <f, p_i>/<p_i,p_i>
function Polynomials.fit(P::Type{<:AbstractOrthogonalPolynomial}, f, deg::Int; var=:x)
    # use least squares fit of orthogonal polynomial family to f
    P(cks(P,f,deg), var)
end


## Some utilities

_quadgk(f, a, b) = quadgk(Wrapper(f), a, b)[1]
const ∫ = _quadgk
_monic(p::AbstractOrthogonalPolynomial) = p/convert(Polynomial,p)[end]
