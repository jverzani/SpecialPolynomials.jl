## References
## Gautschi: https://www.cs.purdue.edu/homes/wxg/Madrid.pdf
##

"""
    OrthogonalPolynomial{T}

Abstract type for an orthogonal polynomial family.

An orthogonal polynomial family, `P_0`, `P_1`, ..., `P_n`, ... has the properties that:

* the degree  of `P_i` is `i`
* the set `P_0`, `P_1`, ..., `P_n` forms a basis for the polynomials of degree `n` or less
* the polynomials are orthogonal with respect to  a weight function over a domain. That is if `i ≠ j` then
`∫_a^b P_i ⋅ P_j ⋅ w dt = 0`

The fact that these form a basis, allows the representation a degree `n` polynomial through a coefficient vector with `n+1` elements.

The polynomials have numerous properties. For purposes of definitions, the three-point recursion property allows the implementation of much of the functionality, including evaluation and conversion to the standard basis (the `Polynomial` type).
"""
abstract type OrthogonalPolynomial{T} <: AbstractSpecialPolynomial{T} end

basis_symbol(::Type{<:AbstractSpecialPolynomial}) = "P"
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
weight_function(::Type{P}) where {P <: OrthogonalPolynomial} = throw(MethodError("Not implemented"))
weight_function(::P) where {P <: OrthogonalPolynomial} = weight_function(P)

"""
    generating_function(p)
    generating_function(::Type{P})

The generating function is a function defined by: `(t,x) -> sum(t^n Pn(x) for n in 0:oo)`.
"""
generating_function(::Type{P}) where {P <: OrthogonalPolynomial} = throw(MethodError("Not implemented"))
generating_function(::P) where {P <: OrthogonalPolynomial} = generating_function(P)

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
An(::Type{P}, n) where {P <: OrthogonalPolynomial} = throw(MethodError("Not implemented"))

"""
    Bn(::Type{P},n)
    Bn(p::P, n)

cf. [`An()`](@ref)
"""
Bn(::Type{P}, n) where {P <: OrthogonalPolynomial} = throw(MethodError("Not implemented"))


"""
    Cn(::Type{P},n)
    Cn(p::P, n)

cf. [`An()`](@ref)
"""
Cn(::Type{P}, n) where {P <: OrthogonalPolynomial} = throw(MethodError("Not implemented"))
An(p::P, n) where {P <: OrthogonalPolynomial} = An(P,n)
Bn(p::P, n) where {P <: OrthogonalPolynomial} = Bn(P,n)
Cn(p::P, n) where {P <: OrthogonalPolynomial} = Cn(P,n)


## For monic polynomials, we have
## π_{n+1} = (x - α(n)) π_{n} - β(n)π_{n-1}
"""
    alpha(::Type{P}, n)
    alpha(p::P, n)

cf. [`An()`](@ref)
"""
alpha(P::Type{<:OrthogonalPolynomial}, n)        = - Bn(P,n) / An(P,n)
alpha(p::P, n) where {P <: OrthogonalPolynomial} = alpha(P, n)

"""
    beta(::Type{P}, n)
    beta(p::P, n)

cf. [`An()`](@ref)
"""
function beta(P::Type{<:OrthogonalPolynomial}, n)
    iszero(n) &&  return _quadgk(weight_function(P), extrema(P)...)
    -Cn(P,n)/An(P,n)/An(P, n-1)
end
beta(p::P, n) where {P <: OrthogonalPolynomial} = beta(P, n)

# typically, P_0 = 1; P_{-1} = 0 and P_1 from recursion, but here we compute P1
P0(::Type{P}, x) where {P <: OrthogonalPolynomial} = one(x)
P0(p::P, x) where {P <: OrthogonalPolynomial} = P0(P,x)
P1(::Type{P}, x) where {P <: OrthogonalPolynomial} = An(P,0)*x .+ Bn(P,0)
P1(p::P, x) where {P <: OrthogonalPolynomial} = P1(P,x)


Polynomials.variable(::Type{P}, var::Polynomials.SymbolLike=:x) where {P <: OrthogonalPolynomial} = P([-Bn(P,0), 1]/An(P,0))
Polynomials.variable(p::P, var::Polynomials.SymbolLike=:x) where {P <: OrthogonalPolynomial} = variable(P)


function Base.extrema(::Type{P}) where {P <: OrthogonalPolynomial}
    dom = domain(P)
    first(dom), last(dom)
end
Base.extrema(p::P) where {P <: OrthogonalPolynomial} = extrema(P)


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

"""
function jacobi_matrix(::Type{P}, n) where {P <: OrthogonalPolynomial}
    LinearAlgebra.SymTridiagonal([alpha(P,i) for i in 0:n-1], [sqrt(beta(P,i)) for i in 1:n-1])
end
jacobi_matrix(p::P, n) where {P <: OrthogonalPolynomial} = jacobi_matrix(P,n)

##  Compute weights and nodes for quadrature
"""
    gauss_nodes_weights(::Type{P}, n)
    gauss_nodes_weights(p::P, n)

Returns a tuple of nodes and weights for Gauss quadrature for the given orthogonal family.

Uses the Jacobi matrix, J_n, for which [Golub Welsch](https://en.wikipedia.org/wiki/Gaussian_quadrature#The_Golub-Welsch_algorithm)
yields:

* nodes -- roots of the orthogonal polynomial π_{n} computed from the eigenvalues of J_n
* weights -- if vs are the normalized eigen vectors of J_n, then  β_0 * [v[1] for v in vs]

"""
function gauss_nodes_weights(p::Type{P}, n) where {P <: OrthogonalPolynomial}
    J = jacobi_matrix(P, n)
    eig = eigen(J, extrema(P)...)
    # Is this necessary?
    nm  = 1 #diag(eig.vectors * eig.vectors')
    wts =  beta(P,0) * (eig.vectors[1,:] ./ nm).^2
    eig.values,  wts
end
gauss_nodes_weights(p::P, n) where {P <: OrthogonalPolynomial} = gauss_nodes_weights(P, n)




## Evaluate an orthogonal polynomial using the three-point recursion representation
## Need to specify for each type (ch::Type)(x) = orthogonal_polyval(ch, x)
## as we can't have abstract type  here
function orthogonal_polyval(ch::P, x::S) where {P <: OrthogonalPolynomial, S}
    S <: Number && x ∉ domain(ch) && throw(ArgumentError("$x outside of domain"))
    T = eltype(ch)
    oS = one(x)
    R = eltype(one(T) * oS)
    #R = promote_type(T, S)
    length(ch) == 0 && return zero(R)
    length(ch) == 1 && return ch[0]


    c0 = ch[end - 1]
    c1 = ch[end]
    @inbounds for i in degree(ch):-1:2
        c0, c1 = ch[i - 2] + c1 * Cn(ch, i-1), c0 + c1 * (An(ch, i-1) * x + Bn(ch, i-1))
    end
    # Now have c0 P_0(x) + c1 P_1 = c0 + c1*x
    return c0 * P0(ch, x) + c1 * P1(ch, x)
end

# conversion to Polynomial is simply evaluating
# the polynomial on `variable(Polynomial)`
function Base.convert(P::Type{<:Polynomial}, ch::OrthogonalPolynomial)
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
function Base.convert(P::Type{J}, p::Polynomial) where {J <: OrthogonalPolynomial}
    d = degree(p)
    R = eltype(one(eltype(p))/1)
    pp =  convert(Polynomial{R}, p)
    qs = zeros(R, d+1)
    for i in d:-1:0
        Pn = convert(Polynomial{R}, Polynomials.basis(P,  i))
        qs[i+1] = lambda =  pp[i]/Pn[i]
        pp = pp - lambda*Pn
    end
    P(qs, p.var)
end



function Polynomials.vander(p::Type{P}, x::AbstractVector{T}, n::Integer) where {P <:  OrthogonalPolynomial, T <: Number}
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
    fn = x -> f(x) * g(x) * weight_function(P)(x)
    a,err =  quadgk(fn, first(dom)+sqrt(eps()), last(dom)-sqrt(eps()))
    a
end

## Compute <p_i, p_i> = \| p \|^2; allows export, but work is in norm2
Base.abs2(::Type{P}, n) where {P <: OrthogonalPolynomial} = norm2(P, n)
Base.abs2(p::P, n) where {P <: OrthogonalPolynomial} = norm2(p, n)

## Compute <p_i, p_i> = \| p \|^2
## Slow default; generally should  be directly expressed for each family
function norm2(::Type{P}, n) where {P <: OrthogonalPolynomial}
    p = Polynomials.basis(P,n)
    innerproduct(P, p, p)
end
norm2(p::P, n) where {P <: OrthogonalPolynomial} = norm2(P, n)

function dot(p::P, q::P) where {P <: OrthogonalPolynomial}
    # use orthogonality  here

    n, m = degree(p), degree(q)
    tot = zero(eltype(P))
    for i in 0:minimum(n,m)
        tot += p[i] * q[i] * norm2(P, i)
    end

    tot
end
dot(::Type{P}, f, i::Int) where {P <: OrthogonalPolynomial} = innerproduct(P, f, Polynomials.basis(P, i))
dot(::Type{P}, f, p::P) where {P <: OrthogonalPolynomial} = innerproduct(P, f, p)

## return ck =  <f,P_k>/<P_k,P_k>, k =  0...n
function cks(::Type{P}, f, n::Int) where {P  <:  OrthogonalPolynomial}

    xs =  eigvals(jacobi_matrix(P,n+1))
    fxs  = f.(xs)

    p0 = P0.(P, xs)
    p1 = P1.(P, xs)

    cks = similar(xs)
    cks[1] = dot(fxs, p0) / norm2(P,0) / pi
    cks[2] = dot(fxs, p1) / norm2(P,1) / pi


    for i in 1:n-1
        p0[:], p1[:] = p1, p1 .* (An(P, i)*xs .+ Bn(P,i)) .+ (Cn(P,i) * p0)
        cks[i+2] =  dot(fxs, p1) / norm2(P,i) / pi
    end
    cks / (n+1)
end



## for an orthogonal family, fit f with a0 p_0 + ... + a_d P_d
## where a_i = <f, p_i>/<p_i,p_i>
function Polynomials.fit(P::Type{<:OrthogonalPolynomial}, f, deg::Int; var=:x)
    # use least squares fit of orthogonal polynomial family to f
    P(cks(P,f,deg), var)
end


## Some utilities
_quadgk(f, a, b) = quadgk(f, a+eps(), b-eps())[1]
_monic(p::OrthogonalPolynomial) = p/convert(Polynomial,p)[end]
