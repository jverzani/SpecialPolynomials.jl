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
# also commonly  written as
# An⋅x⋅P_n = -P_{n+1} + Bn⋅P_n + C_n P_{n-1}
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

"""
    leading_term(::Type{P},n)

Return leading term of `basis(P,n)` in the  standard basis. By default this is generated through the three-point recursion.
"""
function leading_term(::Type{P},n::Int) where {P <: AbstractOrthogonalPolynomial}
    n < 0 && throw(ArgumentError("n must be a non-negative integer"))
    n == 0 && return P0(P,0)
    prod(An(P,i) for i in n-1:-1:0)
end

"""
    monic(p::AbstractOrthogonalPolynomial)

Return `p` as a monic polynomial *when* represented in the standard basis. Retursn the zero polynomial if the degree of `p` is `-1`. 
"""
function monic(p::P) where {P <: AbstractOrthogonalPolynomial}
    n = degree(p)
    n == -1 && return ⟒(P)(0/one(eltype(p)))
    p / leading_term(P, n)
end
    
# return `x`  in the  polynomial basis.
Polynomials.variable(::Type{P}, var::Polynomials.SymbolLike=:x) where {P <: AbstractOrthogonalPolynomial} = P([-Bn(P,0), 1]/An(P,0))
Polynomials.variable(p::P, var::Polynomials.SymbolLike=:x) where {P <: AbstractOrthogonalPolynomial} = variable(P)

# return the domain as a tuple, not an interval object
function Base.extrema(::Type{P}) where {P <: AbstractOrthogonalPolynomial}
    dom = domain(P)
    first(dom), last(dom)
end
Base.extrema(p::P) where {P <: AbstractOrthogonalPolynomial} = extrema(P)


##
## Evaluation
##
## We have three evaluations:
## * p(x) (Clenshaw recurrence  using  An, Bn, Cn)
## * π(x) (Clenshshaw  recurrence  using alpha and beta)
## * (p(x), p'(x)) using recurrence
##

## Evaluate an orthogonal polynomial using the three-point recursion representation and
## Clenshaw recurrence.
## For each type define (ch::Type)(x) = orthogonal_polyval(ch, x)
## as we can't call methods defined for an abstract type
function orthogonal_polyval(ch::P, x::S) where {P <: AbstractOrthogonalPolynomial, S}
    S <: Real && !(first(domain(ch)) <= x <= last(domain(ch))) && throw(ArgumentError("$x outside of domain"))
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

        
## Direct implementations
##             Linearization   Inversion   Connection{αs} integral derivative
## Legendre         ✓              ✓          NA             ✓        ✓         
## Hermite          ✓              ✓          NA             ✓        ✓
## Hermite_e        ✓              ✓          NA             ✓        ✓
## ChebyshevT       ✓              ✓          T <-> U        ✓        ✓
## ChebyshevU       ✓              ✓          T <-> U        ✓        ✓
## Laguerre        p29             ✓          ✓              ✓        ✓
## Jacobi          p29             ✓          ✓              ✓        .
## Gegenbauer       ✓              ✓          ✓              ✓        ✓
## Bessel           X              X          ✓               X        .



##
## conversion and multiplication
##

## conversion of a polynomial `p = ∑ aᵢϕᵢ` of type P into `q=∑ bᵢψᵢ`
##  of type `Q` can be done through polynomial evaluation
## using the scheme: `x = variable(Q); p(x)`.
##
## *However* this requires being able to pass `x` through the
## polynomial `p` and in particular, being able to multiply elements in
## `Q`. So conversion from Q to P is possible if multiplication in Q
## is defined.  Multiplication can be defined if a "linearization"
## formula is available of the form
##
## ψᵢψⱼ = ∑_k α(k,i,j) ψ_k
##
## Alternatively, conversion from `Q` to `P` can be directly done
## should there be a "connection"
##
## ψᵢ = ∑ α_{i,j} ϕ_j
##
## If `Q=Polynomial`, then conversion can be used to define
## multiplication through `x = variable(Q); convert(P, p(x)*q(x))`
## that is, multiplication is done in the standard basis.

##
## Linearization
##

# Let `P` be a polynomial family with basis `ϕᵢ`. A linearization formula is:
#
# ϕ_n * ϕ_m = ∑₀^(m+n) a(l,n,m) ϕ_l
#
# Let `p` have degree `n` and `q` have degree `m`, both polynomials in
# `P`. The linearization formula leads to a product formula `pq=p⋅q`
# through the following.  The `l`th coefficient, `c`, of the product
# pq is given formally by
#
# c = 0
# for i = 0:n
#   for j = 0:m
#      c += p[i]*q[j]*α(l, i, j)
#   end
# end
#
# To utilize this formula, define a `linearization_α(), l, n, m)` function.
#
# Typically, the matrix (i,j) -> α(d, i,j) is sparse so that the sum
# over `n⋅m` terms can be greatly reduced.  For this we turn a
# linearization formula into an iterator over the non-zero terms
#
# A full iteration over the entire range of n*m -(l-1)*l/2 values of
# (p,q) would look like this:
#
# function Base.iterate(o::Linearization{P}, state =  nothing) where {P <: PolynomialFamily}
#
#     l, m, n = o.l, o.m, o.n
#     l  > m + n && return nothing
#    
#     if state == nothing
#
#         k =  0
#         p = min(l, n)
#         q = l +  k - p
#
#     else
#
#         k, p, q, val  = state
#
#         if p == 0 ||  q == m
#             # bump k
#             k  += 1
#             l + k > n+m && return nothing  # at end  of k
#             p = min(l+k, n)
#             q = l + k - p
#         else
#             p -= 1
#             q += 1
#         end
#
#     end
#    
#     val = linearization_α(P, l, p, q)
#                
#     return (p,q,val), (k, p, q, val)
# end
#
#
# For a sparse example, were we to do this for the standard polynomial
# basis where x^n * x^m = x^(n+m) we would have a(l,m,n) = δ_{l, n+m}
# and we would have an iterator for `Linearization{Polynomial}(l,n,m)`
# defined as
#
# function Base.iterate(o::Linearization{<:Polynomial}, state=nothing)
#
#     if state == nothing
#         l,m,n = o.l, o.m, o.n
#         l > m + n && return nothing
#         p = min(l,n)
#         q = l-p
#         return (p,q,1), (p,q,1)
#     else
#         l,m,n = o.l, o.m, o.n
#         p,q,val =  state
#         p == 0 && return nothing
#         q == m && return nothing
#         p -= 1
#         q += 1
#         return (p,q,1), (p,q,1)
#     end
# end
#
# The the double sum could become `sum(p[i]*q[j]*val for (i,j,val) in
# Linearization{P}(l,n,m))`
#
# To visualize, α(2,4,5) might look like this where the iterator would skip the `*` values as the degree is too small
# and skip the `0` values, as `δ_{l,i,j}=0` is `i+j > l`.
# *  *  1  0  0
# *  1  0  0  0
# 1  0  0  0  0
# 0  0  0  0  0
# 0  0  0  0  0
# 0  0  0  0  0
#
# for the legendre polynomials,
# The  matrix α(2, 4, 5) might look like
#
# *  *  x  0  0
# *  x  0  x  0
# x  0  x  0  x
# 0  x  0  x  0
# 0  0  x  0  x
# 0  0  0  x  0
#
# the striped diagonal terms are a result of the basis polynomials
# being even or odd. The diagonal terms are not all complete, as α can
# be 0 for certain values.
#
#  References.
#  These formula  appear in various places, such as
#
# Daniel Duviol Tcheutia, On Connection, Linearization and DuplicationCoefficients of Classical Orthogonal Polynomials (http://hdl.handle.net/123456789/2014071645714)
# for specific formula for the classic orthogonal polynomials  (Jacobi, Laguerre, Hermite, and Bessel).
#
# Koepf and  Schmersau, Representations of Orthogonal Polynomials (https://arxiv.org/pdf/math/9703217.pdf)
#
# Chaggara and Mabrouk, On inversion and connection coefficients for basic hypergeometric polynomials (https://arxiv.org/pdf/1601.06122.pdf)
#

##
struct Linearization{P,V}
    l::Int
    n::Int
    m::Int
end
Linearization{P}(l,m,n) where {P} = Linearization{P, Val{:diagonal}}(l,m,n)


# non-sparse iteration over p,q (n*m steps)
# to hook into this, define `linearization_α(P, k, n, m)` and an interator along the lines of
# function Base.iterate(o::Linearization{P, Val{:pq}}, state =  nothing) where {α, P <: AbstractOrthogonalPolynomial}
#     k = state == nothing ? 1 : state + 1
#     k > o.n + o.m && return nothing
#     return (k, linearization_α(P, k, n, m)), k
# end
#
#
function linearization_product_pq(p::P,  q::Q) where {P <: AbstractOrthogonalPolynomial, Q <: AbstractOrthogonalPolynomial}
    p.var == q.var || throw(ArgumentError("bases must match"))
    ⟒(P) == ⟒(Q) || throw(ArgumentError("Base polynomial  family must match"))
    PP = promote_type(P,Q)

    n,m = length(p)-1, length(q)-1 # n,m = degree(p), degree(q)
    T,S = eltype(p), eltype(q)
    R = eltype(one(T)*one(S)/1)
    cs = zeros(R, n+m+1)

    ## pn*pm = ∑ a(l,m,n) H_l
    for n in eachindex(p)
        p_n = p[n]
        iszero(p_n) && continue
        for m in eachindex(q)
            c_nm = p_n * q[m]
            for (k, val) in Linearization{P, Val{:pq}}(0, n,m)
                cs[1+k] = muladd(c_nm, val, cs[1+k])
            end
        end
    end
    
    ⟒(P)(cs, p.var)
end

# Sparser product, implented by fixing k and finding combinations of (m,n) that produce a
# non-zero α(k,n,m) for P_k through the linearization formula P_n⋅P_m = ∑ α(k,n,m) P_k
#
# compute product using a linearization formula, as implemented in an `iterate` method for the type
#
@inline function linearization_product(pn::P,  pm::Q) where {P <: AbstractOrthogonalPolynomial, Q <: AbstractOrthogonalPolynomial}
    pn.var == pm.var || throw(ArgumentError("bases must match"))
    ⟒(P) == ⟒(Q) || throw(ArgumentError("Base polynomial  family must match"))
    PP = promote_type(P,Q)

    n,m = length(pn)-1, length(pm)-1 # n,m = degree(pn), degree(pm)
    T,S = eltype(pn), eltype(pm)
    R = eltype(one(T)*one(S)/1)
    cs = zeros(R, n+m+1)
    # pn*pm = ∑ a(l,m,n) H_l
    # so ∑_l ∑_p ∑_q a(l,p,q) H_l
    # can  use iterator to keep as simple sum
    # `sum(pn[p] * pm[q] * val for (p,q,val) in Linearization{PP}(l,n,m))`
    for l in 0:n+m
        for (p,q,val) in Linearization{PP}(l,n,m)
            cs[1+l] += pn[p] * pm[q] * val
        end
    end
    ⟒(P)(cs, pn.var)
end

function Base.:*(p1::P, p2::Q) where
    {P <: AbstractOrthogonalPolynomial,
     Q <: AbstractOrthogonalPolynomial}
    
    R = promote_type(P, Q)

    if hasmethod(iterate, (Linearization{R, Val{:diagonal}}, ))
        linearization_product(p1, p2)
    elseif hasmethod(iterate, (Linearization{R, Val{:pq}}, ))
        linearization_product_pq(p1, p2)        
    else
        # gives poor error message if ⟒(R) not available
        convert(⟒(R), convert(Polynomial, p1) * convert(Polynomial, p2))
    end
    
end

## The connection problem relates
##
##    ψ_n(x) = ∑_{k=0}^n α_{n,k} ϕ_k(x)
##
## That is we convert from basis ψ (Q) into basis family ϕ (P) through
## (`convert(P, q::Q)`):
##
## ∑_{i=0}^n  aᵢ ψᵢ = ∑_{i=0}^n aᵢ ∑_{k=0}^i α_{i,k} ϕ_kⱼ
##                  = ∑_{k=0}^n (∑_{i=k}^n aᵢ * α_{i,k}) ϕ_k
##
## If Q is the standard basis this is the inversion problem.
##
## For this implementation, we suppose there is an is an iterator
## α_{n, k} which iterates
##
## (k, α(k,k)), (k+1, α(k+1, k)), ..., (n, α(n, k))
##
## defined by iterate(Connection{ϕ, ψ}(n, k)
##
## When the ϕ family is `Polynomial`, we can just use evaluation
## through the Clenshaw recursion using `x =variable(Polynomial)`. For
## other conversions, we can directly define a `Base.convert` method,
## or we can define a `Connection` object and have `connection` create
## the coefficients
##
## When ψ_n  = ∑_k^n α_{n,k} ϕ_k
##
## the iterator for (n,k) returns
## (i, val) ∈  { (k, α_{k,k}), (k+1, α_{k+1,k}),  ..., (n, α_{n,k}) }
##
## A sample would be to convert polys in Q to polys in P:
# function Base.iterate(o::Connection{P, Q}, state=nothing) where
#     {P <: Polynomial, Q <: Polynomial}
#
#     n,k = o.n, o.k
#     if state  == nothing
#         j = k
#     else
#         j = state
#         j += 1  # likely j += 2 if parity condition
#     end
#     j > n && return nothing
#     val = connection_α(P,Q,j,k)
#     (j,val), j
# end
struct Connection{P,Q}
    n::Int
    k::Int
end

function connection(::Type{P}, q::Q) where
    {P <: Polynomials.AbstractPolynomial,
     Q<:Polynomials.AbstractPolynomial}
    
    n = degree(q)
    T = eltype(one(q)/1)
    cs = zeros(T, n+1)

    for k in 0:n
        for (i,val) in Connection{P,Q}(n, k)
            cs[1+k] = muladd(q[i], val, cs[1+k])
        end
        #cs[1+k] = sum(q[i] * val for (i,val) in Connection{P,Q}(n, k))
    end

    ⟒(P)(cs, q.var)
    
end
                      



# conversion to Polynomial is simply evaluating
# the polynomial on `variable(Polynomial)`
function Base.convert(P::Type{<:Polynomials.StandardBasisPolynomial},
                      ch::AbstractOrthogonalPolynomial)

    x = variable(P)
    return ch(x)
    
end

function Base.convert(P::Type{J}, q::Q) where
    {J <: AbstractOrthogonalPolynomial,
     Q <: Polynomials.StandardBasisPolynomial}
    
    if  hasmethod(iterate, (Connection{P,Q},))
        return connection(P, q)
    else
        # default conversion from the standard basis polynomial into type P
        # brute force this, by matching up leading terms and subtracting
        d = degree(q)
        R = eltype(one(eltype(q))/1)
        QQ = Polynomials.constructorof(Q){R}
        qq =  convert(QQ, q)
        qs = zeros(R, d+1)
        for i in d:-1:0
            Pn = convert(QQ, basis(P,  i))
            qs[i+1] = lambda =  qq[i]/Pn[i]
            qq -= lambda*Pn
        end
        ⟒(P)(qs, q.var)
    end
    
end

function Base.convert(::Type{P}, q::Q) where
    {P <: AbstractOrthogonalPolynomial,
     Q <: AbstractSpecialPolynomial}

    P == Q && return q
    
    if hasmethod(iterate, (Connection{P,Q},))
        connection(P, q)
    else
        # else, try conversion through Polynomial
        # Does not give good error messages though
        convert(P, convert(Polynomial, q))
    end
    
end


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

##
## Gauss nodes and weights
##


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

# Trait to indicate if computation of nodes and weights is is O(n) or O(n^2)
has_fast_gauss_nodes_weights(p::P)  where {P <: AbstractOrthogonalPolynomial} = has_fast_gauss_nodes_weights(P)
has_fast_gauss_nodes_weights(::Type{P})  where {P <: AbstractOrthogonalPolynomial} = false








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


##
## Fitting
##
##
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
`p(xᵢ) = f(xᵢ)` and `∑ p_n(xᵢ) p_m(xᵢ) wᵢ = K_m δ_{nm}` and
`p(x) = ∑ᵢ cᵢ Pᵢ(x)`.
Using this:
`∑ᵢ f(xᵢ)P_k(xᵢ) wᵢ` `= ∑ᵢ (∑_j c_j P_j(xᵢ)) P_k(xᵢ) wᵢ =`
`∑_j c_j (K_k δ_{ik}) = c_k K_k`, So 
`c_k = (1/K_k) ∑ᵢ f(xᵢ)P_k(xᵢ) wᵢ`
"""
function cks(::Val{:interpolating}, ::Type{P}, f, n::Int) where {P  <:  AbstractOrthogonalPolynomial}
    
    xs, ws = gauss_nodes_weights(P, n)
    return [sum(f(xⱼ) * basis(P, k)(xⱼ) * wⱼ for (xⱼ,wⱼ) in zip(xs,ws)) / norm2(P,k) for k in 0:n]
end


"""
    cks(::Val{:lsq}, ::Type{P}, f, n::Int)

Fit `f` with a polynomial `∑ᵢⁿ cᵢ Pᵢ` chosen so `<f-p,f-p>_w` is as small as possible. Using the normal equations, the coefficients are found to be
`c_k = <f,P_k>_w / <P_k,P_k>_w. For some families an approximation to the inner product, `<f,P_k>_w` may be used.


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



