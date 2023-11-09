"""
    glaser_liu_rokhlin_gauss_nodes(π::P, x0=pqr_start(P,degree(π));
                                         symmetry=pqr_symmetry(P),
                                         m=5)

Find Gauss nodes and weights of the  basis  polynomial  `π = basis(P,n)`. The nodes  are the  roots of the  polynomial and the weights are computed  from known formulas.

The method from  [Glaser,  Liu, and Rokhlin](DOI: 10.1137/06067016X) is  used. This is  an  O(n)  method (whereas the method  basedon the Jacobi  matrix is O(n^2)). This method fits easily  into the framework  provided  through  the `AbstractCCOP` types.  The [FastGaussQuadrature](https://github.com/JuliaApproximation/FastGaussQuadrature.jl) package provides even more efficient  algorithms and  pays attention to numerical  issues, examples of  which  can  be found in  [Hale and Townsend](https://core.ac.uk/download/pdf/9403469.pdf). The `FastGuassQuadrature` package is used when it is loaded.


An example:

XXX not the case anymore, pass in a basis
```
function  gauss_nodes_weights(P::Type{<:Jacobi{α,β}}, n)  where {α,β}
    # we don't  have a  good starting point  unless α=β
    if α == β
        xs, ws = glaser_liu_rokhlin_gauss_nodes(basis(MonicJacobi{α,β},n))
        λ = kn(P,n)^2
        xs, ws/λ
    else
        J = jacobi_matrix(P, n)
        eig = eigen(J, extrema(P)...)
        wts =  π̃βn(P,0) * (eig.vectors[1,:]).^2
        eig.values,  wts
    end

end


pqr_start(P::Type{MonicJacobi{α,α}}, n) where {α} =  0.0 # zero(eltype(P))
pqr_symmetry(P::Type{<:MonicJacobi{α,α}}) where {α} = true
function pqr_weight(P::Type{<:MonicJacobi{α,α}}, n, x, dπx) where {α}
    # https://www.chebfun.org/publications/HaleTownsend2013a.pdf
    # Eqn  (1.4)
    # Cnαβ should be  computed  using asymptotic formula  for larger n (§3.2.3)
    # XXX TThere i ss ome constant  that  makes this not work....
    β = α
    Cnαβ = 2^(α+β+1)
    Cnαβ *= gamma(n + α +  1) * gamma(n + β + 1) / gamma(n + α + β + 1)
    Cnαβ /= gamma(1 + n)
    val = Cnαβ / (1-x^2) / dπx^2
    val
end
```

"""
glaser_liu_rokhlin_gauss_nodes()

"""
    pqr

Compute  `p,q,r,p',q',r'` where `p(x)y'' + q(x)y' + r(x)  = 0`. Uses
the `a,b,c,d,e` values characterizing the  system.

"""
function pqr(B::Type{<:AbstractCCOPBasis}, n, x)
    a, b, c, d, e = abcde(B)

    p  = a * x^2 + b * x + c
    dp = 2a * x + b
    q  = d * x + e
    dq = d
    r  = -(a * n * (n - 1) + d * n)
    dr = 0
    (p, q, r, dp, dq, dr)
end

"""
    pqr_start(p::P, n)

A starting  value for finding the gauss nodes
"""
pqr_start(::Type{B}, n) where {B<:AbstractCCOPBasis} = 0

"""
    pqr_symmetry(p::P)

Boolean to specify if symmetry should apply to output
"""
pqr_symmetry(::Type{B}) where {B<:AbstractCCOPBasis} = false

"""
    pqr_weight(p::P, x, dx)

Compute weight from `x`, `dπx` -- `x` a node and `dπx` the derivative's value at the node
"""
pqr_weight(p::P, n, x, dx) where {P} = throw(ArgumentError("Not implemented"))

# We just use clenshaw to  evaluate  `p(x), dp(x)` for Newton's method below,
# but implementing the following can be more efficient
# ## Evaluate basis(P,n)  and its  derivative using the three-point recursion representation
# function orthogonal_polyval_derivative(p::P, x::S) where {P <: AbstractOrthogonalPolynomial, S}

#     T = eltype(p)
#     oS = one(x)
#     R = eltype(one(T) * oS)

#     length(p) == 0 && return (zero(R), zero(R))
#     d = length(p) - 1

#     # P0,  P_{-1} case
#     pn,  pn_1 = P0(P,x), zero(x)
#     dpn, dpn_1  = dP0(P,x), zero(x)

#     ptot = p[0]*P0(P,x)
#     dptot  = p[0]*dP0(P,x)

#     ## we compute  forward here, finding pi, dp_i,  p_{i+1}, dp_{i+1}, ...
#     for i in 0:d-1
#         an, bn, cn=  An(p,i), Bn(p,i), Cn(p,i)
#         dpn_1, dpn = dpn,  an*pn + (an*x + bn)*dpn + cn*dpn_1
#         pn_1,  pn =  pn, (an*x + bn)*pn + cn*pn_1
#         ptot += p[i+1]*pn
#         dptot += p[i+1]*dpn
#     end
#     return (ptot, dptot)
# end

# run Newton's method to find zero of p
function newton(p::P, x0::S) where {B<:AbstractCCOPBasis,P<:AbstractUnivariatePolynomial{B},S}
    maxsteps = 25
    dp = derivative(p)
    while maxsteps > 0
        px, dpx = p(x0), dp(x0) #orthogonal_polyval_derivative(p, x0)
        Δ = px / dpx
        isnan(Δ) && return (x0, dpx)
        if isinf(Δ)
            x0 += sqrt(eps())  # give a nudge from minimum or maximum
            continue
        end
        x0 -= Δ
        min(abs(Δ), abs(px)) <= sqrt(eps(S)) / 100 && return x0, dp(x0)# orthogonal_polyval_derivative(p, x0)[2]
        maxsteps -= 1
    end
    @warn "Newton's method did not converge from $x0"
    NaN
end

# Second-order Runge-Kutta with global truncation error of h^2
function RK(t0, x0, F, h, n)
    k0 = h * F(t0, x0)
    local x1
    for i in 1:n
        t1 = t0 + h
        k1 = h * F(t0 + h, x0 + k0)
        x1 = x0 + 1 / 2 * (k0 + k1)
        t0, x0, k0 = t1, x1, k1
    end
    x1
end

function prufer(Π::Type{B}, n) where {B<:AbstractCCOPBasis}
    (θ, x) -> begin
        dom = domain(Π)
        a, b = first(dom) + eps(), last(dom) - eps()
        x = clamp(x, a, b)
        p, q, r, dp, dq, dr = pqr(Π, n, x)
        -inv(sqrt(r / p) + (dr * p - dp * r + 2r * q) / (2r * p) * sin(2θ) / 2)
    end
end

## Compute gauss node using  algorithm from
## A Fast Algorithm for the Calculation of the Roots of Special Functions
## by Glaser, Liu, Rokhlin
## DOI: 10.1137/06067016X
## specialized to the orthogonal polynomial case
function glaser_liu_rokhlin_gauss_nodes(
    π::P,
    x0=pqr_start(B, degree(π));
    symmetry=pqr_symmetry(B),
    m=5,
) where {B<:AbstractCCOPBasis,P<:AbstractUnivariatePolynomial{B}}
    n = degree(π)
    F = prufer(B, n)
    x = float(x0)

    # step 1 to get initial might be newton or might be RK
    if symmetry && iseven(n)
        # a max at π(x0)
        x = RK(0, x0, F, -(pi / 2) / m, m)
    end
    x1, dπ1 = newton(π, x)
    rts, dπrts = [x1], [dπ1]
    N = symmetry ? div(n + 1, 2) : n
    for i in 2:N
        xx = RK(pi / 2, rts[end], F, -pi / m, m)
        x, dπx = newton(π, xx)
        push!(rts, x)
        push!(dπrts, dπx)
    end
    if symmetry
        if iseven(n)
            rts = vcat(-reverse(rts), rts)
            dπrts = vcat(-reverse(dπrts), dπrts)
        else
            rts = vcat(-reverse(rts[2:end]), rts)
            dπrts = vcat(-reverse(dπrts[2:end]), dπrts)
        end
    end
    # get weights
    weights = pqr_weight.(B, n, rts, dπrts)
    rts, weights
end
