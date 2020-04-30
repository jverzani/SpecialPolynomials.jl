##################################################
#  Implement  a method from  A. Glaser, X. Liu, and V. Rokhlin. "A fast algorithm for the calculation of the roots of special functions." SIAM J. Sci. Comput., 29 (2007), 1420-1438.
#
# See https://github.com/JuliaApproximation/FastGaussQuadrature.jl for a more thorough, definitive treatment of finding gauss nodes  and weights.

"""
    pqr(x)

For a family of orthogonal  polynomials satisfying the differential equation  `p(x)P''_n + q(x)P'_n + r(x)P_n=0`,
return the function `(x,n) -> (p,q,r,dp, dq, dr).

"""
pqr(p::P) where {P} = throw(ArgumentError("pqr not defined for polynomials of type $P"))

"""
   pqr_scale(p::P)

To scale a polynomial family so that π_n(x) = s_n(x) P_n(x) we need a function  to compute
`s_{n+1}(x)/s_n(x)`, `s_{n+1}(x)/s_{n-1}(x)`, and the  two derivatives in `x`. This returns them.
"""
pqr_scale(p::P) where {P} = (x,n) -> (one(x), one(x), zero(x), zero(x))

"""
    pqr_start(p::P, n)

A starting  value for finding the gauss nodes
"""
pqr_start(p::P, n) where {P} = 0

"""
    pqr_symmetry(p::P)

Boolean to specify if symmetry should apply to output
"""
pqr_symmetry(p::P) where {P} = false

"""
    pqr_weight(p::P, x, dx)

Compute weight from x, dπx, x a node and dπx the derivative's value at the node
"""
pqr_weight(p::P, n, x, dx) where {P} = throw(MethodError())



# run Newton's method to find zero of p
function newton(p::P, x0::S) where  {P <: AbstractOrthogonalPolynomial, S}
    maxsteps = 25
    while maxsteps > 0
        px, dpx = orthogonal_polyval_derivative(p, x0)
        Δ = px/dpx
        isnan(Δ) && return (x0, dpx)
        if isinf(Δ)
            x0 += sqrt(eps())  # give a nudge from minimum or maximum
            continue
        end
        x0 -= Δ
        abs(px) <= 1e3*eps(S) && return  x0, orthogonal_polyval_derivative(p, x0)[2]
        maxsteps -= 1
    end
    @warn "Newton's method did not converge from $x0"
    NaN
end

# Second-order Runge-Kutta with global truncation error of h^2
function RK(t0,  x0, F, h, n)
    k0 = h*F(t0,x0)
    local x1
    for i in 1:n
        t1 = t0 + h
        k1 = h * F(t0+h,x0+k0)
        x1 = x0 + 1/2 * (k0+k1)
        t0, x0, k0 = t1, x1, k1
    end
    x1
end

## Compute gauss node using  algorithm from  
## A Fast Algorithm for the Calculation of the Roots of Special Functions
## by Glaser, Liu, Rokhlin
## DOI: 10.1137/06067016X
## specialized to the orthogonal polynomial case 
function glaser_liu_rokhlin_gauss_nodes(π, x0=pqr_start(π,degree(π)); symmetry=pqr_symmetry(π), m=5)
    n = degree(π)
    dom = domain(π)
    a, b = first(dom)+eps(), last(dom)-eps()
    F = (θ, x) -> begin
        x = clamp(x, a, b)
        p,q,r,dp,dq,dr = pqr(π)(x, n)
        -inv(sqrt(r/p) + (dr*p - dp*r + 2r*q)/(2r*p) * sin(2θ)/2)
    end
    x = float(x0)

    # step 1 to get initial might be newton or might be RK
    if symmetry && iseven(n)
        # a max at π(x0)
        x = RK(0, x0, F, -(pi/2)/m, m)
    end
    x1, dπ1 = newton(π, x)
    rts, dπrts = [x1], [dπ1]
    N = symmetry ? div(n+1,2) : n
    for i in 2:N
        xx = RK(pi/2, rts[end], F, -pi/m, m)
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
    weights = pqr_weight.(π, n, rts, dπrts)
    rts,  weights
end
