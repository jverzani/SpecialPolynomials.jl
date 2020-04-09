"""
    Newton{S,T}

A [Newton](https://en.wikipedia.org/wiki/Newton_polynomial)
interpolating polynomial uses a basis `1`, `(x-x_0)`,
`(x-x_0)(x-x_1)`, ..., `(x-x0)(x-x1)⋅⋅⋅(x-x_{n-1})` and coefficients
(in forward form) `f[x_0]`, `f[x_0,x_1]`, ...,`f[x_0,...,x_n]`. The
Newton class stores the nodes (after sorting) and the Newton tableau
used to generate the coefficients on fitting.

The easiest way to construct an instance is with `fit`, as in:

```jldoctest
julia> xs = [1,2,3,4]; f(x)= x^3 - 2x + 1;

julia> p = fit(Newton, xs, f)
Newton(5.0⋅p_1(x) + 6.0⋅p_2(x) + 1.0⋅p_3(x))

julia> p.(xs) == f.(xs)  # p interpolates
true

julia> convert(Polynomial, p)
Polynomial(1.0 - 2.0*x + 1.0*x^3)
```

"""
struct Newton{N, S<:Number, T <: Number} <: AbstractInterpolatingPolynomial{T}
    xs::Vector{S}
    tableau::Matrix{T}
    var::Symbol
    function Newton(xs::Vector{S}, nt::AbstractMatrix{T}, var::Symbol=:x) where {S,T}
        N = length(xs)
        xs = sort(unique(xs))
        N != length(xs) && throw(ArgumentError("xs must be unique"))
        new{N,S,T}(xs, nt,  var)
    end
end

export Newton

basis_symbol(::Type{<:Newton}) = "p"

## Boilerplate code reproduced here, as there are two type parameters
Base.convert(::Type{P}, p::P) where {P <: Newton} =  p
Base.promote_rule(::Type{Newton{N, S, T}}, ::Type{U}) where {N, S, T, U <: Number} = Newton{N, S, promote_type(T,U)}

Polynomials.domain(::Type{<:Newton}) = Polynomials.Interval(-Inf, Inf)
Polynomials.variable(p::Newton{S, T}, var::Polynomials.SymbolLike=:x) where {S, T} = fit(Newton, p.xs, x -> x, var=var)

"""
    newton_tableau(f, xs::Vector)
    newton_tableau(xs, ys)

Construct the
[matrix](https://en.wikipedia.org/wiki/Divided_differences#Matrix_form)
form of Newton's divided differences. This assumes the xs are unique. 

"""
newton_tableau(f, xs::Vector{S}) where {S} = newton_tableau(xs, f.(xs)/one(S))
function newton_tableau(xs::Vector{S}, ys::Vector{T}) where {S,T}
    xs = sort(unique(xs))
    n = length(xs)
    M = diagm(0=>ys/one(S))
    for i in 2:n
        for j in 1:n-i+1
            k = i + j -1
            M[j, k] = (M[j+1,k] - M[j,k-1]) / (xs[j-1+1+(i-1)] - xs[j-1+1])
        end
    end
    M
end

"""
    update_newton_tableau(x, fx, xs, nt)

Update the newton tableau, `nt`, generated from the `xs` for some
function `f`. The point `(x,fx)` is a new `x` value (from the `xs`)
and `fx` the function value.

"""
function update_newton_tableau(x::T, fx::S, xs, nt::Matrix{R}) where {T, S, R}
    n = size(nt, 1)
    M = zeros(R, n+1, n+1)
    for i in 1:n
        for j in i:n
            M[i,j] = nt[i,j]
        end
    end
    M[end,end] = fx
    xxs = vcat(xs, x)
    for i  in n:-1:1
        M[i, n+1] = (M[i+1, n+1] - M[i,n]) / (xxs[end] - xxs[end-(n-i+1)])
    end

    xxs, M
end

## We use the forward difference formula
## the coefficients are the top row of the upper
## triangular tableau
Polynomials.coeffs(p::Newton) = p.tableau[1,:]

## Evaluation
## p(x) = f[x0] + f[x0,x1](x-x0) + ... + f[x0,x1,...,x_n](x-x0)(x-x1)...(x-x_{n-1})
## use nested evaluation (http://pages.cs.wisc.edu/~amos/412/lecture-notes/lecture08.pdf)
function (p::Newton{N,S,T})(x) where {N,S,T}
    xs, cs = p.xs, coeffs(p)
    tot = cs[end]
    for j in N-1:-1:1
        tot = cs[j] + (x - xs[j])*tot
    end
    tot
end

## compute (p,dp)  at x:
## http://pages.cs.wisc.edu/~amos/412/lecture-notes/lecture08.pdf
function _dual(p::Newton{N, S,T}, x) where {N,S,T}
    xs, cs = p.xs, coeffs(p)
    tot = cs[end]
    totp = zero(S)
    for j in N-1:-1:1
        delta = (x - xs[j])
        tot, totp = cs[j] + delta * tot, tot + delta * totp
    end
    tot, totp
end

    

function Base.:+(p1::Newton{N,S,T}, p2::Newton{M,S1,T1}) where {N,S,T,M,S1,T1}
    p1.var == p2.var || throw(ArgumentError("p1 and p2 must share the same variable"))
    if N == M && p1.xs == p2.xs
        tableau = p1.tableau + p2.tableau
        return Newton(p1.xs, tableau, p1.var)
    else
        p,q = N >= M ? (p1, p2) : (p2, p1)
        qq = fit(Newton, p.xs, q) # stretch than add
        ## https://en.wikipedia.org/wiki/Divided_differences
        ## T_{f+g}x = T_fx + T_gx
        p + qq
    end

end

function Base.:*(p1::Newton{N,S,T}, p2::Newton{M,S1,T1}) where {N,S,T,M,S1,T1}
    ## Basic idea
    ## find larger of xs, expand to xxs
    ## refit p1 and p2 on this
    ## use  product of tableaus
    p,q = N >= M ? (p1, p2) : (p2, p1)
    xs = p.xs
    new_xs = _new_nodes(xs, q.xs)

    xxs = vcat(xs, new_xs)
    pp = fit(Newton, xxs, p)
    qq = fit(Newton, xxs, q)

    ## https://en.wikipedia.org/wiki/Divided_differences
    ## T_{f*g}(x) = T_f(x) * T_g(x)

    Newton(pp.xs, pp.tableau * qq.tableau, p.var)
    
end

## Scalar
function Base.:+(p::Newton{N,T,S}, c::Number) where {N,T,S}
    tableau = copy(p.tableau)
    for i in 1:N
        tableau[i,i] += c
    end
    Newton(p.xs, tableau, p.var)
end

function Base.:*(p::P, c::Number) where {P <: Newton}
    Newton(p.xs, c*p.tableau, p.var)
end

Base.:-(p::P) where {P<:Newton} = p*(-1)


function  Base.convert(Q::Type{<:Polynomial},  p::Newton{N, S,T})  where  {N,S, T}
    p(variable(Polynomial{S}))
end

# return derivative
# to evaluate derivative at a value
# use  _dual
# These are not efficient due to the need to  compute a tableau
function Polynomials.derivative(p::Newton{N,S,T}, order::Int=1) where {N,S,T}
    order < 0 && throw(ArgumentError("order must be positive"))
    order == 0 && return p
    N==1 &&  return fit(Newton, p.xs, zero)
    
    x = variable(Polynomial{S})
    q,qq = _dual(p, x)
    dp = fit(Newton, p.xs[1:end-1], qq)
    if order > 1
        dp = derivative(dp, order-1)
    end
    dp
    
end

# integrate
function Polynomials.integrate(p::Newton{N,S,T}, C::Number=0) where {N,S,T}
    q = convert(Polynomial, p)
    Q = integrate(q, C)
    if length(p.xs) >= 2
        new_x = p.xs[end] + (p.xs[end] - p.xs[end-1])
        xxs = vcat(p.xs, new_x)
    else
        xxs = vcat(p.xs, p.xs[1]+1)
    end
    fit(Newton, xxs, Q, var=p.var)
    
end


function Polynomials.fit(P::Type{<:Newton},
                         xs::AbstractVector{S},
                         ys::AbstractVector{T};
                         var = :x,) where {S, T}
    ind = sortperm(xs)
    nt = newton_tableau(xs[ind], ys[ind])

    Newton(xs, nt, var)
end


function Polynomials.fit(P::Type{<:Newton}, xs::AbstractVector{S}, f; var=:x) where {S}
    xs = sort(xs)
    nt = newton_tableau(f, xs)
    Newton(xs, nt, var)
end
