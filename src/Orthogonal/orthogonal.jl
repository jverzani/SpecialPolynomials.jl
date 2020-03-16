abstract type OrthogonalPolynomial{T} <: AbstractSpecialPolynomial{T} end



"""
    weight_function(p)
    weight_function(::Type{P})

For orthogonal polynomials, a function w with `integrate Bn Bm w dt = 0` when n and m are not equal
"""
weight_function(::P) where {P <: OrthogonalPolynomial} = weight_function(P)
weight_function(::Type{P}) where {P <: OrthogonalPolynomial} = throw(MethodError("Not implemented"))

"""
    generating_function(p)
    generating_function(::Type{P})

(t,x) -> sum(t^n Pn(x) for n in 0:oo)
"""
generating_function(::P) where {P <: OrthogonalPolynomial} = generating_function(P)
generating_function(::Type{P}) where {P <: OrthogonalPolynomial} = throw(MethodError("Not implemented"))

# cf. https://en.wikipedia.org/wiki/Orthogonal_polynomials#Recurrence_relation
# We have this () for orthogonal polynomials
# Pn = (An*x + Bn) * P_{n-1} + Cn P_{n-2}
# which is used to express in terms of c0 Po + c1 P1
An(::Type{P}, n) where {P <: OrthogonalPolynomial} = throw(MethodError("Not implemented"))
Bn(::Type{P}, n) where {P <: OrthogonalPolynomial} = throw(MethodError("Not implemented"))
Cn(::Type{P}, n) where {P <: OrthogonalPolynomial} = throw(MethodError("Not implemented"))
P0(::Type{P}, x) where {P <: OrthogonalPolynomial} = throw(MethodError("Not implemented"))
P1(::Type{P}, x) where {P <: OrthogonalPolynomial} = throw(MethodError("Not implemented"))

## Evaluate an orthogonal polynomial
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
    length(ch) == 2 && return ch[0] + ch[1]*x

    c0 = ch[end - 1]
    c1 = ch[end]
    @inbounds for i in degree(ch):-1:2
        c0, c1 = ch[i - 2] + c1 * Cn(P, i), c0 + c1 * (An(P, i) * x + Bn(P, i))
    end
    # Now have c0 P_0(x) + c1 P_1 = c0 + c1*x
    return c0 * P0(P, x) + c1 * P1(P, x)
end

function Base.convert(P::Type{<:Polynomial}, ch::OrthogonalPolynomial)
    if length(ch) == 1
        return P(ch.coeffs, ch.var)
    end
    T = eltype(ch)
    Q = typeof(ch)
    ##
    ## P_n = (An x + Bn)P_{n-1} + CnP_{n-2}
    ## so ... c0 P_{n-1} + c1 P_n
    ## c1 = c0 + (An x+ Bn)
    ## c0 = p[i-1] + Cn
    c0 = P(ch[end - 1], ch.var)
    c1 = P(ch[end], ch.var)
    x = variable(P)
    @inbounds for i in degree(ch):-1:2
        c0, c1 = ch[i - 2] + c1 * Cn(Q, i), c0 + c1 * (An(Q, i) * x + Bn(Q, i))
    end
    return c0 *  P0(Q,  x) + c1 * P1(Q, x)
end



function Polynomials.vander(p::Type{P}, x::AbstractVector{T}, n::Integer) where {P <:  OrthogonalPolynomial, T <: Number}
    A = Matrix{T}(undef, length(x), n + 1)
    A[:, 1] .= P0(P, one(T))
    if n > 0
        A[:, 2] .= P1(P, x)
        @inbounds for i in 3:n + 1
            A[:, i] .= A[:, i - 1] .* (An(P, i)*x .+ Bn(P,i)) .+ (Cn(P,i) * A[:, i - 2])
        end
    end
    return A
end

### This is too slow...
function innerproduct(P::Type{<: OrthogonalPolynomial}, f, g)
    dom = domain(P)
    fn = x -> f(x) * g(x) * weight_function(P)(x)
    a,err =  quadgk(fn, first(dom), last(dom))
    a
end

dot(f, p::P) where {P <: OrthogonalPolynomial} = innerproduct(P, f, p)
function dot(p::P, q::P) where {P <: OrthogonalPolynomial}
    # use orthogonality  here

    n, m = degree(p), degree(q)
    tot = zero(eltype(P))
    for i in 0:minimum(n,m)
        tot += p[i] * q[i] * innerproduct(P, Polynomials.basis(P, i), Polynomials.basis(P, i))
    end

    tot
end

Base.abs2(p::P) where {P <: OrthogonalPolynomial} = innerproduct(P, p, p)

normdot(f, p::P) where {P <: OrthogonalPolynomial} = dot(f, p)/abs2(p)

function Polynomials.fit(P::Type{<:OrthogonalPolynomial}, f, deg::Int; var=:x)
    # use least squares fit of orthogonal polynomial family to f

    cs = [normdot(f, Polynomials.basis(P,i)) for i  in 0:deg]
    P(cs, var)
end
