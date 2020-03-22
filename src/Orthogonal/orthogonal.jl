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
An(p::P, n) where {P <: OrthogonalPolynomial} = An(P,n)
Bn(p::P, n) where {P <: OrthogonalPolynomial} = Bn(P,n)
Cn(p::P, n) where {P <: OrthogonalPolynomial} = Cn(P,n)

P0(::Type{P}, x) where {P <: OrthogonalPolynomial} = one(x)
P0(p::P, x) where {P <: OrthogonalPolynomial} = one(x)
P1(::Type{P}, x) where {P <: OrthogonalPolynomial} = An(P,0)*x + Bn(P,0)
P1(p::P, x) where {P <: OrthogonalPolynomial} = An(p,0)*x + Bn(p,0)

function alphas_betas(::Type{P}, n) where {P <: OrthogonalPolynomial}
    ani = An(P,1)
    alpha0 = -Bn(P,1)/ani
    T = eltype(alpha0)
    alphas = zeros(T, n)
    betas = zeros(T, n-1)
    alphas[1] = alpha0
    for i in 1:n-1
        ani_1 = ani
        ani = An(P,i)
        alphas[i+1] = -Bn(P,i)/ani
        betas[i] = sqrt(-Cn(P,i)/ani/ani_1)
    end

    (alphas, betas)
end


# Jacobi matrix
# roots(basis(P,n)) = eigvals(jacobi_matrix(P,n)), but this is more stable
function jacobi_matrix(::Type{P}, n) where {P <: OrthogonalPolynomial}
    LinearAlgebra.SymTridiagonal(alphas_betas(P, n)...)
end

## Compute <p_i, p_i> = \| p \|^2; allows export, but work is in norm2
Base.abs2(::Type{P}, n) where {P <: OrthogonalPolynomial} = norm2(P, n)

## Compute <p_i, p_i> = \| p \|^2
## Slow default; generally should  be directly expressed for each family
function norm2(::Type{P}, n) where {P <: OrthogonalPolynomial}
    p = Polynomials.basis(P,n)
    innerproduct(P, p, p)
end

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
        c0, c1 = ch[i - 2] + c1 * Cn(ch, i), c0 + c1 * (An(ch, i) * x + Bn(ch, i))
    end
    # Now have c0 P_0(x) + c1 P_1 = c0 + c1*x
    return c0 * P0(ch, x) + c1 * P1(ch, x)
end

function Base.convert(P::Type{<:Polynomial}, ch::OrthogonalPolynomial)
    if length(ch) == 1
        return P(ch.coeffs, ch.var)
    end
    T = eltype(ch)
    ##
    ## P_n = (An x + Bn)P_{n-1} + CnP_{n-2}
    ## so ... c0 P_{n-1} + c1 P_n
    ## c1 = c0 + (An x+ Bn)
    ## c0 = p[i-1] + Cn
    c0 = P(ch[end - 1], ch.var)
    c1 = P(ch[end], ch.var)
    x = variable(P)
    @inbounds for i in degree(ch):-1:2
        c0, c1 = ch[i - 2] + c1 * Cn(ch, i), c0 + c1 * (An(ch, i) * x + Bn(ch, i))
    end
    return c0 *  P0(ch,  x) + c1 * P1(ch, x)
end

## This is going  to be slow!
## Use orthogonality:
## x^n = sum( λ(P,n,i) P_i for i in 0:n)
## so <x^n, P_j>/<P_j,P_j> =  λ(P,n,j), the other terms being 0
function Base.convert(P::Type{J}, p::Polynomial) where {J <: OrthogonalPolynomial}
    d = degree(p)
    R = eltype(one(eltype(p))/1)
    qs = zeros(R, d+1)
    for i in 0:d
        qs[i+1] = sum(p[j] * λ(J, j, i) for j in i:d)
    end
    J(qs, p.var)
end

# find <x^n, P_n>
# * must define abs2(J,j) =  <P_j, P_j>
function λ(J::Type{<: OrthogonalPolynomial}, n, j)
   p = convert(Polynomial, Polynomials.basis(J, j))
   innerproduct(J, x->x^n,  p) / norm2(J, j)
end


function Polynomials.vander(p::Type{P}, x::AbstractVector{T}, n::Integer) where {P <:  OrthogonalPolynomial, T <: Number}
    A = Matrix{T}(undef, length(x), n + 1)
    A[:, 1] .= P0(P, one(T))
    if n > 0
        A[:, 2] .= P1(P, x)
        @inbounds for i in 3:n + 1
            A[:, i] .= A[:, i - 1] .* (An(p, i)*x .+ Bn(p,i)) .+ (Cn(p,i) * A[:, i - 2])
        end
    end
    return A
end

### This is too slow...
function innerproduct(P::Type{<: OrthogonalPolynomial}, f, g)
    dom = domain(P)
    fn = x -> f(x) * g(x) * weight_function(P)(x)
    a,err =  quadgk(fn, first(dom)+sqrt(eps()), last(dom)-sqrt(eps()))
    a
end

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

## return ck =  <f,P_k>/<P_k,P_k>
ck(::Type{P}, f, k::Int, n::Int) where {P  <:  OrthogonalPolynomial} = dot(P,f,k)/norm2(P,k)

function cks(::Type{P}, f, n::Int) where {P  <:  OrthogonalPolynomial}
      P_n1 = Polynomials.basis(P, n+1)
    xs::Vector{Float64} = roots(P_n1)
      fxs  = f.(xs)

      p0 = P0.(P, xs)
      p1 = P1.(P, xs)

      cks = similar(xs)
      cks[1] = dot(fxs, p0) / norm2(P,0) / pi
      cks[2] = dot(fxs, p1) / norm2(P,1) / pi


      for i in 2:n
          p0[:], p1[:] = p1, p1 .* (An(P, i)*xs .+ Bn(P,i)) .+ (Cn(P,i) * p0)
          cks[i+1] =  dot(fxs, p1) / norm2(P,i) / pi
      end
      cks / (n+1)
end



## for an orthogonal family, fit f with a0 p_0 + ... + a_d P_d
## where a_i = <f, p_i>/<p_i,p_i>
function Polynomials.fit(P::Type{<:OrthogonalPolynomial}, f, deg::Int; var=:x)
    # use least squares fit of orthogonal polynomial family to f
    cs = [ck(P,f,k,deg) for k in 0:deg]
    P(cs, var)
end
