# generic classical discrete orthogonal polynomial, CDOP

abstract type AbstractCDOP{T,X,N} <: AbstractCOP{T,X,N} end
"""
     AbstractCDOP{T,X,N} <: AbstractCOP{T,X,N}

Following [Koepf  and Schmersau](https://arxiv.org/pdf/math/9703217.pdf), a family `y(x)=p_n(x)=k_x⋅x^n +  ...`  
for  `n  ∈  {0, 1,…}, k_n ≠ 0` of polynomials is a family of classic *discrete* orthogonal polynomials if it  is  a
solution of a differential equation

`(a⋅x²+b⋅x+c) ⋅ Δ∇y + (d⋅x + e) ⋅ ∇' + λᵢ⋅ y = 0`,

where  `Δy(x) = y(x+1) - y(x)` and `∇y(x) = y(x) - y(x-1)`.

A family is characterized by the 5 coefficients: `a,b,c,d,e`.
Let `σ = (a⋅x²+b⋅x+c)`, `τ = (d⋅x + e).`

As in the classical-continuous-orthogonal-polynomial case
[`AbstractCCOP`](@ref), from these 5 values the cofficients in the
there-point recursion, and other structural equations can be
represented. These allow polynomial multiplication, integration,
differentiation, conversion, etc. to be defined generically.

[Koekoek and Swarttouw](https://arxiv.org/pdf/math/9602214.pdf)
present an encyclopedia of formula characterizing families of
orthogonal polynomials. 

For example, on p29 they give  formulas for Hahn polynomials through:

`n(n+α+β+1)y(x) = B(x)y(x+1) -[B(x)+D(x)]y(x) + D(x)y(x-1)`,  with  explicit values  for  `B` and `D`. Reexpressing gives:
`BΔy(x) - D∇y(x) -λ y(x)  = 0`. From the rexpressed Eqn (4) for Koepf & Schemersau we have the identification:
`σ+τ =  B; σ=D`,  so  `τ=B-D`. From this `a,b,c,d,e` can be  gleaned.

The above, is termed the eigevalue equation (e.g. [Goertz and Offner](https://arxiv.org/pdf/1609.07291.pdf)), as it can be reexpressed as 

`Δ(D(x)⋅ω(x)⋅∇yᵢ(x) = λᵢ⋅ω(x)⋅yᵢ(x)`




"""
AbstractCDOP

## Implemented  families
##                    a     b     c      d      e
##                    -------------------------------
## Charlier{μ}        0     1     0     -1      μ
## Meixner{γ,μ}       0     1     0    μ-1     γ⋅μ
## Krawchouk{p,𝐍}     0    p-1    0     1      p/𝐍
## Hahn{α,β,𝐍}        1 -(β+𝐍+1)  0   α+β+2  -𝐍⋅(α+1)

# subtypes  to keep track of number of parameters
# passed to  @registerN macros
abstract type AbstractCDOP0{T,X,N} <: AbstractCDOP{T,X,N} end
abstract type AbstractCDOP1{α,T,X,N} <: AbstractCDOP{T,X,N} end
abstract type AbstractCDOP2{α,β,T,X,N} <: AbstractCDOP{T,X,N} end
abstract type AbstractCDOP3{α,β,γ,T,X,N} <: AbstractCDOP{T,X,N} end

⟒(P::Type{<:AbstractCDOP1{α}}) where {α} = constructorof(P){α}
⟒(P::Type{<:AbstractCDOP2{α,β}}) where {α,β} = constructorof(P){α,β}
⟒(P::Type{<:AbstractCDOP3{α,β,γ}}) where {α,β,γ} = constructorof(P){α,β,γ}

# for conversion to base case
# \upstigma[tab]
ϛ(P::Type{<:AbstractCDOP}) = FallingFactorial

# compose with FallingFactorial
function Base.convert(::Type{Q}, p::P) where {Q<:AbstractCDOP,P<:AbstractCDOP}
    T = eltype(P)
    _convert_cop(Q, _convert_cop(FallingFactorial{T}, p))
end

Δₓ(p::AbstractCDOP) = p(variable(p) + 1) - p(variable(p))
∇ₓ(p::AbstractCDOP) = p(variable(p)) - p(variable(p) - 1)

function innerproduct(P::Type{<:AbstractCDOP}, f, g)
    a, b = extrema(P)
    (isinf(a) || isinf(b)) &&
        throw(ArgumentError("Only defined for discrete polynomials  with a finite domain"))
    ω = weight_function(P)
    sum(f(i) * g(i) * ω(i) for i in a:b)
end

## Structural Equations
## https://dlmf.nist.gov/18.22
## we have σ⋅Δ∇p + τ⋅∇p + λp =  0
## Using parameterization   of  above, τ = C, σ = A - C
function Ãn(P::Type{<:AbstractCDOP}, a, b, c, d, e, n::Int)
    one(eltype(P))
end

function B̃n(P::Type{<:AbstractCDOP}, a, b, c, d, e, n::Int)
    S = eltype(P)

    num = n * (d + 2b) * (d + a * n - a) + e * (d - 2a)
    den = (2a * n - 2a + d) * (d + 2a * n)

    iszero(den) && return B̃n(P, Val(n))

    val = one(S) * num / den

    val
end

function C̃n(P::Type{<:AbstractCDOP}, a, b, c, d, e, n::Int)
    S = eltype(P)

    num = one(S)
    num *= (n - 1) * (d + a * n - a)
    num *= a * n * d - d * b - a * d + a^2 * n^2 - 2a^2 * n + 4c * a + a^2 + 2e * a - b^2
    num += -d * b * e + d^2 * c + a * e^2

    num *= (a * n + d - 2a) * n
    den = (d - a + 2a * n) * (d + 2a * n - 3a) * (2a * n - 2a + d)^2

    iszero(den) && return C̃n(P, Val(n))

    val = -one(S) * num / den

    val
end

# αn, βn,γn
# σ⋅pn' = [αn, βn,γn] ⋅ [p_{n+1},p_n,p_{n-1}]
#αn(P::Type{<:AbstractCDOP}, n::Int) = α̃(P,n) / k1k0(P,n) 
function α̃n(P::Type{<:AbstractCDOP}, n::Int)
    a, b, c, d, e = abcde(P)
    S = eltype(P)

    num = a * n
    val = one(S) * num

    return val
end

function β̃n(P::Type{<:AbstractCDOP}, n::Int)
    a, b, c, d, e = abcde(P)
    S = eltype(P)

    num =
        -n *
        (d + a * n - a) *
        (2 * a * n * d - a * d - d * b + 2e * a - 2a^2 * n + 2a^2 * n^2)
    den = (2a * n - 2a + d) * (d + 2a * n)

    iszero(den) && return βn(P, Val(n))

    val = (one(S) * num) / den

    return val
end

#γn(P::Type{<:AbstractCDOP}, n::Int) =  γ̃n(P,n) * k1k0(P,n-1)
function γ̃n(P::Type{<:AbstractCDOP}, n::Int)
    a, b, c, d, e = abcde(P)
    S = eltype(P)

    num = (n - 1) * (d + a * n - a)
    num *= (a * n * d - d * b - a * d + a^2 * n^2 - 2a^2 * n + 4c * a + a^2 + 2e * a - b^2)
    num += -d * b * e + d^2 * c + a * e^2
    num *= (d + a * n - a) * (a * n + d - 2a) * n
    den = (d - a + 2a * n) * (d + 2a * n - 3a) * (2a * n - 2a + d)^2
    iszero(den) && return γn(P, Val(n))

    val = one(S) * num / den

    return val
end

function abcdeᴵ(P::Type{<:AbstractCDOP})
    a, b, c, d, e = abcde(P).a, abcde(P).b, abcde(P).c, abcde(P).d, abcde(P).e
    NamedTuple{(:a, :b, :c, :d, :e)}((a, b, c, d + 2a, d + e + a + b))
end

function ẫn(P::Type{<:AbstractCDOP}, n::Int)
    a, b, c, d, e = abcde(P)
    S = eltype(P)

    val = one(S)
    val /= n + 1
    return val
end

function b̂̃n(P::Type{<:AbstractCDOP}, n::Int)
    a, b, c, d, e = abcde(P)
    S = eltype(P)

    num = -2 * a * n * (d + a * n - a) - d * b + a * d - d^2 + 2e * a
    den = (2a * n - 2a + d) * (d + 2a * n)
    iszero(den) && return b̂̃n(P, Val(n))

    val = one(S) * num / den
    return val
end

function ĉ̃n(P::Type{<:AbstractCDOP}, n::Int)
    a, b, c, d, e = abcde(P)
    S = eltype(P)

    num =
        (n - 1) *
        (d + a * n - a) *
        (a * n * d - d * b - a * d + a^2 * n^2 - 2a^2 * n + 4c * a + a^2 + 2e * a - b^2)
    num += -d * b * e + d^2 * c + a * e^2
    num *= a * n
    den = (d - a + 2a * n) * (d + 2a * n - 3a) * (2a * n - 2a + d)^2
    iszero(den) && return ĉ̃n(P, Val(n))

    val = (one(S) * num) / den
    return val
end

function ⊗(::Type{<:AbstractCDOP}, p::P, q::Q) where {P<:AbstractCDOP,Q<:AbstractCDOP}

    #@assert  ⟒(P) == ⟒(Q)
    #@assert eltype(p) == eltype(q)
    isconstant(p) && return q * constantterm(p)
    isconstant(q) && return p * constantterm(q)
    assert_same_variable(p, q)

    convert(⟒(P), convert(FallingFactorial, p) * convert(FallingFactorial, q))

    #    R = eltype(p)
    #    convert(⟒(P){R}, convert(Polynomial, p) * convert(Polynomial, q))
end
