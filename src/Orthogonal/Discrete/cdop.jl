# generic classical discrete orthogonal polynomial, CDOP

abstract type AbstractCDOPBasis <: AbstractCOPBasis end

#abstract type AbstractCDOP{T,X} <: AbstractCOP{T,X} end
"""
     AbstractCDOPBasis{T,X} <: AbstractCOPBasis{T,X}

Following [Koepf  and Schmersau](https://arxiv.org/pdf/math/9703217.pdf), a family `y(x)=p_n(x)=k_x⋅x^n +  ...`
for  `n  ∈  {0, 1,…}, k_n ≠ 0` of polynomials is a family of classical *discrete* orthogonal polynomials if it  is  a
solution of a difference equation

`(a⋅x²+b⋅x+c) ⋅ Δ∇y + (d⋅x + e) ⋅ ∇' + λᵢ⋅ y = 0`,

where  `Δy(x) = y(x+1) - y(x)` and `∇y(x) = y(x) - y(x-1)`.

A family is characterized by the 5 coefficients: `a,b,c,d,e`.
Let `σ = (a⋅x²+b⋅x+c)`, `τ = (d⋅x + e).`

As in the classical-continuous-orthogonal-polynomial case
[`AbstractCCOPBasis`](@ref), from these 5 values the coefficients in the
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

The above, is termed the eigenvalue equation (e.g. [Goertz and Offner](https://arxiv.org/pdf/1609.07291.pdf)), as it can be reexpressed as

`Δ(D(x)⋅ω(x)⋅∇yᵢ(x) = λᵢ⋅ω(x)⋅yᵢ(x)`




"""
AbstractCDOPBasis

# type for dispatch
const AbstractCDOPPolynomial =
    AbstractUnivariatePolynomial{<:AbstractCDOPBasis,T,X} where {T,X}

#AbstractCDOP

## Implemented  families
##                    a     b     c      d      e
##                    -------------------------------
## Charlier{μ}        0     1     0     -1      μ
## Meixner{γ,μ}       0     1     0    μ-1     γ⋅μ
## Krawchouk{p,𝐍}     0    p-1    0     1      p/𝐍
## Hahn{α,β,𝐍}        1 -(β+𝐍+1)  0   α+β+2  -𝐍⋅(α+1)

# subtypes  to keep track of number of parameters
# # passed to  @registerN macros
# abstract type AbstractCDOP0{T,X} <: AbstractCDOP{T,X} end
# abstract type AbstractCDOP1{α,T,X} <: AbstractCDOP{T,X} end
# abstract type AbstractCDOP2{α,β,T,X} <: AbstractCDOP{T,X} end
# abstract type AbstractCDOP3{α,β,γ,T,X} <: AbstractCDOP{T,X} end

# ⟒(P::Type{<:AbstractCDOP1{α}}) where {α} = constructorof(P){α}
# ⟒(P::Type{<:AbstractCDOP2{α,β}}) where {α,β} = constructorof(P){α,β}
# ⟒(P::Type{<:AbstractCDOP3{α,β,γ}}) where {α,β,γ} = constructorof(P){α,β,γ}

# for conversion to base case
# \upstigma[tab]
ϛ(P::Type{<:AbstractCDOPPolynomial}) = FallingFactorial

# compose with FallingFactorial
function Base.convert(
    ::Type{Q},
    p::P,
) where {Q<:AbstractCDOPPolynomial,P<:AbstractCDOPPolynomial}
    T = eltype(P)
    _convert_cop(Q, _convert_cop(FallingFactorial{T}, p))
end

Δₓ(p::AbstractCDOPPolynomial) = p(variable(p) + 1) - p(variable(p))
∇ₓ(p::AbstractCDOPPolynomial) = p(variable(p)) - p(variable(p) - 1)

# XXX Make extrema(B) only
function innerproduct(B::Type{<:AbstractCDOPBasis}, f, g)
    a, b = extrema(Polynomials.domain(B))
    (isinf(a) || isinf(b)) &&
        throw(ArgumentError("Only defined for discrete polynomials  with a finite domain"))
    ω = weight_function(B)
    sum(f(i) * g(i) * ω(i) for i in a:b)
end

## Structural Equations
## https://dlmf.nist.gov/18.22
## we have σ⋅Δ∇p + τ⋅∇p + λp =  0
## Using parameterization   of  above, τ = C, σ = A - C
function Ãn(B::Type{<:AbstractCDOPBasis}, a, b, c, d, e, n::Int)
    1
end

function B̃n(B::Type{<:AbstractCDOPBasis}, a, b, c, d, e, n::Int)
    S = Int #eltype(P)

    num = n * (d + 2b) * (d + a * n - a) + e * (d - 2a)
    den = (2a * n - 2a + d) * (d + 2a * n)

    iszero(den) && return B̃n(P, Val(n))

    val = one(S) * num / den

    val
end

function C̃n(B::Type{<:AbstractCDOPBasis}, a, b, c, d, e, n::Int)
    S = Int#eltype(P)

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
#αn(P::Type{<:AbstractCDOPPolynomial}, n::Int) = α̃(P,n) / k1k0(P,n)
function α̃n(B::Type{<:AbstractCDOPBasis}, n::Int)
    a, b, c, d, e = abcde(B)
    S = Int#eltype(P)

    num = a * n
    val = one(S) * num

    return val
end

function β̃n(B::Type{<:AbstractCDOPBasis}, n::Int)
    a, b, c, d, e = abcde(B)
    S = Int#eltype(P)

    num =
        -n *
        (d + a * n - a) *
        (2 * a * n * d - a * d - d * b + 2e * a - 2a^2 * n + 2a^2 * n^2)
    den = (2a * n - 2a + d) * (d + 2a * n)

    iszero(den) && return β̃n(P, Val(n))

    val = (one(S) * num) / den

    return val
end
β̃n(B::Type{<:AbstractCDOPBasis}, n::Val) = throw(ArgumentError("Not defined"))

#γn(B::Type{<:AbstractCDOPBasis}, n::Int) =  γ̃n(B,n) * k1k0(B,n-1)
function γ̃n(B::Type{<:AbstractCDOPBasis}, n::Int)
    a, b, c, d, e = abcde(B)
    S = Int #eltype(P)

    num = (n - 1) * (d + a * n - a)
    num *= (a * n * d - d * b - a * d + a^2 * n^2 - 2a^2 * n + 4c * a + a^2 + 2e * a - b^2)
    num += -d * b * e + d^2 * c + a * e^2
    num *= (d + a * n - a) * (a * n + d - 2a) * n
    den = (d - a + 2a * n) * (d + 2a * n - 3a) * (2a * n - 2a + d)^2
    iszero(den) && return γ̃n(P, Val(n))

    val = one(S) * num / den

    return val
end
γ̃n(B::Type{<:AbstractCDOPBasis}, n::Val) = throw(ArgumentError("Not defined"))

function abcdeᴵ(B::Type{<:AbstractCDOPBasis})
    a, b, c, d, e = abcde(B).a, abcde(B).b, abcde(B).c, abcde(B).d, abcde(B).e
    NamedTuple{(:a, :b, :c, :d, :e)}((a, b, c, d + 2a, d + e + a + b))
end

function ẫn(B::Type{<:AbstractCDOPBasis}, n::Int)
    a, b, c, d, e = abcde(B)
    S = Int #eltype(P)

    val = one(S)
    val /= n + 1
    return val
end

function b̂̃n(B::Type{<:AbstractCDOPBasis}, n::Int)
    a, b, c, d, e = abcde(B)
    S = Int#eltype(P)

    num = -2 * a * n * (d + a * n - a) - d * b + a * d - d^2 + 2e * a
    den = (2a * n - 2a + d) * (d + 2a * n)
    iszero(den) && return b̂̃n(B, Val(n))

    val = one(S) * num / den
    return val
end

function ĉ̃n(B::Type{<:AbstractCDOPBasis}, n::Int)
    a, b, c, d, e = abcde(B)
    S = Int#eltype(P)

    num =
        (n - 1) *
        (d + a * n - a) *
        (a * n * d - d * b - a * d + a^2 * n^2 - 2a^2 * n + 4c * a + a^2 + 2e * a - b^2)
    num += -d * b * e + d^2 * c + a * e^2
    num *= a * n
    den = (d - a + 2a * n) * (d + 2a * n - 3a) * (2a * n - 2a + d)^2
    iszero(den) && return ĉ̃n(B, Val(n))

    val = (one(S) * num) / den
    return val
end

for P in Polynomials.ZeroBasedDensePolynomialContainerTypes
    @eval begin
        function Base.:*(
            p::P,
            q::Q,
        ) where {B<:AbstractCDOPBasis,X,T,P<:$P{B,T,X},S,Q<:$P{B,S,X}}
            p′, q′ = _convert_cop.(FallingFactorial, (p, q))
            _convert_cop(⟒(P), p′ * q′)
        end
    end
end

# function ⊗(::Type{<:AbstractCDOPPolynomial}, p::P, q::Q) where {P<:AbstractCDOPPolynomial,Q<:AbstractCDOPPolynomial}

#     #@assert  ⟒(P) == ⟒(Q)
#     #@assert eltype(p) == eltype(q)
#     isconstant(p) && return q * constantterm(p)
#     isconstant(q) && return p * constantterm(q)
#     assert_same_variable(p, q) || throw(ArgumentError("`p` and `q` have different indeterminate"))

#     convert(⟒(P), convert(FallingFactorial, p) * convert(FallingFactorial, q))

#     #    R = eltype(p)
#     #    convert(⟒(P){R}, convert(Polynomial, p) * convert(Polynomial, q))
# end
