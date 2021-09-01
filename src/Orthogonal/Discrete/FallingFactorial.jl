# represent the   basis x^(mbar) = x⋅(x-1)⋯(x-m+1) ∈ Πm
"""
    FallingFactorial{T}

Construct  a  polynomial with   respect to the basis `x⁰̲,  x¹̲, x²̲, …` where 
`xⁱ̲ = x  ⋅  (x-1) ⋅  (x-2)  ⋯ (x-i+1)` is the falling Pochhammer  symbol.  See 
[Falling factorial](https://en.wikipedia.org/wiki/Falling_and_rising_factorials)  for several  facts
about this  polynomial basis.

In [Koepf and Schmersau](https://arxiv.org/pdf/math/9703217.pdf)
connection coefficients between the falling factorial polynomial
system and classical discrete orthogonal polynomials are given.

## Examples

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> p = basis(FallingFactorial, 3)
FallingFactorial(1.0⋅x³̲)

julia> x = variable(Polynomial)
Polynomials.Polynomial(1.0*x)

julia> p(x) ≈ x*(x-1)*(x-2)
true
```

"""
struct FallingFactorial{T<:Number,X} <: AbstractSpecialPolynomial{T,X}
    coeffs::Vector{T}
    function FallingFactorial{T,X}(coeffs::Vector{T}) where {T<:Number,X}
        length(coeffs) == 0 && return new{T,X}(zeros(T, 1))
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{T,X}(coeffs[1:last])
    end
    function FallingFactorial{T}(
        coeffs::Vector{T},
        var::Polynomials.SymbolLike=:x,
    ) where {T<:Number}
        FallingFactorial{T,Symbol(var)}(coeffs)
    end
end

Polynomials.@register FallingFactorial
export FallingFactorial

# modified from  Polynomials.unicode_exponent
function unicode_exponent(io, var, j)
    print(io, var)
    a = ("⁻", "", "", "⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹")
    for i in string(j)
        print(io, a[Int(i) - 44])
    end
end

# show as 1⋅(x)₀ + 2⋅(x)₁ + 3⋅(x)₂
function Polynomials.showterm(
    io::IO,
    ::Type{P},
    pj::T,
    var,
    j,
    first::Bool,
    mimetype,
) where {N,T,P<:FallingFactorial}
    iszero(pj) && return false
    !first && print(io, " ")
    print(io, Polynomials.hasneg(T) && Polynomials.isneg(pj) ? "- " : (!first ? "+ " : ""))
    print(io, "$(abs(pj))⋅")
    #    print(io,"$(var)")
    unicode_exponent(io, var, j)
    print(io, "̲")
    return true
end

function Polynomials.evalpoly(x, p::FallingFactorial)
    d = degree(p)
    d <= 0 && return p[0] * one(x)

    xⁱ̲ = one(x)
    tot = p[0] * xⁱ̲
    for i in 1:d
        xⁱ̲ *= x - (i - 1)
        tot = muladd(xⁱ̲, p[i], tot)
    end

    tot
end

function Base.one(::Type{P}) where {P<:FallingFactorial}
    T, X = eltype(P), Polynomials.indeterminate(P)
    ⟒(P){T,X}(ones(T, 1))
end

function Polynomials.variable(::Type{P}) where {P<:FallingFactorial}
    T, X = eltype(P), Polynomials.indeterminate(P)
    ⟒(P){T,X}([zero(T), one(T)])
end

k0(P::Type{<:FallingFactorial}) = one(eltype(P))

Base.:*(p::FallingFactorial, q::FallingFactorial) =
    convert(FallingFactorial, convert(Polynomial, p) * convert(Polynomial, q))

## Connect FallingFactorial with AbstractCDOP
function Base.convert(::Type{P}, q::Q) where {P<:FallingFactorial,Q<:AbstractCDOP}
    _convert_cop(P, q)
end
function Base.convert(::Type{P}, q::Q) where {P<:AbstractCDOP,Q<:FallingFactorial}
    _convert_cop(P, q)
end

##  Koepf and Schmersa Thm 6 connection for P, Q=x^n̲
function connection_m(
    ::Type{P},
    ::Type{Q},
    m,
    n,
) where {P<:AbstractCDOP,Q<:FallingFactorial}
    a, b, c, d, e = abcde(P)

    c₀ = (a * n + a * m - a + d) * (n - m)

    c₁ = (m + 1) * (a * n^2 - 2a * m^2 - a * n - a * m + n * d - 2d * m - b * m - d - e)

    c₂ = -(m + 1) * (m + 2) * (a * m^2 + 2a * m + d * m + b * m + a + d + b + c + e)

    (c₀, c₁, c₂)
end

##  Koepf and Schmersa Thm 7, p25 connection for P=x^n̲, Q
function connection_m(
    ::Type{P},
    ::Type{Q},
    m,
    n,
) where {P<:FallingFactorial,Q<:AbstractCDOP}
    ā, b̄, c̄, d̄, ē = abcde(Q)

    c₀ =
        (2m * ā + ā + d̄) *
        (2m * ā + 3ā + d̄) *
        (2m * ā + 2ā + d̄)^2 *
        (2m * ā + d̄) *
        (n - m)

    c₁ = (2m * ā + ā + d̄) * (2m * ā + 3ā + d̄) * (2m * ā + 2ā + d̄) * (m + 1)
    c₁a = 2m^2 * n * ā^2 - 2m^2 * ā^2
    c₁a +=
        m^2 * ā * d̄ + 2m^2 * ā * b̄ + 2m * n * ā^2 + 2m * n * ā * d̄ - 2m * ā^2 -
        m * ā * d̄ +
        2m * ā * b̄ +
        m * d̄^2
    c₁a +=
        2m * d̄ * b̄ + n * ā * d̄ + 2n * ā * ē - n * d̄ * b̄ - ā * d̄ +
        d̄ * b̄ +
        d̄ * ē
    c₁ *= c₁a

    c₂ = (m + 1) * (2m * ā + d̄)
    c₂a = m^4 * ā^3 + 4m^3 * ā^3 + 2m^3 * ā^2 * d̄ + 6m^2 * ā^3 + 6m * ā^2 * d̄
    c₂a +=
        4m^2 * ā^2 * c̄ + 2m^2 * ā^2 * ē + m^2 * ā * d̄^2 - m^2 * ā * d̄ * b̄ -
        m^2 * ā * b̄^2 +
        4m * ā^3 +
        6m^2 * ā^2 * d̄
    c₂a +=
        8m * ā^2 * c̄ + 4m * ā^2 * ē + 2m * ā * d̄^2 - 2m * ā * d̄ * b̄ +
        4m * ā * d̄ * c̄ +
        2m * ā * d̄ * ē
    c₂a +=
        -2m * ā * b̄^2 - m * d̄^2 * b̄ - m * d̄ * b̄^2 +
        ā^3 +
        2ā^2 * d̄ +
        4ā^2 * c̄ +
        2ā^2 * ē +
        ā * d̄^2
    c₂a +=
        -ā * d̄ * b̄ + 4ā * d̄ * c̄ + 2 * ā * d̄ * ē - ā * b̄^2 + ā * ē^2 -
        d̄^2 * b̄ + d̄^2 * c̄
    c₂a += -d̄ * b̄^2 - d̄ * b̄ * ē
    c₂ *= c₂a * (m + 2) * (m * ā + n * ā + ā + d̄)

    (c₀, c₁, c₂)
end

function Base.convert(
    ::Type{P},
    q::Q,
) where {P<:Polynomials.StandardBasisPolynomial,Q<:FallingFactorial}
    q(variable(P))
end
function Base.convert(
    ::Type{P},
    q::Q,
) where {P<:FallingFactorial,Q<:Polynomials.StandardBasisPolynomial}
    connection(P, q)
end

# stirling of 2nd kind
@memoize function sterling²(n::Int, k::Int)
    (iszero(n) && iszero(k)) && return 1
    (iszero(n) || iszero(k)) && return 0

    k * sterling²(n - 1, k) + sterling²(n - 1, k - 1)
end

# xʲ = ∑ᵢ sterling²(j,i) (x)ᵢ [wikipedia](https://en.wikipedia.org/wiki/Falling_and_rising_factorials)
function Base.iterate(
    o::Connection{P,Q},
    state=nothing,
) where {P<:FallingFactorial,Q<:Polynomials.StandardBasisPolynomial}
    n, k = o.n, o.k
    if state == nothing
        j = k
    else
        j = state
        j += 1
    end
    j > n && return nothing
    val = sterling²(j, k)
    (j, val), j
end
