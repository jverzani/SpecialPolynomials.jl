# some utility functions
const Γ = gamma

function generalized_binomial(α::T, n::S) where {T,S <: Integer}
    R = Base.promote_op(/,T,S)
    out = one(R)
    for i in n:-1:1
        out *= α/i
        α -= 1
    end
    out
end



## x⁽ⁱ⁾ is falling (x)⋅(x-1)⋯(x-n+1)
function Pochhammer(::Val{:falling}, z::T,n::Int) where {T}
    iszero(n) && return one(T)    
    out = z
    for i in 1:n-1
        z -= 1
        out *= z
    end
    out
end

## (x)ᵢ is rising (x)⋅(x+1)⋯(x+n-1) = Γ(z+n)/Γ(z)
function Pochhammer(::Val{:rising}, z::T,n::Int) where {T}
    iszero(n) && return one(T)
    out = z
    for i in 1:n-1
        z += 1
        out *= z
    end
    out
end
# default is (x)ᵢ
Pochhammer(z,n::Int) = Pochhammer(Val(:rising),  z, n)

# (x)_n/n!
Pochhammer_factorial(z, n) = Pochhammer_factorial(Val(:rising), z, n)

function Pochhammer_factorial(::Val{:rising}, z, n)
    iszero(n) && return one(z)/1
    prod((z+i)/(n-i) for i in 0:n-1)
end

## Alternative to HypergeometricFunctions.jl so that
## variables can be inserted into the numerator
"""
    pFq(as, bs, z; [maxevals=1000])

Compute the generalized [hypergeometric](https://en.wikipedia.org/wiki/Generalized_hypergeometric_function) function. The `HypergeometricFunction.jl` package would normally be used for such calculations, but this one accepts polynomial values for the `as` and `z` values.

# Example

From [mathworld](https://mathworld.wolfram.com/HypergeometricFunction.html)

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> import SpecialPolynomials: pFq, Pochhammer

julia> pFq((1/3,2/3), 5/6, 27/32) ≈ 8/5
true

julia> pFq([1/4, 1/2], [3/4], 80/81; maxevals=2000) ≈ 9/5
true

julia> x = variable()
Polynomials.Polynomial(x)

julia> n = 5
5

julia> pFq((-n,n+1), 1, (1-x)/2) ≈ basis(Legendre,n)(x)
true

julia> α, β, n = 1/2, 1/2, 5;

julia> Pochhammer(α+1,n)/factorial(n) * pFq((-n,n+α+β+1), α+1, (1-x)/2) ≈ basis(Jacobi{α,β}, n)(x)
true
```
"""
function pFq(as, bs, z; maxevals=1000)
    n = 1
    acc = trm = one(z)

    p,q = length(as), length(bs)

    # element in `as` is negative integer or 0 && converges
    # p ≤ q && converges
    # p = q + 1 |x| < 1 && converge
    # p > q + 1 && diverges
    while n < maxevals
        a = isempty(as) ? 1 : prod(as)
        b = isempty(bs) ? 1 : prod(bs)

        iszero(a) && return acc
        iszero(b) && return Inf  # check on b

        trm *= a /b * z  /n
        acc += trm

        as = plus_1(as)

        bs = plus_1(bs)
        n += 1
    end

    if (p > q + 1 || (p == q+1 && abs(z) < 1))
        return acc
    else
        NaN*acc
    end
end
## tuples
@inline plus_1(as) = map(x->x+1, as)
        
    
## Try to speed up  quadgk by not specializing on F
mutable struct Wrapper
    F
end
(F::Wrapper)(x) = F.F(x)

_quadgk(f, a, b) = quadgk(Wrapper(f), a, b)[1]
const ∫ = _quadgk


checked_div(a, b) = (iszero(a) && iszero(b)) ? a : a/b
