module FastGaussQuadratureExt

using SpecialPolynomials
using FastGaussQuadrature

SpecialPolynomials.gauss_nodes_weights(p::Type{P}, n) where {P<:Chebyshev} =
    gausschebyshev(n)

SpecialPolynomials.gauss_nodes_weights(p::Type{P}, n) where {α,β,P<:Jacobi{α,β}} =
    gaussjacobi(n, α, β)

SpecialPolynomials.eval_basis(::Type{P}, n, x::Float64) where {P<:Legendre} = fastlegendre(n, x)
SpecialPolynomials.eval_basis(::Type{P}, n, x::Union{Int,Int32,Float16,Float32}) where {P<:Legendre} =
    fastlegendre(n, float(x))
SpecialPolynomials.gauss_nodes_weights(p::Type{P}, n) where {P<:Legendre} =
    gausslegendre(n)



# From
# O(1) Computation of Legendre polynomials and Gauss-Legendre nodes and weights for parallel computing
# Bogaert, Ignace and Michiels, Bart and Fostier, Jan
# SIAM JOURNAL ON SCIENTIFIC COMPUTING
# 34, 2012
# http://dx.doi.org/10.1137/110855442
#
# as mentioned https://github.com/JuliaMath/SpecialFunctions.jl/issues/124
#
# Only Float64
#
# Gives a speed up over Clenshaw, which is expected, as this is O(1), not O(n):
#
# n= 1000 Clenshaw might be faster
# n= 10_000 10x faster
# n= 100_000 10x faster
# n= 1_000_000 500x faster
# n= 10_000_000 10,000x faster
# (measurement of basis(P,n)(x) compared to Basis(P,n)(x)

import FastGaussQuadrature: gamma

# 4 regions of evaluation
# specialized for Float64 values
function fastlegendre(n, x)
    n == 0 && return one(x)
    n == 1 && return x
    n < 100 && return direct_32(n, x)
    x < 0 && return (iseven(n) ? 1 : -1) * fastlegendre(n, -x)
    θ = acos(x)
    (n + 1) * sin(θ) > 25 && return asy_32(n, x)
    asy_33(n, x)
end

## ---- n < 100 case

# precompute first 99
const xᵢs = [gausslegendre(n)[1][1:(n ÷ 2)] for n in 1:99]
const xᵢs⁻¹ = [inv.(xᵢ) for xᵢ in xᵢs]
const xᵢs2⁻¹ = [inv.(1 .- xᵢ .^ 2) for xᵢ in xᵢs]

# precompute
const Cn = zeros(Float64, 99)
Cn[1] = 1.0
for i in 1:49
    n = 2i
    Cn[n] = sqrt(pi) / gamma(i + 1) / gamma((1 - n) / 2)
    Cn[n + 1] = (n + 1) * Cn[n]
end

# use cosine or sine
function direct_32(n, x)
    if abs(x) > 0.707
        if iseven(n)
            direct_sine(Val(true), n, x)
        else
            direct_sine(Val(false), n, x)
        end
    else
        if iseven(n)
            direct_cosine(Val(true), n, x)
        else
            direct_cosine(Val(false), n, x)
        end
    end
end

#x = cos(theta) 1 - x^2 = sintheta0
function direct_sine(::Val{true}, n, x)
    tot = one(Float64)
    xs = xᵢs2⁻¹[n]
    for aᵢ in xs
        tot *= 1 - (1 - x^2) * aᵢ
    end
    tot
end

direct_sine(::Val{false}, n, x) = x * direct_sine(Val(true), n, x)

function direct_cosine(p::Val{true}, n, x)
    tot = Cn[n]
    xs = xᵢs⁻¹[n]
    for aᵢ in xs
        tot *= (1 - (x * aᵢ)^2)
    end
    tot
end
direct_cosine(p::Val{false}, n, x) = x * direct_cosine(Val(true), n, x)

## ----
const Cl0s = [gamma(n + 1) / gamma(n + 3 / 2) for n in 1:9]
τ(x) = evalpoly(
    1 / x,
    (
        1 / 1,
        0,
        -1 / 64,
        0,
        21 / 8192,
        0,
        -671 / 524288,
        0,
        180323 / 80323134,
        0,
        -20898423 / 8589934592,
        -7426362705 / 1099511627776,
    ),
)

function Cl0(l)
    1 <= l < 10 && return Cl0s[l]
    1 / sqrt(l + 3 / 4) * τ(l + 3 / 4)
end

function asy_32(n, x)
    M = 17
    θ = acos(x)
    val = sqrt((2 / pi) * (1 / sin(θ)))
    tot = 0.0
    sθ = sin(θ)
    λ = Cl0(n)
    for m in 0:(M - 1)
        α = (n + m + 1 / 2) * θ - (m + 1 / 2) * pi / 2
        tot += λ * cos(α)
        λ *= (m + 1 / 2)^2 / 2 / (m + 1) / (n + m + 3 / 2) / sθ
    end
    val * tot
end

## ---
h(n, y) = y^n * besselj(n, y)

const cns = [
    [1 / 8, -1 / 12],
    [11 / 384, -7 / 160, 1 / 160],
    [173 / 15360, -101 / 3584, 671 / 80640, -61 / 120960],
    [22931 / 3440640, -90497 / 3870720, 217 / 20480, -1261 / 967680, 1261 / 29030400],
    [
        1319183 / 247726080,
        -10918993 / 454164480,
        +1676287 / 113541120,
        -7034857 / 2554675200,
        +1501 / 8110080,
        -79 / 20275200,
    ],
    [
        233526463 / 43599790080,
        -1396004969 / 47233105920,
        +2323237523 / 101213798400,
        -72836747 / 12651724800,
        +3135577 / 5367398400,
        -1532789 / 61993451520,
        +66643 / 185980354560,
    ],
]

function fn(i, y)
    n = 2i
    cs = cns[i]
    tot = 0.0
    Δ = i - 1
    for j in i:(2i)
        tot += cs[j - Δ] * h(j, y)
    end
    tot
end

function asy_33(n, x)
    v = n + 1 / 2
    θ = acos(x)
    y = v * θ

    # Eq (3.19)
    val = besselj(0, y)

    λ = 1 / v^2
    for i in 1:6
        val += fn(i, y) * λ
        λ /= v^2
    end
    val
end

end
