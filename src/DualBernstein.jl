#=
struct DualBernstein{𝐍,α,β,T,X} <: AbstractBernstein{T,X}
    coeffs::Vector{T}
    function DualBernstein{𝐍,α,β,T,X}(coeffs::AbstractVector{T}) where {𝐍,α, β,T,X}
        # coeffs is length 𝐍 + 1 make sure
        N = length(coeffs)
        (N > 𝐍 + 1) && throw(ArgumentError("Wrong length for coefficients"))
        if N < 𝐍
            coeffs = vcat(coeffs, zeros(𝐍 - N + 1))
        end
        return new{𝐍,α,β,T,X}(coeffs)
    end
end
function DualBernstein{𝐍,α,β}(coeffs::AbstractVector{T}; var=:x) where {𝐍,α, β,T}
        N = findlast(!iscoeffzero, coeffs)
        N == nothing && return new{𝐍,α,β,T,X}(zeros(T, 0))
        (N > 𝐍 + 1) && throw(ArgumentError("Wrong length for coefficients"))
        return DualBernstein{𝐍,α,β,T,var}(coeffs[1:N])
    end

=#
struct DualBernsteinBasis{𝐍,α,β} <: AbstractBasis end

"""

    DualBernstein{N, α, β, T, X}

Given the inner product `<f,g> = ∫₀¹ (1-x)ᵅ xᵝ f(x) g(x) dx` and the Bernstein polynomials `Bᵢⁿ(x)` the dual Bernstein polynomials satisfy `<Bᵢⁿ, Dⱼⁿ(x;α,β)> = δᵢⱼ`

## Example

```julia
julia> n,α,β = 5, 0, 0;

julia> D = DualBernstein{n,α,β}
DualBernstein{5, 0, 0}

julia> bi = basis(D, 3)
DualBernstein(1.0⋅βᵅᵝ₅,₃(x))

julia> bi(0.2)
7.910400000000036
```

The Bernstein-Bezier form of `f` minimizes the value of the least-square error for the `α-β` norm.

```julia
julia> f(x) = sinpi(x);

julia> n, α, β = 5, 1/2, 1/2
(5, 0.5, 0.5)

julia> B, D = Bernstein{n}, DualBernstein{n,α,β};

julia> Iₖ = [ip(f, basis(D,k), α, β) for k in 0:n];

julia> pₙ = B(Iₖ);

julia> λ = x -> f(x) - pₙ(x);

julia> SpecialPolynomials.innerproduct(ShiftedJacobi{α, β}, λ, λ)
0.00017514530881540565
```

## Reference

The implementation follows that of [Chudy and Woźny](https://arxiv.org/abs/2004.09801).
"""
DualBernstein{𝐍,α,β} = MutableDensePolynomial{DualBernsteinBasis{𝐍,α,β}} where {𝐍,α,β}
export DualBernstein

Polynomials._typealias(
    ::Type{P},
) where {P<:DualBernstein{𝐍,α,β}} where {𝐍,α,β} = "DualBernsteinᵅᵝ"
export DualBernstein

Polynomials.@registerN DualBernstein 𝐍
Polynomials.:⟒(::Type{<:DualBernstein{𝐍,α,β}}) where {𝐍,α,β} = Bernstein{𝐍,α,β}

function basis(
    P::Type{<:DualBernstein{𝐍,α,β}},
    k::Int,
    _var::Polynomials.SymbolLike=:x;
    var=_var,
) where {𝐍,α,β}
    @assert k <= 𝐍
    T = eltype(P)
    zs = zeros(T, 𝐍 + 1)
    zs[k + 1] = 1
    P{T,var}(zs)
end

function basis(::P, k::Int) where {𝐍,α,β,P<:DualBernstein{𝐍,α,β}}
    basis(P, k)
end

function Polynomials.showterm(
    io::IO,
    ::Type{DualBernstein{N,α,β,T,X}},
    pj::T,
    var,
    j,
    first::Bool,
    mimetype,
) where {N,α,β,T,X}
    iscoeffzero(pj) && return false
    !first && print(io, " ")
    print(io, Polynomials.hasneg(T) && Polynomials.isneg(pj) ? "- " : (!first ? "+ " : ""))
    print(io, "$(abs.(pj))⋅β")
    print(io, "ᵅᵝ")

    Polynomials.unicode_subscript(io, N)
    print(io, ",")
    Polynomials.unicode_subscript(io, j)
    print(io, "($var)")
    return true
end

# XXX convert between families
Polynomials.domain(::Type{<:DualBernstein}) = Polynomials.Interval(0, 1)
#XXXPolynomials.constantterm(p::DualBernstein) = p[0] # β.,ᵢ = 0 + ... when i > 0
Base.one(P::Type{DualBernstein{N,α,β,T}}, var::Polynomials.SymbolLike) where {N,α,β,T} =
    DualBernstein{N,α,β,T,Symbol(var)}(ones(T, N + 1))
function Base.one(
    P::Type{<:DualBernstein{N,α,β}},
    var::Polynomials.SymbolLike=:x,
) where {N,α,β}
    N′ = max(N, 0)
    one(DualBernstein{N′,α,β,eltype(P),Symbol(var)})
end

#=
XXXfunction Polynomials.variable(P::Type{DualBernstein{𝐍,α,β,T,X}}) where {𝐍,α,β,T,X}
    𝐍 <= 0 && throw(ArgumentError("Need  𝐍 ≥ 1"))
    R = typeof(one(T) / 1)
    DualBernstein{𝐍,α,β,R,X}([1 - (i * one(T)) / 𝐍 for i in 𝐍:-1:0])
end
=#

function classical_hypergeometric(::Type{<:DualBernstein{n,α,β}}, i, x) where {n,α,β}
    (α ≤ -1 || β ≤ -1) && throw(ArgumentError("α and β must be > -1"))

    S(m, k, a, b, z) =
        Pochhammer(b+1, m) * sum(
            Pochhammer(-m, j)*Pochhammer(-m-a, j)/factorial(j)/Pochhammer(b+1, j)*z^j for
            j in 0:k;
            init=0.0,
        )

    σ = α + β + 1
    K = gamma(α+1)*gamma(β+1)/gamma(σ+1)

    # cf. Corollary 2.2
    out = inv(binomial(n, i))
    out *= (-1)^(n-i) * Pochhammer(σ+1, n)
    out /= K
    out /= Pochhammer(α+1, n)
    out /= Pochhammer(β+1, n)
    out *= ((x-1)/x)^i
    tmp = basis(ShiftedJacobi{α,β+1}, n)(x) * S(n, i, α+1, β, x/(x-1))
    tmp -= basis(ShiftedJacobi{α+1,β}, n)(x) * S(n, i-1, α, β+1, x/(x-1))
    out *= tmp
    out
end

# we only have evaluation
Polynomials.evalpoly(x, p::DualBernstein) = ChudyWozny_eval(p, x)

# We use Algorithms 1 and 2 of Chudy and Wozny to evaluate
# a polynomial
function ChudyWozny_eval(p::DualBernstein{𝐍,α,β}, x) where {𝐍,α,β}
    D = AllDualBer(𝐍, α, β, x)
    coeffs = p.coeffs
    sum(pᵢ*Dᵢ for (pᵢ, Dᵢ) in zip(coeffs, D); init=zero(x))
end

𝐽(n::Int, x::Any) = n ÷ 2
function 𝐽(n::Int, x::Real)
    ps = (
        0.08401156564574855,
        1.62239798468112882,
        −2.37126334791787779,
        1.58084223194525186,
    )
    pₓ = evalpoly(x, ps)
    round(Int, n * pₓ)
end

# their algorithms
function DualBer(n, α, β, x, j, K)
    α1 = α + 1
    β1 = β + 1
    n1 = n + α + 1
    x1x = (x - 1) / x

    C = (iseven(n) ? -1 : 1) * K/n1 * prod(1 + β1/(j+α1) for j in 0:(n - 1); init=1)
    R₁ = n1 * basis(ShiftedJacobi{α,β1}, n)(x)
    R₂ = x1x * (n + β1) * basis(ShiftedJacobi{α1,β}, n)(x)
    Δ = 1 # offset
    D = zeros(typeof(-C*R₁), j + 1)
    D[0 + Δ] = -C * R₁
    for i in 1:j
        p = i - n - 1
        q = i / p
        C = C * (p - α1) / (i + β)
        Dᵢ = q * x1x * D[(i - 1) + Δ] - C * (R₁ + q*R₂)
        D[i + Δ] = Dᵢ
    end
    D
end

function AllDualBer(n, α, β, x)
    α1 = α + 1
    β1 = β + 1
    K = gamma(α1 + β1) / gamma(α1) / gamma(β1)
    J = 𝐽(n, x)
    fwd = DualBer(n, α, β, x, J, K)
    D = zeros(eltype(fwd), n + 1)
    for (i, x) in enumerate(fwd)
        D[i] = x
    end
    if n - J - 1 > 0
        bckwd = DualBer(n, β, α, 1-x, n - J - 1, K)
        for (i, x) in enumerate(bckwd)
            D[end - i + 1] = x
        end
    end
    D
end
