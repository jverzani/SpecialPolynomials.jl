# special case these once FastGaussQuadrature is loaded
gauss_nodes_weights(p::Type{P}, n) where {P<:Chebyshev} =
    gausschebyshev(n)

gauss_nodes_weights(p::Type{P}, n) where {α,β,P<:Jacobi{α,β}} =
    FastGaussQuadrature.gaussjacobi(n, α, β)

gauss_nodes_weights(p::Type{P}, n) where {P<:Legendre} =
    gausslegendre(n)

gauss_nodes_weights(p::Type{P}, n) where {P<:Hermite} =
    gausshermite(n)

gauss_nodes_weights(p::Type{P}, n) where {α, P<:Laguerre{α}} =
    gausslaguerre(n, α)

include("Orthogonal/FastLegendre.jl")

eval_basis(::Type{P}, n, x::Float64) where {P<:Legendre} =
    FastLegendre.fastlegendre(n, x)

eval_basis(::Type{P}, n, x::Union{Int,Int32,Float16,Float32}) where {P<:Legendre} =
    FastLegendre.fastlegendre(n, float(x))
