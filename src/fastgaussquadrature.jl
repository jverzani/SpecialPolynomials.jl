
gauss_nodes_weights(p::Type{P}, n) where {P<:Chebyshev} =
    FastGaussQuadrature.gausschebyshev(n)

gauss_nodes_weights(p::Type{P}, n) where {α,β,P<:Jacobi{α,β}} =
    FastGaussQuadrature.gaussjacobi(n, α, β)

include("Orthogonal/FastLegendre.jl")

eval_basis(::Type{P}, n, x::Float64) where {P<:Legendre} = FastLegendre.fastlegendre(n, x)
eval_basis(::Type{P}, n, x::Union{Int,Int32,Float16,Float32}) where {P<:Legendre} =
    FastLegendre.fastlegendre(n, float(x))
gauss_nodes_weights(p::Type{P}, n) where {P<:Legendre} =
    FastGaussQuadrature.gausslegendre(n)
