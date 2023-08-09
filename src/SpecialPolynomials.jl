module SpecialPolynomials

using LinearAlgebra
#import LinearAlgebra: dot, ⋅
import SpecialFunctions: gamma

using Polynomials
import Polynomials:
    basis, isconstant, constantterm, assert_same_variable, StandardBasisPolynomial, ⟒, constructorof,
    AbstractBasis, AbstractDenseUnivariatePolynomial,
    MutableDensePolynomial, ImmutableDensePolynomial, MutableDenseViewPolynomial
export basis

#import Intervals
#import Intervals: Open, Closed, Unbounded, bounds_types
import Polynomials: domain, Interval, Open, Closed, Unbounded, bounds_types
using QuadGK
using Memoize

using HypergeometricFunctions

include("utils.jl")
include("abstract.jl")

include("Orthogonal/abstract.jl")
include("Orthogonal/orthogonal.jl")
include("Orthogonal/cop.jl")
include("Orthogonal/ccop.jl")
include("Orthogonal/Bessel.jl")
#=
include("Orthogonal/Chebyshev.jl")
include("Orthogonal/Hermite.jl")
include("Orthogonal/Laguerre.jl")
include("Orthogonal/Gegenbauer.jl")
include("Orthogonal/Jacobi.jl")
include("Orthogonal/Legendre.jl")
include("Orthogonal/WeightFunction.jl")

include("Orthogonal/Discrete/discrete-orthogonal.jl")
include("Orthogonal/Discrete/cdop.jl")
include("Orthogonal/Discrete/FallingFactorial.jl")
include("Orthogonal/Discrete/Hahn.jl")
include("Orthogonal/Discrete/Meixner.jl")
include("Orthogonal/Discrete/Krawchouk.jl")
include("Orthogonal/Discrete/Charlier.jl")
include("Orthogonal/Discrete/DiscreteChebyshev.jl")

include("Orthogonal/connection.jl")
include("Orthogonal/glaser-liu-rokhlin.jl")
include("Orthogonal/comrade.jl")

include("Interpolating/interpolating.jl")
include("Interpolating/Lagrange.jl")
include("Interpolating/Newton.jl")

include("Bernstein.jl")

using Requires
function __init__()
    @static if !isdefined(Base, :get_extension)
        @require FastTransforms = "057dd010-8810-581a-b7be-e3fc3b93f78c" include(
            "fasttransforms.jl",
        )
        @require FastGaussQuadrature = "442a2c76-b920-505d-bb47-c5924d526838" include(
            "fastgaussquadrature.jl",
        )
    end
end
=#
end # module
