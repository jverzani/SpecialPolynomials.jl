module SpecialPolynomials

import LinearAlgebra: dot
import SpecialFunctions: gamma

using Polynomials
using QuadGK
using Memoize
using LinearAlgebra

include("abstract.jl")

include("Orthogonal/orthogonal.jl")
include("Orthogonal/ChebyshevTT.jl")
include("Orthogonal/ChebyshevU.jl")
include("Orthogonal/Gegenbauer.jl")
include("Orthogonal/GeneralizedLaguerre.jl")
include("Orthogonal/Hermite.jl")
include("Orthogonal/HermiteProb.jl")
include("Orthogonal/Jacobi.jl")
include("Orthogonal/Laguerre.jl")
include("Orthogonal/Legendre.jl")
include("Orthogonal/ShiftedLegendre.jl")
include("Orthogonal/WeightFunction.jl")


include("Orthogonal/Discrete/discrete-orthogonal.jl")
include("Orthogonal/Discrete/DiscreteChebyshev.jl")
include("Orthogonal/Discrete/Krawtchouk.jl")



include("Interpolating/interpolating.jl")
include("Interpolating/Lagrange.jl")
include("Interpolating/Newton.jl")

include("Bernstein.jl")

end # module
