module SpecialPolynomials

import LinearAlgebra: dot
import SpecialFunctions: gamma

using Polynomials
using QuadGK
using Memoize
using LinearAlgebra

include("utils.jl")
include("abstract.jl")

include("Orthogonal/orthogonal.jl")
include("Orthogonal/glaser-liu-rokhlin.jl")

include("Orthogonal/Bessel.jl")
include("Orthogonal/ChebyshevT.jl")
include("Orthogonal/ChebyshevU.jl")
include("Orthogonal/Hermite.jl")
include("Orthogonal/Laguerre.jl")
include("Orthogonal/Gegenbauer.jl")
include("Orthogonal/Jacobi.jl")
include("Orthogonal/Legendre.jl")
include("Orthogonal/WeightFunction.jl")

include("Orthogonal/ccop.jl")

include("Orthogonal/Discrete/discrete-orthogonal.jl")
include("Orthogonal/Discrete/DiscreteChebyshev.jl")
include("Orthogonal/Discrete/Krawtchouk.jl")



include("Interpolating/interpolating.jl")
include("Interpolating/Lagrange.jl")
include("Interpolating/Newton.jl")

include("Bernstein.jl")

end # module
