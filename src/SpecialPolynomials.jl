module SpecialPolynomials

import LinearAlgebra: dot
import SpecialFunctions: gamma

using Polynomials
using QuadGK

include("abstract.jl")

include("Orthogonal/orthogonal.jl")
include("Orthogonal/ChebyshevTT.jl")
include("Orthogonal/ChebyshevU.jl")
include("Orthogonal/Hermite.jl")
include("Orthogonal/Jacobi.jl")
include("Orthogonal/Laguerre.jl")
include("Orthogonal/Legendre.jl")

include("Bernstein.jl")
end # module
