using SpecialPolynomials; const SP=SpecialPolynomials
using Polynomials
import LinearAlgebra: eigvals, det, dot, norm
import SpecialFunctions: gamma
using Test

@testset "Orthogonal" begin include("Orthogonal.jl") end
@testset "Orthogonal compare with tables" begin include("Orthogonal-compare.jl") end
#@testset "ChebyshevT" begin include("ChebyshevT.jl") end
#@testset "ChebyshevU" begin include("ChebyshevU.jl") end
#@testset "Laguerre{α}" begin include("GeneralizedLaguerre.jl") end
#@testset "Gegenbauer{α}" begin include("Gegenbauer.jl") end

@testset "Lagrange" begin include("Lagrange.jl") end
@testset "Newton" begin include("Newton.jl") end

@testset "Bernstein" begin include("Bernstein.jl") end
