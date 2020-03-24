using SpecialPolynomials; const SP=SpecialPolynomials
using Polynomials
import LinearAlgebra: eigvals
using Test

@testset "Orthogonal" begin include("Orthogonal.jl") end
@testset "Orthogonal tables" begin include("Orthogonal-compare.jl") end
@testset "ChebyshevT" begin include("ChebyshevT.jl") end
@testset "ChebyshevU" begin include("ChebyshevU.jl") end

@testset "Bernstein" begin include("Bernstein.jl") end
