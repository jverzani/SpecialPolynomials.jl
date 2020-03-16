using SpecialPolynomials
using Polynomials
using Test

@testset "ChebyshevU" begin include("ChebyshevU.jl") end
@testset "Legendre" begin include("Legendre.jl") end
@testset "Bernstein" begin include("Bernstein.jl") end
