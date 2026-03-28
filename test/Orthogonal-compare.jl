# Compare with tabulated values,  classical hypergeometric formulas
using SpecialPolynomials: monic

@testset "Bessel" begin
    T = Float64

    @testset for α in (3 / 2, 1 / 2, 1, 2)
        P = Bessel{α,T}
        β = 2

        p0, p1, p2, p3, p4 = basis.(P, 0:4)

        # Krall and Frink
        @testset for x in range(0, 1, length=5)
            @test p0(x) ≈ 1
            @test p1(x) ≈ 1 + α * (x / β)
            @test p2(x) ≈ 1 + 2(α + 1) * (x / β) + (α + 1) * (α + 2) * (x / β)^2
            @test p3(x) ≈
                  1 +
                  3(α + 2) * (x / β) +
                  3(α + 2) * (α + 3) * (x / β)^2 +
                  (α + 2) * (α + 3) * (α + 4) * (x / β)^3
            @test p4(x) ≈
                  1 +
                  4(α + 3) * (x / β) +
                  6(α + 3) * (α + 4) * (x / β)^2 +
                  4(α + 3) * (α + 4) * (α + 5) * (x / β)^3 +
                  (α + 3) * (α + 4) * (α + 5) * (α + 6) * (x / β)^4
        end
    end

    x = variable(Polynomial)
    @testset for α in (1 / 4, 1 / 2, 3 / 4, 1, 2)
        P = Bessel{α}
        @testset for n in 0:5
            p = basis(P, n)
            @test p(x) ≈ SP.classical_hypergeometric(Polynomials.basistype(P), n, x)
        end
    end
end

@testset "ChebyshevT" begin
    p0 = Chebyshev([1])  # 1
    p1 = Chebyshev([0, 1]) # x
    p2 = Chebyshev([0, 0, 1]) #
    p3 = Chebyshev([0, 0, 0, 1]) #
    p4 = Chebyshev([0, 0, 0, 0, 1]) #

    @testset for x in range(-1, 1, length=5)
        @test p0(x) ≈ 1
        @test p1(x) ≈ x
        @test p2(x) ≈ 2x^2 - 1
        @test p3(x) ≈ 4x^3 - 3x
        @test p4(x) ≈ 8x^4 - 8x^2 + 1
    end
end

@testset "ChebyshevU" begin
    p0 = ChebyshevU([1])  # 1
    p1 = ChebyshevU([0, 1]) # x
    p2 = ChebyshevU([0, 0, 1]) #
    p3 = ChebyshevU([0, 0, 0, 1]) #
    p4 = ChebyshevU([0, 0, 0, 0, 1]) #

    @testset for x in range(-1, 1, length=5)
        @test p0(x) ≈ 1
        @test p1(x) ≈ 2x
        @test p2(x) ≈ 4x^2 - 1
        @test p3(x) ≈ 8x^3 - 4x
        @test p4(x) ≈ 16x^4 - 12x^2 + 1
    end
end

@testset "Hermite" begin
    p0 = Hermite([1])  # 1
    p1 = Hermite([0, 1]) # x
    p2 = Hermite([0, 0, 1]) #
    p3 = Hermite([0, 0, 0, 1]) #
    p4 = Hermite([0, 0, 0, 0, 1]) #

    @testset for x in range(-1, 1, length=5)
        @test p0(x) ≈ 1
        @test p1(x) ≈ 2x
        @test p2(x) ≈ 4x^2 - 2
        @test p3(x) ≈ 8x^3 - 12x
        @test p4(x) ≈ 16x^4 - 48x^2 + 12
    end

    x = variable(LaurentPolynomial)
    P = Hermite
    @testset for n in 2:5
        p = basis(P, n)
        @test p(x) ≈ SP.classical_hypergeometric(Polynomials.basistype(P), n, x)
    end
end

@testset "ChebyshevHermite" begin
    p0 = ChebyshevHermite([1])  # 1
    p1 = ChebyshevHermite([0, 1]) # x
    p2 = ChebyshevHermite([0, 0, 1]) #
    p3 = ChebyshevHermite([0, 0, 0, 1]) #
    p4 = ChebyshevHermite([0, 0, 0, 0, 1]) #

    @testset for x in range(-1, 1, length=5)
        @test p0(x) ≈ 1
        @test p1(x) ≈ x
        @test p2(x) ≈ x^2 - 1
        @test p3(x) ≈ x^3 - 3x
        @test p4(x) ≈ x^4 - 6x^2 + 3
    end

    x = variable(LaurentPolynomial)
    P = ChebyshevHermite
    @testset for n in 2:5
        p = basis(P, n)
        @test p(x) ≈ SP.classical_hypergeometric(Polynomials.basistype(P), n, x) atol = 1e-8
    end
end

@testset "Legendre" begin
    p0 = Legendre([1])  # 1
    p1 = Legendre([0, 1]) # x
    p2 = Legendre([0, 0, 1]) #
    p3 = Legendre([0, 0, 0, 1]) #
    p4 = Legendre([0, 0, 0, 0, 1]) #

    @testset for x in range(-1, stop=1, length=5)
        @test p0(x) ≈ 1
        @test p1(x) ≈ x
        @test p2(x) ≈ 1 / 2 * (3x^2 - 1)
        @test p3(x) ≈ 1 / 2 * (5x^3 - 3x)
        @test p4(x) ≈ 1 / 8 * (35x^4 - 30x^2 + 3)
    end
end

@testset "Shifted Legendre" begin
    T = Float64
    P = ShiftedLegendre{T}

    p0 = P([1])  # 1
    p1 = P([0, 1]) # x
    p2 = P([0, 0, 1]) #
    p3 = P([0, 0, 0, 1]) #
    p4 = P([0, 0, 0, 0, 1]) #

    @testset for x in range(0, stop=1, length=5)
        @test p0(x) ≈ 1
        @test p1(x) ≈ 2x - 1
        @test p2(x) ≈ 6x^2 - 6x + 1
        @test p3(x) ≈ 20x^3 - 30x^2 + 12x - 1
        @test p4(x) ≈ 70x^4 - 140x^3 + 90x^2 - 20x + 1
    end
end

@testset "Laguerre" begin
    P = Laguerre{0}
    p0 = P([1])  # 1
    p1 = P([0, 1]) # x
    p2 = P([0, 0, 1]) #
    p3 = P([0, 0, 0, 1]) #
    p4 = P([0, 0, 0, 0, 1]) #

    @testset for x in range(0, stop=5, length=5)
        @test p0(x) ≈ 1
        @test p1(x) ≈ -x + 1
        @test p2(x) ≈ 1 / 2 * (x^2 - 4x + 2)
        @test p3(x) ≈ 1 / 6 * (-x^3 + 9x^2 - 18x + 6)
        @test p4(x) ≈ 1 / 24 * (x^4 - 16x^3 + 72x^2 - 96x + 24)
    end

    x = variable(Polynomial)
    @testset for α in (0, 1 / 4, 1 / 2, 3 / 4)
        P = Laguerre{α}
        @testset for n in 2:5
            p = basis(P, n)
            @test p(x) ≈ SP.classical_hypergeometric(Polynomials.basistype(P), n, x)
        end
    end
end

@testset "Gegenbauer" begin
    @testset for α in (1 / 4, 1 / 3, 1 / 2, 1, 2)
        P = Gegenbauer{α,Float64}
        p0 = P([1])  # 1
        p1 = P([0, 1]) # x
        p2 = P([0, 0, 1]) #
        p3 = P([0, 0, 0, 1]) #
        p4 = P([0, 0, 0, 0, 1]) #

        @testset for x in range(-1, stop=1, length=5)
            @test p0(x) ≈ 1
            @test p1(x) ≈ 2α * x
            @test p2(x) ≈ -α + 2α * (1 + α) * x^2
            @test p3(x) ≈ -2α * (1 + α) * x + 4 / 3 * α * (1 + α) * (2 + α) * x^3
        end
    end

    x = variable(Polynomial)
    @testset for α in (1 / 4, 1 / 2, 3 / 4, 1)
        P = Gegenbauer{α}
        @testset for n in 2:5
            p = basis(P, n)
            @test p(x) ≈ SP.classical_hypergeometric(Polynomials.basistype(P), n, x)
        end
    end

    # Jacobi{α - 1/2, α-1/2 } = Gegenbauer{α}
    # Legendre = Gegenbauer{1/2}
    # ChebyshevU = Gegenbauer{1}
    @testset for (alpha, Q) in
                 ((1/4, Jacobi{1/4 - 1/2,1/4 - 1/2}), (1, ChebyshevU), (1/2, Legendre))
        P = Gegenbauer{alpha}
        @testset for i in 2:6
            p = basis(P, i)
            q = basis(Q, i)
            monic(p)(x) ≈ monic(q)(x)
        end
    end
end

@testset "Jacobi" begin
    T = Float64
    x = variable()

    # ChebyshevT = J(-1/2, -1/2)
    # ChebyshevU = J(1/2, 1/2)
    # Legendre = J(0,0)
    @testset for (alpha_beta, Q) in (
        ((-1 / 2, -1 / 2), Chebyshev),
        ((1 / 2, 1 / 2), ChebyshevU),
        ((0, 0), Legendre),
    )
        P = Jacobi{alpha_beta...}
        @testset for i in 2:6
            p = basis(P, i)
            q = basis(Q, i)
            monic(p)(x) ≈ monic(q)(x)
        end
    end

    x = variable(Polynomial)
    @testset for (α, β) in
                 ((1 / 2, 1 / 2), (-1 / 2, 1 / 2), (1 / 2, -1 / 2), (-1 / 2, -1 / 2))
        P = Jacobi{α,β}
        @testset for n in 2:5
            p = basis(P, n)
            @test p(x) ≈ SP.classical_hypergeometric(Polynomials.basistype(P), n, x)
        end
    end

    ## We have `jacobi_eval` for α and β outside [-1, ∞)
    α, β = 1 / 2, -1 / 2
    @testset for n in 2:5
        x0 = 1 / 2
        @test SP.jacobi_eval(α, β, n, x0) ≈ basis(Jacobi{α,β}, n)(x0)
    end
end

## CDOP

@testset "Charlier" begin
    x = variable()
    @testset for μ in (1 / 4, 1 / 2, 1, 2)
        P = Charlier{μ}
        @test basis(P, 0)(x) ≈ one(x)
        @test basis(P, 1)(x) ≈ 1 - x / μ
        @test basis(P, 2)(x) ≈ (x^2 + μ^2 - x * (1 + 2μ)) / μ^2
    end

    @testset for μ in (1 / 4, 1 / 2, 1, 2)
        P = Charlier{μ}
        @testset for i in 0:5
            @test basis(P, i)(x) ≈
                  SP.classical_hypergeometric(Polynomials.basistype(P), i, x)
        end
    end
end

@testset "Hahn" begin
    x = variable()
    @testset for α in (1 / 4, 1 / 2)
        @testset for β in (1 / 4, 1 / 2)
            𝐍 = 5

            P = Hahn{α,β,𝐍}
            @testset for i in 0:𝐍
                @test basis(P, i)(x) ≈
                      SP.classical_hypergeometric(Polynomials.basistype(P), i, x)
            end
        end
    end
end

@testset "Krawchouk" begin
    x = variable()
    @testset for p in (1 / 4, 1 / 2, 1, 2)
        @testset for N in (3, 6)
            P = Krawchouk{p,N}
            @test basis(P, 0)(x) ≈ one(x)
            @test basis(P, 1)(x) ≈ -N * p + x
            @test basis(P, 2)(x) ≈ 1 / 2 * (N^2 * p^2 + x * (2p + x - 1) - N * p * (p + 2x))
        end
    end

    @testset for p in (1 / 4, 1 / 2, 1, 2)
        @testset for N in (3, 6)
            P = Krawchouk{p,N}
            @testset for i in 0:N
                @test basis(P, i)(x) ≈
                      SP.classical_hypergeometric(Polynomials.basistype(P), i, x)
            end
        end
    end
end

@testset "Meixner" begin
    x = variable()
    @testset for γ in (1 / 4, 1 / 2, 3 / 4)
        @testset for μ in (1 / 4, 1 / 2, 3 / 4)
            P = Meixner{γ,μ}
            @testset for i in 0:4
                @test basis(P, i)(x) ≈
                      SP.classical_hypergeometric(Polynomials.basistype(P), i, x)
            end
        end
    end
end
