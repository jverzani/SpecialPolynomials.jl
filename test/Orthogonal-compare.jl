# Compare with tabulated values

@testset "ChebyshevT" begin

    p0 = Chebyshev([1])  # 1
    p1 = Chebyshev([0,1]) # x
    p2 = Chebyshev([0,0,1]) #
    p3 = Chebyshev([0,0,0,1]) #
    p4 = Chebyshev([0,0,0,0,1]) #

    for x in range(-1, 1, length=5)
        @test p0(x) ≈ 1
        @test p1(x) ≈ x
        @test p2(x) ≈ 2x^2 - 1
        @test p3(x) ≈ 4x^3  - 3x
        @test p4(x) ≈ 8x^4 - 8x^2 + 1
    end

end

@testset "ChebyshevU" begin

    p0 = ChebyshevU([1])  # 1
    p1 = ChebyshevU([0,1]) # x
    p2 = ChebyshevU([0,0,1]) #
    p3 = ChebyshevU([0,0,0,1]) #
    p4 = ChebyshevU([0,0,0,0,1]) #

    for x in range(-1, 1, length=5)
        @test p0(x) ≈ 1
        @test p1(x) ≈ 2x
        @test p2(x) ≈ 4x^2 - 1
        @test p3(x) ≈ 8x^3 - 4x
        @test p4(x) ≈ 16x^4 - 12x^2 +  1
    end

end


@testset "Hermite" begin

    p0 = Hermite([1])  # 1
    p1 = Hermite([0,1]) # x
    p2 = Hermite([0,0,1]) #
    p3 = Hermite([0,0,0,1]) #
    p4 = Hermite([0,0,0,0,1]) #

    for x in range(-1, 1, length=5)
        @test p0(x) ≈ 1
        @test p1(x) ≈ 2x
        @test p2(x) ≈ 4x^2 - 2
        @test p3(x) ≈ 8x^3 - 12x
        @test p4(x) ≈ 16x^4 - 48x^2 + 12
    end

end

@testset "ChebyshevHermite" begin

    p0 = ChebyshevHermite([1])  # 1
    p1 = ChebyshevHermite([0,1]) # x
    p2 = ChebyshevHermite([0,0,1]) #
    p3 = ChebyshevHermite([0,0,0,1]) #
    p4 = ChebyshevHermite([0,0,0,0,1]) #

    for x in range(-1, 1, length=5)
        @test p0(x) ≈ 1
        @test p1(x) ≈ x
        @test p2(x) ≈ x^2 - 1
        @test p3(x) ≈ x^3 - 3x
        @test p4(x) ≈ x^4 - 6x^2 + 3
    end

end

@testset "Legendre" begin

    p0 = Legendre([1])  # 1
    p1 = Legendre([0,1]) # x
    p2 = Legendre([0,0,1]) #
    p3 = Legendre([0,0,0,1]) #
    p4 = Legendre([0,0,0,0,1]) #

    for x in range(-1, 1, length=5)
        @test p0(x) ≈ 1
        @test p1(x) ≈ x
        @test p2(x) ≈ 1/2*(3x^2-1)
        @test p3(x) ≈ 1/2*(5x^3-3x)
        @test p4(x) ≈ 1/8*(35x^4 - 30x^2 + 3)
    end

end


@testset "Shifted Legendre" begin

    T = Float64
    P = ShiftedLegendre{T}

    p0 = P([1])  # 1
    p1 = P([0,1]) # x
    p2 = P([0,0,1]) #
    p3 = P([0,0,0,1]) #
    p4 = P([0,0,0,0,1]) #

    for x in range(0, 1, length=5)
        @test p0(x) ≈ 1
        @test p1(x) ≈ 2x - 1
        @test p2(x) ≈ 6x^2 - 6x + 1
        @test p3(x) ≈ 20x^3 - 30x^2+ 12x - 1
        @test p4(x) ≈ 70x^4 - 140x^3 + 90x^2 - 20x + 1
    end

end



@testset "Laguerre" begin

    p0 = Laguerre([1])  # 1
    p1 = Laguerre([0,1]) # x
    p2 = Laguerre([0,0,1]) #
    p3 = Laguerre([0,0,0,1]) #
    p4 = Laguerre([0,0,0,0,1]) #

    for x in range(0, 5, length=5)
        @test p0(x) ≈ 1
        @test p1(x) ≈ -x + 1
        @test p2(x) ≈ 1/2*(x^2 - 4x + 2)
        @test p3(x) ≈ 1/6*(-x^3 + 9x^2 -18x + 6)
        @test p4(x) ≈ 1/24 * (x^4 - 16x^3 + 72x^2 - 96x + 24)
    end


end


@testset "Gegenbauer" begin

    for α in (1/4, 1/3, 1/2, 1, 2)
        P = Gegenbauer{α, Float64}
        p0 = P([1])  # 1
        p1 = P([0,1]) # x
        p2 = P([0,0,1]) #
        p3 = P([0,0,0,1]) #
        p4 = P([0,0,0,0,1]) #

        for x in range(-1, 1, length=5)
            @test p0(x) ≈ 1
            @test p1(x) ≈ 2α*x
            @test p2(x) ≈ -α + 2α*(1+α)*x^2
            @test p3(x) ≈ -2α * (1 + α) *x + 4/3*α*(1+α) * (2 + α)*x^3
        end
    end
end

@testset "Jacobi" begin

    T = Float64
    ps  = [0,0,0,0,0,1]
    # ChebyshevT = J(-1/2, -1/2)
    for (alpha_beta, P) in  (((-1/2, -1/2), Chebyshev{T}),
                             ((1/2,1/2), ChebyshevU{T}),
                             ((0,0), Legendre{T})
                             )
        for x  in range(-1, 1, length=6)
            p1 = Jacobi{alpha_beta..., T}(ps)
            q1 = P(ps)
            @test SP._monic(p1)(x) ≈ SP._monic(q1)(x)
        end
    end

end
