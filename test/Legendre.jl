@testset "Evaluation" begin

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
