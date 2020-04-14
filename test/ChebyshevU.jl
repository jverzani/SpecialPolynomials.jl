
@testset "Conversion" begin

    for _ in 1:10
        ps = rand(1.0:8, 5)
        p = Chebyshev(ps)
        @test -1 == degree(convert(Chebyshev, convert(ChebyshevU, p)) - p)
    end

    for _ in 1:10
        ps = rand(1.0:8, 5)
        p = ChebyshevU(ps)
        @test -1 == degree(convert(ChebyshevU, convert(Chebyshev, p)) - p)
    end

end

@testset "Evaluation" begin

    p0 = ChebyshevU([1])  # 1
    p1 = ChebyshevU([0,1]) # 2x
    p2 = ChebyshevU([0,0,1]) # 4x^2 - 1
    p3 = ChebyshevU([0,0,0,1]) # 8x^3 - 4x
    p4 = ChebyshevU([0,0,0,0,1]) # 16x^4 - 12x^2 + 1

    x = 1/2
    @test p0(x) ≈ 1
    @test p1(x) ≈ 2x
    @test p2(x) ≈ 4x^2 - 1
    @test p3(x) ≈ 8x^3 - 4x
    @test p4(x) ≈ 16x^4 - 12x^2 + 1
    @test ChebyshevU([1,2,3,4])(x) ≈ 1p0(x) + 2p1(x) + 3p2(x) + 4p3(x)


end



@testset "roots" begin

    for n in 4:2:8
        p = ChebyshevU(vcat(zeros(n),1))
        @test cos.((n:-1:1)/(n+1) * pi) ≈ roots(p)
    end

end
