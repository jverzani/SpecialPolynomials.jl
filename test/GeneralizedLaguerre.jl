@testset "Convert" begin

    for (α, β) in ((1/2, 1/4), (1/4, 1/2), (1,2), (2,1))
        P, Q = GeneralizedLaguerre{α, Float64}, GeneralizedLaguerre{β, Float64}
        for _ in 1:10
            ps = rand(1.0:8, 5)
            p = P(ps)
            @test -1 == degree(convert(P, convert(Q, p)) - p)
        end
    end

end
