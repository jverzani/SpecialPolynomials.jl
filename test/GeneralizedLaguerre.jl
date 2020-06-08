@testset "Convert" begin

    x = variable()
    for (α, β) in ((1/2, 1/4), (1/4, 1/2), (1,2), (2,1))
        P, Q = Laguerre{α}, Laguerre{β}
        for _ in 1:10
            p = P(rand(1:5,5)); q=P(rand(1:5,5))
            @test convert(P, convert(Q, p))(x) ≈ p(x)
            @test convert(Q, convert(P,q))(x) ≈ q(x)
        end
    end
    
end
