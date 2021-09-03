@testset "Special cases" begin
    T = Float64
    f(l) = n -> gamma(l + 1 / 2) / gamma(2l) * gamma(n + 2l) / gamma(n + l + 1 / 2)

    for (P, Q, fn) in (
        (Gegenbauer{1 / 2,T}, Legendre{T}, n -> 1.0),
        (Gegenbauer{1 / 4,T}, Jacobi{1 / 4 - 1 / 2,1 / 4 - 1 / 2,T}, f(1 / 4)),
        (Gegenbauer{3 / 4,T}, Jacobi{3 / 4 - 1 / 2,3 / 4 - 1 / 2,T}, f(3 / 4)),
        (Gegenbauer{1,T}, Jacobi{1 - 1 / 2,1 - 1 / 2,T}, f(1)),
    )
        for i in 1:5
            @test convert(Polynomial, Polynomials.basis(P, i)) ≈
                  fn(i) * convert(Polynomial, Polynomials.basis(Q, i))
        end
    end
end
