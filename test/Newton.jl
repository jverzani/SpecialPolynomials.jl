@testset "Construction" begin

    # basic constructor
    xs = [1, 2, 3]
    f(x) = x^2
    q = fit(Newton, xs, f)
    @test all(f.(xs) .≈ q.(xs))

    # test tableau
    @test all(q.tableau .≈ [1 3 1; 0 4 5; 0 0 9])

    # test one, zero, variable
    x = variable()
    @test one(q)(x) == one(x)
    @test iszero(zero(q))
    @test variable(q)(x) == x
end

@testset "Conversion" begin
    xs = [1, 2, 3]
    f(x) = x^2
    p = fit(Newton, xs, f)
    q = convert(Polynomial, p)
    @test coeffs(q) == [0, 0, 1]
end

@testset "Evaluation" begin
    xs = sort(rand(5))
    ys = sin.(xs)
    p = fit(Newton, xs, sin)
    @test p.(xs) ≈ sin.(xs)
    us = rand(5)
    @test maximum(abs.([p(x) - sin(x) for x in us])) <= 1e-1
end

@testset "Arithmetic" begin
    p = fit(Newton, [1, 2, 3], x -> x^2)
    q = fit(Newton, [2, 3, 4, 5], x -> x^3)

    us = rand(5)

    @test all(p.(us) + q.(us) .≈ (p + q).(us))
    @test p.(us) .* q.(us) ≈ (p * q).(us)
end

@testset "Elementwise operations" begin
    p = fit(Newton, [1, 2, 3, 4], x -> x^2)
    p + 1
    -p
    2p
end

@testset "roots" begin
    q = fit(Newton, [2, 3, 4, 5], x -> (x - 1) * (x - 1 / 2) * (x + 1 / 2))
    @test all(roots(q) .≈ [-1 / 2, 1 / 2, 1])
end
