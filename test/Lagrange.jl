@testset "Construction" begin

    # basic constructor
    xs = rand(5)
    ys = rand(5)
    p = Lagrange(xs, ys)
    @test all(p.(xs) .≈ ys)

    # test weights
    @test all(p.ws .≈ SP.lagrange_barycentric_weights(xs))
    zs = rand(5)
    q = Lagrange(xs, p.ws, zs)
    @test all(q.(xs) .≈ zs)

    # fit is alternate
    q = fit(Lagrange, xs, ys)
    us = rand(5)
    @test all(p.(us) .≈ q.(us))

    # test one, zero, variable
    x = variable()
    q = fit(Lagrange, [1, 2, 3], [2, 3, 1])
    @test iszero(zero(q))
    @test convert(Polynomial, one(q)) == one(x)
    @test convert(Polynomial, variable(q)) == x
end

@testset "Conversion" begin
    xs = [1, 2, 3]
    ys = xs .^ 2
    p = fit(Lagrange, xs, ys)
    q = convert(Polynomial, p)
    @test coeffs(q) == [0, 0, 1]

    pp = fit(Lagrange, xs, q)
    @test coeffs(pp) == coeffs(p)
end

@testset "Evaluation" begin
    xs = rand(5)
    ys = sin.(xs)
    p = Lagrange(xs, ys)
    @test all(p.(xs) == ys)   # shortcuts  for xs
    us = rand(5)
    @test all(abs.(p.(us) - sin.(us)) .<= 1e-1)

    xs, ws = SP.lagrange_barycentric_nodes_weights(Chebyshev, 64)
    ys = sin.(xs)
    p = Lagrange(xs, ws, ys)
    us = 0.1 .+ 0.8 * rand(5) # avoid edges...
    @test all(abs.(sin.(us) - p.(us)) .<= 1e-3)
end

@testset "Arithmetic" begin
    p = fit(Lagrange, [1, 2, 3], x -> x^2)
    q = fit(Lagrange, [2, 3, 4, 5], x -> x^3)

    us = rand(5)

    @test maximum(abs, p.(us) + q.(us) - (p + q).(us)) <= 1e-5
    @test maximum(abs, p.(us) .* q.(us) - (p * q).(us)) <= 1e-2
end

@testset "Elementwise operations" begin
    p = fit(Lagrange, [1, 2, 3], x -> x^2)
    p + 1
    -p
    2p
end

@testset "roots" begin
    q = fit(Lagrange, [2, 3, 4, 5], x -> (x - 1) * (x - 1 / 2) * (x + 1 / 2))
    @test all(roots(q) .≈ [-1 / 2, 1 / 2, 1])
end
