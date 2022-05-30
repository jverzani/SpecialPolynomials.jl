## Note:
## This is the `ChebyshevT.jl` test suite from the `Polynomials.jl` package by Miles Lucas

@testset "Construction" begin
    @testset for coeff in [
        Int64[1, 1, 1, 1],
        Float32[1, -4, 2],
        ComplexF64[1 - 1im, 2 + 3im],
        [3 // 4, -2 // 1, 1 // 1],
        # AbstractVector
        view([1.0, 2.0], :),
        3:6
    ]
        p = @inferred Chebyshev(coeff)
        @test p.coeffs == coeff
        @test coeffs(p) == coeff
        @test degree(p) == length(coeff) - 1
        @test Polynomials.indeterminate(p) == :x
        @test length(p) == length(coeff)
        @test size(p) == size(coeff)
        @test size(p, 1) == size(coeff, 1)
        @test typeof(p).parameters[1] == eltype(coeff)
        @test eltype(p) == eltype(coeff)
    end
end

@testset "Other Construction" begin
    # Leading 0s
    p = Chebyshev([1, 2, 0, 0])
    @test p.coeffs == [1, 2]
    @test length(p) == 2

    # Different type
    p = Chebyshev{Float64}(ones(Int32, 4))
    @test p.coeffs == ones(Float64, 4)

    p = Chebyshev(30)
    @test p.coeffs == [30]

    p = zero(Chebyshev{Int})
    @test p.coeffs == Int[]

    p = one(Chebyshev{Int})
    @test p.coeffs == [1]

    pNULL = Chebyshev(Int[])
    @test iszero(pNULL)
    @test degree(pNULL) == -1

    p0 = Chebyshev([0])
    @test iszero(p0)
    @test degree(p0) == -1
end

@testset "Roots $i" for i in 1:5
    roots = cos.(range(-π, 0, length=2i + 1)[2:2:end])
    target = Chebyshev(vcat(zeros(i), 1))
    res = fromroots(ChebyshevT, roots) .* 2^(i - 1)
    @test res == target
end

@testset "Roots" begin
    r = [-1, 1, 0]
    c = fromroots(ChebyshevT, r)
    @test c == Chebyshev([0, -0.25, 0, 0.25])
    @test roots(c) ≈ r

    r = [1im, -1im]
    c = fromroots(ChebyshevT, r)
    @test c ≈ Chebyshev([1.5 + 0im, 0 + 0im, 0.5 + 0im])
    @test roots(c) ≈ r
end

@testset "Values" begin
    c1 = Chebyshev([2.5, 2.0, 1.5]) # 1 + 2x + 3x^2
    x = -1:0.5:1
    expected = [2.0, 0.75, 1.0, 2.75, 6.0]
    @test c1.(x) == expected
end

@testset "Conversions" begin
    c1 = Chebyshev([2.5, 2.0, 1.5])
    @test c1 ≈ Polynomial([1, 2, 3])
    p = convert(Polynomial{Float64}, c1)
    @test p == Polynomial{Float64}([1, 2, 3])
    @test convert(ChebyshevT, p) == c1
end

@testset "Companion" begin
    c_null = Chebyshev(Int[])
    c_1 = Chebyshev([1])
    @test_throws ArgumentError companion(c_null)
    @test_throws ArgumentError companion(c_1)
    for i in 1:5
        coef = vcat(zeros(i), 1)
        c = Chebyshev(coef)
        @test size(companion(c)) == (i, i)
    end
    c = Chebyshev([1, 2])
    @test companion(c)[1, 1] == -0.5
end

@testset "Vander" begin
    x = 0:0.1:1
    v = vander(ChebyshevT, x, 5)
    @test size(v) == (length(x), 6)
    @inbounds for i in 1:6
        coef = vcat(zeros(i - 1), 1)
        c = Chebyshev(coef)
        @test v[:, i] ≈ c.(x)
    end
end

@testset "Arithmetic $i, $j" for i in 0:5, j in 0:5
    # multiplication
    target = zeros(i + j + 1)
    target[end] += 0.5
    target[abs(i - j) + 1] += 0.5
    c1 = Chebyshev(vcat(zeros(i), 1))
    c2 = Chebyshev(vcat(zeros(j), 1))
    @test c1 * c2 ≈ Chebyshev(target)

    # divrem
    target = c1 + c2
    quo, rem = divrem(target, c1)
    res = quo * c1 + rem
    @test res ≈ target
end

@testset "Mapdomain" begin
    x = -30:20
    mx = mapdomain(ChebyshevT, x)
    @test extrema(mx) == (-1, 1)

    x = 0.5:0.01:0.6
    mx = mapdomain(ChebyshevT, x)
    @test extrema(mx) == (-1, 1)
end

@testset "Arithmetic" begin
    # multiplication
    c1 = Chebyshev([1, 2, 3])
    c2 = Chebyshev([3, 2, 1])
    @test c1 * c2 == Chebyshev([6.5, 12, 12, 4, 1.5])

    c1 = Chebyshev([1, 2, 3])
    c2 = Chebyshev([3, 2, 1])
    d, r = divrem(c1, c2)
    @test d.coeffs ≈ [3]
    @test r.coeffs ≈ [-8, -4]

    c2 = Chebyshev([0, 1, 2, 3])
    d, r = divrem(c2, c1)

    @test d.coeffs ≈ [0, 2]
    @test r.coeffs ≈ [-2, -4]

    # GCD
    c1 = Chebyshev([1, 2, 3])
    c2 = Chebyshev([3, 2, 1])
    @test gcd(c1, c2) ≈ Chebyshev(6)
end

@testset "integrals derivatives" begin
    c1 = one(Chebyshev{Int})
    @test integrate(c1) == Chebyshev([0, 1])
    for k in [-3, 0, 2]
        @test integrate(c1, k) == Chebyshev([k, 1])
    end

    for i in 0:4
        scl = i + 1
        p = Polynomial(vcat(zeros(i), 1))
        target = Polynomial(vcat(i, zeros(i), 1 / scl))
        cheb = convert(ChebyshevT, p)
        cint = integrate(cheb, i)
        res = convert(Polynomial, cint)
        @test res ≈ target
        @test derivative(cint) == cheb
    end

    ## issue #56
    coeffs  = [3 // 4, -2 // 1, 1 // 1]
    p = Chebyshev(coeffs)
    @test derivative(p) ≈ convert(Chebyshev, derivative(convert(Polynomial, p)))
    @test derivative(p, 2) ≈ convert(Chebyshev, derivative(convert(Polynomial, p), 2))

end

@testset "z-series" for i in 0:5
    # c to z
    input = vcat(2, ones(i))
    target = vcat(fill(0.5, i), 2, fill(0.5, i))
    zs = Polynomials._c_to_z(input)
    @test zs == target
    c = Polynomials._z_to_c(zs)
    @test c == input

    # div
    z1 = [0.5, 0.0, 0.5]
    z2 = [0.5, 0.0, 0.5]
    q, r = Polynomials._z_division(z1, z2)
    @test q == [1]
    @test r == [0]
end

@testset "Fit" begin
    f(x) = exp(x)
    p = fit(Chebyshev{Float64}, f, 16)
    @test maximum(abs.(p(x) - f(x) for x in range(0, 1, length=100))) <= 100 * eps()

    f(x) = sin(6x) - sin(60 * exp(x)) # chebfun example
    p = fit(Chebyshev{Float64}, f, 150)
    @test maximum(abs.(p(x) - f(x) for x in range(0, 1, length=100))) <= 1000 * eps()
end
