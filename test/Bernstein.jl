@testset "Construction" begin
    @testset for coeff in [
        Int64[1, 1, 1, 1],
        Float32[1, -4, 2],
        ComplexF64[1 - 1im, 2 + 3im],
        [3 // 4, -2 // 1, 1 // 1],
        # AbstractVector
        view([1,2,3], :),
        2:6,
    ]
        p = Bernstein(coeff)
        @test p.coeffs == coeff
        @test coeffs(p) == coeff
        @test degree(p) <= length(coeff) - 1
        @test Polynomials.indeterminate(p) == :x
        @test length(p) == length(coeff)
        @test size(p) == size(coeff)
        @test size(p, 1) == size(coeff, 1)
        @test typeof(p).parameters[end - 1] == eltype(coeff)
        @test eltype(p) == eltype(coeff)
    end

    # issue 84
    @test repr(basis(Bernstein{2}, 1)) == "Bernstein(1⋅β₂,₁(x))"
end

@testset "Other Construction" begin
    # Leading 0s are  trimmed
    p = Bernstein([1, 2, 0, 0])
    @test p.coeffs == [1, 2]

    # Different type
    p = Bernstein{3,Float64}(ones(Int32, 4))
    @test p.coeffs == ones(Float64, 4)

    p = Bernstein{1,Int}(30)
    @test p.coeffs == [30, 30]

    p = zero(Bernstein{1,Int})
    @test p.coeffs == Int[]

    p = zero(Bernstein{5,Int})
    @test length(p.coeffs) == 0

    p = one(Bernstein{1,Int})
    @test p.coeffs == ones(Int, 2)

    @test iszero(Bernstein{2,Int}(Int[]))

    p0 = Bernstein([0])
    @test iszero(p0)
    @test degree(p0) == -1
end

@testset "Roots" begin
    r = [1 / 4, 1 / 2, 3 / 4]
    c = fromroots(Bernstein, r)
    @test roots(c) ≈ r

    r = [-1im, 1 + 2im]
    c = fromroots(Bernstein, r)
    @test roots(c) ≈ r
end

@testset "Values" begin
    c1 = Bernstein([2.5, 2.0, 1.5]) # 1 + 2x + 3x^2
    x = [1 / 4, 1 / 2, 3 / 4]
    expected = [2.25, 2, 1.75]
    @test c1.(x) ≈ expected
end

@testset "Conversions" begin
    c1 = Bernstein([2.5, 2.0, 1.5])
    @test c1 ≈ Polynomial([2.5, -1.0])
    p = convert(Polynomial{Float64}, c1)
    @test p == Polynomial{Float64}([2.5, -1.0])
    @test convert(Bernstein{2,Float64}, p) == c1
end

@testset "Mapdomain" begin
    x = -30:20
    mx = mapdomain(Bernstein, x)
    @test extrema(mx) == (0, 1)

    x = 0.5:0.01:0.6
    mx = mapdomain(Bernstein, x)
    @test extrema(mx) == (0, 1)
end

@testset "Arithmetic" begin
    # multiplication
    c1 = Bernstein([1, 2, 3])
    c2 = Bernstein([3, 2, 1])
    @test convert(Polynomial, c1 * c2) ≈ convert(Polynomial, c1) * convert(Polynomial, c2)

    c1 = Bernstein([1, 2, 3])
    c2 = Bernstein([1, 2, 3, 4])
    p1, p2 = convert(Polynomial, c1), convert(Polynomial, c2)
    @test p1 * p2 ≈ convert(Polynomial, c1 * c2)

    c1 = Bernstein([1, 2, 3])
    c2 = Bernstein([3, 2, 1])
    d, r = divrem(c1, c2)
    @test d.coeffs ≈ [-1.0]
    @test r.coeffs ≈ [4.0]
    @test d * c2 + r ≈ c1

    c2 = Bernstein([0, 1, 2, 3])
    d, r = divrem(c2, c1)

    @test d.coeffs ≈ [1.5]
    @test r.coeffs ≈ [-1.5]
    z = d * c1 + r - c2
    @test coeffs(z) == zeros(length(z))

    # GCD
    c1 = Bernstein(Float64[1 // 1, 2, 3])
    c2 = Bernstein(Float64[3 // 1, 2, 1])
    @test gcd(c1, c2) ≈ Bernstein([4.0])

    for i in 0:5, j in 0:5
        # multiplication
        target = zeros(i + j + 1)
        target[end] += 1
        c1 = Bernstein(vcat(zeros(i), 1))
        c2 = Bernstein(vcat(zeros(j), 1))
        @test c1 * c2 ≈ Bernstein(target)

        # divrem
        target = c1 + c2
        quo, rem = divrem(target, c1)
        res = quo * c1 + rem
        @test res ≈ target
    end
end

@testset "addition" begin
    n = 4
    B = Bernstein{n}
    ps = [2,1,0,0,0]
    u,v = B(ps), sum(p*b for (p,b) ∈ zip(ps, basis.(B,0:n)))
    us,vs = getindex.(u,0:n), getindex.(v,0:n)
    us ≈ vs

    m = 5
    p,q = basis(B,3), basis(Bernstein{m},3)
    x = 0.2
    @test (p+q)(x) ≈ p(x) + q(x)
end

@testset "integrals derivatives" begin
    c1 = one(Bernstein{1,Int})
    @test integrate(c1) == Bernstein([0, 0.5, 1])
    for k in [-3, 0, 2]
        @test integrate(c1, k) == Bernstein(k .+ [0, 0.5, 1])
    end

    for i in 0:4
        scl = i + 1
        p = Polynomial(vcat(zeros(i), 1))
        pint = integrate(p, i)
        cheb = convert(Bernstein, p)
        cint = integrate(cheb, i)
        res = convert(Polynomial, cint)
        @test res ≈ pint
        @test derivative(cint) ≈ cheb
    end

    for i in 1:10
        p = Bernstein{5,Float64}(rand(1:5, 6))
        pp = p - (integrate ∘ derivative)(p)
        qq = p - (derivative ∘ integrate)(p)
        @test degree(truncate(pp, atol=1e-13)) <= 0
        @test degree(truncate(qq, atol=1e-13)) < 0
    end
end
