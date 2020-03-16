@testset "Construction" begin

    pNULL = Legendre([0])
    p0 = Legendre([1])
    p1 = Legendre([0,1])
    p2 = Legendre([0,0,1])
    p3 = Legendre([0,0,0,1])
    p4 = Legendre([0,0,0,0,1])
    p5 = Legendre([0,0,0,0,0,1])
    p6 = Polynomials.basis(Legendre{Int}, 6)

    # basic relation
    x = variable(p0)
    ps =  (p1, p2, p3, p4, p5, p6)
    for n in 2:5
        @test (n+1)*ps[n+1] ≈ (2n+1)*x*ps[n] - n*ps[n-1]
    end

    ## variable
    p = Legendre([1,2,3])
    x =  variable(p)
    @test convert(Polynomial,p) == convert(Polynomial, -0.5 + 2.0x + 4.5x^2)



end

@testset "Conversion" begin

    for _ in 1:10
        ps = rand(1.0:8, 5)
        p = Legendre(ps)
        q = convert(Legendre, convert(Polynomial, p)) - p
        q = round(q, digits=13)
        @test -1 == degree(q)
    end


end

@testset "Evaluation" begin

    p0 = Legendre([1])  # 1
    p1 = Legendre([0,1]) # x
    p2 = Legendre([0,0,1]) #
    p3 = Legendre([0,0,0,1]) #
    p4 = Legendre([0,0,0,0,1]) #

    x = 1/2
    @test p0(x) ≈ 1
    @test p1(x) ≈ x
    @test p2(x) ≈ 1/2*(3x^2-1)
    @test p3(x) ≈ 1/2*(5x^3-3x)
    @test p4(x) ≈ 1/8*(35x^4 - 30x^2 + 3)
    @test Legendre([1,2,3,4])(x) ≈ 1p0(x) + 2p1(x) + 3p2(x) + 4p3(x)


end


@testset "Arithmetic" begin

    p0 = Legendre([1])  # 1
    p1 = Legendre([0,1]) #
    p2 = Legendre([0,0,1])
    p3 = Legendre([0,0,0,1]) #
    p4 = Legendre([0,0,0,0,1]) #

    @test 1p0 + 2p1 + 3p2 + 4p3 + 5p4 ≈ Legendre([1,2,3,4,5])

    for p in (p0, p1, p2, p3, p4)
        for q in (p0, p1, p2, p3, p4)
            @test convert(Polynomial, p*q) ≈ convert(Polynomial, p) * convert(Polynomial, q)
        end
    end

    for p in (p0, p1, p2, p3, p4)
        for q in (p0, p1, p2, p3, p4)
            @test convert(Polynomial, p+q) ≈ convert(Polynomial, p) + convert(Polynomial, q)
        end
    end

end

@testset "Elementwise operations" begin

    xs = [1,2,3,4]
    p = Legendre(xs)
    @test -p == Legendre(-xs)
    @test p + 1 == Legendre([1,0,0,0] + xs)
    @test 2p == Legendre(2xs)


end



@testset "divrem" begin

    ps, qs = rand(1:6, 4), rand(1:6, 3)
    p, q = Legendre(ps), Legendre(qs)
    a, b =  divrem(p, q)
    z  = q  *  a +  b -  p
    @test_broken degree(round(z, digits=13)) == -1

end


@testset "Derivatives and integrals" begin


    p = Legendre([1.0, 2, 3, 4, 3, 2, 1])
    q = convert(Polynomial, p)
    @test derivative(p) ≈ convert(Legendre, derivative(q))


    p = Legendre([1.0, 2, 3, 4, 3, 2, 1])
    q = convert(Polynomial, p)
    @test degree(round( integrate(p) - convert(Legendre, integrate(q)), digits=13)) <= 0

    for _ in 1:10
        ps = rand(1:10, 5)
        p = Legendre(ps)
        q = (derivative ∘ derivative)(p) - derivative(p, 2)  |> p->round(p, digits=13)
        @test degree(q) == -1

        q = (derivative ∘ integrate)(p) - p |> p->round(p, digits=13)
        @test degree(q) <= 0

        q = (integrate ∘ derivative)(p) - p |> p->round(p, digits=13)
        @test degree(q) <= 0
    end
end


@testset "roots" begin

    for n in 5:15
        p = Polynomials.basis(Legendre, n)
        @test all(isreal.(roots(p)))
    end

end


@testset "fitting" begin

    f(x) = exp(-2pi*x) * sinpi(x)
    xs = range(0, 1, length=10)
    ys = f.(xs)
    q = fit(Legendre, xs, ys, domain=(-1,1))
    @test maximum(abs.(q.(xs) - ys)) <= sqrt(eps())

    q = fit(Legendre, xs, ys, 4)
    @test degree(q) <= 4

end
