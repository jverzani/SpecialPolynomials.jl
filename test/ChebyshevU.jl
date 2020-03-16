@testset "Construction" begin

    # basic recurrence
    pNULL = ChebyshevU([0])
    p0 = ChebyshevU([1])
    p1 = ChebyshevU([0,1])
    p2 = ChebyshevU([0,0,1])
    p3 = ChebyshevU([0,0,0,1])
    p4 = ChebyshevU([0,0,0,0,1])
    p5 = ChebyshevU([0,0,0,0,0,1])
    p6 = Polynomials.basis(ChebyshevU{Int}, 6)

    # basic relation
    x = variable(p0)
    ps =  (p0, p1, p2, p3, p4, p5, p6)
    for n in 3:6
        @test ps[n+1] ≈ 2x*ps[n] - ps[n-1]
    end

    ## variable
    p = ChebyshevU([1,2,3])
    x =  variable(p)
    @test convert(Polynomial,p) == convert(Polynomial, -2 + 4*x + 12*x^2)



end

@testset "Conversion" begin

    for _ in 1:10
        ps = rand(1.0:8, 5)
        p = ChebyshevU(ps)
        @test -1 == degree(convert(ChebyshevU, convert(Polynomial, p)) - p)
    end

    for _ in 1:10
        ps = rand(1.0:8, 5)
        p = ChebyshevT(ps)
        @test -1 == degree(convert(ChebyshevT, convert(ChebyshevU, p)) - p)
    end

    for _ in 1:10
        ps = rand(1.0:8, 5)
        p = ChebyshevU(ps)
        @test -1 == degree(convert(ChebyshevU, convert(ChebyshevT, p)) - p)
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


@testset "Arithmetic" begin

    p0 = ChebyshevU([1])  # 1
    p1 = ChebyshevU([0,1]) # 2x
    p2 = ChebyshevU([0,0,1]) # 4x^2 - 1
    p3 = ChebyshevU([0,0,0,1]) # 8x^3 - 4x
    p4 = ChebyshevU([0,0,0,0,1]) # 16x^4 - 12x^2 + 1

    @test 1p0 + 2p1 + 3p2 + 4p3 + 5p4 == ChebyshevU([1,2,3,4,5])

    for p in (p0, p1, p2, p3, p4)
        for q in (p0, p1, p2, p3, p4)
            @test convert(Polynomial, p*q) == convert(Polynomial, p) * convert(Polynomial, q)
        end
    end

    for p in (p0, p1, p2, p3, p4)
        for q in (p0, p1, p2, p3, p4)
            @test convert(Polynomial, p+q) == convert(Polynomial, p) + convert(Polynomial, q)
        end
    end

end

@testset "Elementwise operations" begin

    xs = [1,2,3,4]
    p = ChebyshevU(xs)
    @test -p == ChebyshevU(-xs)
    @test p + 1 == ChebyshevU([1,0,0,0] + xs)
    @test 2p == ChebyshevU(2xs)


end



@testset "divrem" begin

    ps, qs = rand(1:6, 4), rand(1:6, 3)
    p, q = ChebyshevU(ps), ChebyshevU(qs)
    a, b =  divrem(p, q)
    z  = q  *  a +  b -  p
    @test degree(round(z, digits=13)) == -1

end


@testset "Derivatives and integrals" begin


    p = ChebyshevU([1.0, 2, 3, 4, 3, 2, 1])
    q = convert(Polynomial, p)
    @test derivative(p) == convert(ChebyshevU, derivative(q))


    p = ChebyshevU([1.0, 2, 3, 4, 3, 2, 1])
    q = convert(Polynomial, p)
    @test degree(round( integrate(p) - convert(ChebyshevU, integrate(q)), digits=13)) <= 0

    for _ in 1:10
        ps = rand(1:10, 5)
        p = ChebyshevU(ps)
        q = (derivative ∘ derivative)(p) - derivative(p, 2)  |> p->round(p, digits=13)
        @test degree(q) == -1

        q = (derivative ∘ integrate)(p) - p |> p->round(p, digits=13)
        @test degree(q) <= 0

        q = (integrate ∘ derivative)(p) - p |> p->round(p, digits=13)
        @test degree(q) <= 0
    end
end


@testset "roots" begin

    for n in 4:2:8
        p = ChebyshevU(vcat(zeros(n),1))
        @test cos.((n:-1:1)/(n+1) * pi) ≈ roots(p)
    end

end


@testset "fitting" begin

    f(x) = exp(-2pi*x) * sinpi(x)
    xs = range(0, 1, length=10)
    ys = f.(xs)
    q = fit(ChebyshevU, xs, ys, domain=(-1,1))
    @test maximum(abs.(q.(xs) - ys)) <= sqrt(eps())

    q = fit(ChebyshevU, xs, ys, 4)
    @test degree(q) <= 4

end
