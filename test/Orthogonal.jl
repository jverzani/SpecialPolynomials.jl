## General test for orthogonal polynomials

T = Float64
Ps = (ChebyshevTT{T},
      ChebyshevU{T},
      Laguerre{T},
      Hermite{T},
      Jacobi{1/2, 1/2, T},
      Jacobi{0,0,T},
      Jacobi{0,1,T},
      Legendre{T},
      Gegenbauer{1/2,T},
      GeneralizedLaguerre{1/2, T}
      )


@testset "Construction" begin

    for P in Ps

        # basic recurrence
        pNULL = P([0])
        p0 = P([1])
        p1 = P([0,1])
        p2 = P([0,0,1])
        p3 = P([0,0,0,1])
        p4 = P([0,0,0,0,1])
        p5 = P([0,0,0,0,0,1])
        p6 = Polynomials.basis(P, 6)

        # basic relation
        x = variable(p0)
        ps =  (p0, p1, p2, p3, p4, p5, p6)
        for n in 3:6
            @test ps[n+1] ≈ (SP.An(P,n-1) * x + SP.Bn(P,n-1)) * ps[n] + SP.Cn(P,n-1)* ps[n-1]
        end

        for x in range(max(-1,first(domain(P))), min(1, last(domain(P))), length=10)
            @test P([1,2,3,4,5,6])(x) ≈ 1p0(x) + 2p1(x) + 3p2(x) + 4p3(x) + 5p4(x) + 6p5(x)
        end

    end


end

@testset "Conversion" begin

    for P in Ps


        for _ in 1:10
            ps = rand(1.0:8, 5)
            p = P(ps)
            q = convert(P, convert(Polynomial, p)) - p
            q = round(q, digits=10)
            @test -1 == degree(q)
        end

    end
end



@testset "Arithmetic" begin

    for P in Ps
        p0 = P([1])
        p1 = P([0,1])
        p2 = P([0,0,1])
        p3 = P([0,0,0,1])
        p4 = P([0,0,0,0,1])

        @test 1p0 + 2p1 + 3p2 + 4p3 + 5p4 ≈ P([1,2,3,4,5])

        for p in (p0, p1, p2, p3, p4)
            for q in (p0, p1, p2, p3, p4)
                @test convert(Polynomial, p+q) ≈ convert(Polynomial, p) + convert(Polynomial, q)
            end
        end


        for p in (p0, p1, p2, p3, p4)
            for q in (p0, p1, p2, p3, p4)
                @test convert(Polynomial, p*q) ≈ convert(Polynomial, p) * convert(Polynomial, q)
            end
        end

    end
end

@testset "Elementwise operations" begin

    for P in Ps
        xs = [1,2,3,4]
        p = P(xs)
        @test -p == P(-xs)
        @test p + 1 == P([1,0,0,0] + xs)
        @test 2p == P(2xs)
    end

end



@testset "divrem" begin

    for P in Ps
        ps, qs = rand(1:6, 4), rand(1:6, 3)
        p, q = P(ps), P(qs)
        a, b =  divrem(p, q)
        z  = q  *  a +  b -  p
        if P <: Legendre
            @test_broken degree(round(z, digits=13)) == -1
        else
            @test degree(round(z, digits=13)) == -1
        end
    end

end


@testset "Derivatives and integrals" begin
    _truncate(x) = truncate(x, atol=sqrt(eps(T)))
    for P in Ps
        p = P([1.0, 2, 3, 4, 3, 2, 1])
        q = convert(Polynomial, p)
        @test derivative(p) ≈ convert(P, derivative(q))


        p = P([1.0, 2, 3, 4, 3, 2, 1])
        q = convert(Polynomial, p)
        @test degree(_truncate( integrate(p) - convert(P, integrate(q)))) <= 0

        for _ in 1:10
            ps = rand(1:10, 5)
            p = P(ps)
            q = (derivative ∘ derivative)(p) - derivative(p, 2)
            @test degree(_truncate(q)) == -1

            q = (derivative ∘ integrate)(p) - p
            @test degree(_truncate(q)) <= 0

            q = (integrate ∘ derivative)(p) - p
            @test degree(_truncate(q)) <= 0
        end
    end

end

@testset "roots" begin

    for P in Ps
        for n in 6:2:12
            rts = roots(Polynomials.basis(P, n))
            evals = eigvals(SpecialPolynomials.jacobi_matrix(P, n))
            @test all(rts .≈ evals)
        end
    end

end


@testset "fitting" begin

    f(x) = exp(-2pi*x) * sinpi(x)
    xs = range(0, 1, length=10)
    ys = f.(xs)

    for P in Ps
        if !(P <: Laguerre || P <: GeneralizedLaguerre || P <: Hermite)
            q = fit(P, xs, ys, domain=domain(P))
            @test maximum(abs.(q.(xs) - ys)) <= sqrt(eps())
        end

        q = fit(P, xs, ys, 4)
        @test degree(q) <= 4
    end
end

@testset "quadrature" begin

    f(x) = x^7
    n = 4

    for P in Ps
        !all(isfinite.(extrema(P))) && continue
        q = sum(f(tau)*w for (tau, w)  in  zip(SP.gauss_nodes_weights(P,n)...))
        p = SP._quadgk(x -> f(x) * SP.weight_function(P)(x), extrema(P)...)
        @test abs(p - q)  <= sqrt(eps(T))
     end

end
