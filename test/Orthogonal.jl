## General test for orthogonal polynomials

T = Float64
Ps = (Chebyshev{T},
      ChebyshevU{T},
      Laguerre{T},
      Hermite{T},
      ChebyshevHermite{T},
      Jacobi{1/2, 1/2, T},
      Jacobi{0,0,T},
      Jacobi{0,1,T},
      Legendre{T},
      Gegenbauer{1/2,T},
      GeneralizedLaguerre{1/2, T}
#      ,DiscreteChebyshev{12,T}
#      ,Krawtchouk{12, 5/2, T}
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

        p = P([1,2,3,4,5,6])
        for x in range(max(-1,first(domain(P))), stop=min(1, last(domain(P))), length=10)
            @test p(x) ≈ 1p0(x) + 2p1(x) + 3p2(x) + 4p3(x) + 5p4(x) + 6p5(x)
        end


        ## Different type
        @test eltype(coeffs(P(ones(Int32, 4)))) ==  T
        @test degree(P(1) -1) ==  -1
        @test degree(zero(P)) == -1
        @test degree(one(P)) == 0
        @test degree(variable(P)) == 1

        # variable and P() constructor for `x` in basis
        @test degree(variable(P)) == 1
        @test variable(P)(1) == 1
        @test degree(P()) == 1
        @test P()(1) == 1
        @test variable(P, :y) == P(:y)


    end




    
end

@testset "Conversion" begin

    for P in Ps


        for _ in 1:10
            ps = rand(1.0:8, 5)
            p = P(ps)
            q = convert(P, convert(Polynomial, p)) - p
            truncate!(q, atol=sqrt(eps(T)))
            @test -1 == degree(q)
        end

    end
end


@testset "Arithmetic" begin

    x = variable(Polynomial{Float64})

    for P in Ps
        pNULL = zero(P)
        p0 = P([1])
        p1 = P([0,1])
        p2 = P([0,0,1])
        p3 = P([0,0,0,1])
        p4 = P([0,0,0,0,1])

        @test 1p0 + 2p1 + 3p2 + 4p3 + 5p4 ≈ P([1,2,3,4,5])

        for p in (pNULL, p0, p1, p2, p3, p4)
            for q in (pNULL, p0, p1, p2, p3, p4)
                @test convert(Polynomial, p+q) ≈ convert(Polynomial, p) + convert(Polynomial, q)
                @test convert(Polynomial, p*q) ≈ convert(Polynomial, p) * convert(Polynomial, q)                
            end
        end

        for _ in 1:10
            p,q = P(rand(5)), P(rand(5))
            for op in (+,*)
                @test (op(p,q))(x) ≈ op(p(x), q(x))
            end
        end
        
    end
end

@testset "Elementwise operations" begin

    for P in Ps
        xs = [1,2,3,4]
        p = P(xs)
        @test -p == P(-xs)
        ys = copy(xs)
        ys[1] += 1
        @test p + 1 == P(ys)
        @test 2p == P(2xs)
    end

end

# too slow
# @testset "Orthogonality" begin
#     n = 5
#     for P in Ps
#         @show P
#         @time for i in 2:n
#             for j in i+1:n
#                 val = SP.innerproduct(P, Polynomials.basis(P, i), Polynomials.basis(P,j))
#                 @test abs(val)  <= 1e-4
#             end
#         end
#     end
# end

@testset "divrem" begin

    for P in Ps
        ps, qs = rand(1:6, 4),rand(1:6, 3)
        p, q = P(ps), P(qs)
        a, b =  divrem(p, q)
        z  = q  *  a +  b -  p
        @test degree(truncate!(z, atol=sqrt(eps()))) == -1
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
        @test maximum(abs, [integrate(p, a, a+1/2) - integrate(q, a, a+1/2) for a in range(0, stop=1/2, length=10)]) <= sqrt(eps())

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
            @test maximum(abs, sort(rts) - sort(evals)) <= 1e-6
        end
    end
    
end


@testset "fitting" begin

    f(x) = exp(-2pi*x) * cospi(x)

    for P in Ps
        P <: SP.OrthogonalPolynomial || continue
        dom = domain(P)
        (isinf(first(dom)) || isinf(last(dom))) && continue
        q = fit(P, f, 10)
        @test maximum(abs, q(x) -  f(x) for x in range(0, stop=1/2, length=10)) <= 1e-1
    end

end

@testset "quadrature" begin

    f(x) = x^7
    n = 4

    for P in Ps
        P <: SP.OrthogonalPolynomial || continue
        !all(isfinite.(extrema(P))) && continue
        q = sum(f(tau)*w for (tau, w)  in  zip(SP.gauss_nodes_weights(P,n)...))
        p = SP.innerproduct(P, f, one)
        @test abs(p - q)  <= sqrt(eps(T))
     end

end
