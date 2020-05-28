## General test for orthogonal polynomials

T = Float64
Ps = (
    Chebyshev,
    ChebyshevU,
    Laguerre{0},
    Laguerre{1/2},
    Hermite,
    ChebyshevHermite,
    Jacobi{1/2, 1/2},
    Jacobi{1/4, 3/4},
    Jacobi{-1/2, 1/2},
    Jacobi{1/2, -1/2},      
    Legendre,
    Gegenbauer{1/2},
    Bessel{3/2}, # Bessel{1} is an issue
    Bessel{1/2}
      )

@testset "Construction" begin
    
    for P in Ps
        x = variable(P)
        # basic recurrence
        for n in 1:6
            @test basis(P,n+1) ≈ (SP.An(P,n) * x + SP.Bn(P,n)) * basis(P,n) - SP.Cn(P,n)*basis(P,n-1)
        end

        # leading term and ratios
        x = variable(Polynomial)
        for n = 0:5
            @test basis(P,n)(x)[end] ≈ SP.kn(P,n)
        end
        for n = 0:5
            @test SP.k1k0(P,n) ≈ SP.kn(P,n+1)/SP.kn(P,n)
        end
        for n = 1:5
            @test SP.k1k_1(P,n) ≈ SP.kn(P,n+1)/SP.kn(P,n-1)
        end
        
        # basis
        x = variable(P)
        p0, p1, p2, p3, p4, p5, p6 = ps =  basis.(P,0:6)
        p = P([1,2,3,4,5,6])
        for x in range(max(-1,first(domain(P))), stop=min(1, last(domain(P))), length=10)
            @test p(x) ≈ 1p0(x) + 2p1(x) + 3p2(x) + 4p3(x) + 5p4(x) + 6p5(x)
        end


        ## Different type
        T = Float64
        @test eltype(P(ones(Int32, 4))) ==  Int32
        @test eltype(P{T}(ones(Int32, 4))) ==  T
        @test degree(P(1)-1) ==  -1
        @test degree(zero(P)) == -1
        @test degree(one(P)) == 0
        @test degree(variable(P)) == 1

        # variable and P() constructor for `x` in basis
        @test degree(variable(P)) == 1
        @test variable(P)(1) ≈ 1
        @test degree(P()) == 1
        @test P()(1) ≈ 1
        @test variable(P, :y) == P(:y)

        # test  ignore constant's symbol
        @test zero(P, :x) == zero(P, :y)
        @test one(P, :x) == one(P, :y)
        @test !(variable(P, :x) == variable(P,:y))

        @test !(zero(P, :x) === zero(P, :y))
        @test !(one(P, :x) === one(P, :y))
        @test !(variable(P, :x) === variable(P,:y))

        @test zero(P, :x) ≈ zero(P, :y)
        @test one(P, :x) ≈ one(P, :y)
        @test (variable(P, :x) ≈ variable(P, :x))
        @test_throws ErrorException variable(P, :x) ≈ variable(P, :y)

    end

end

@testset "Structural equations" begin
    
    for P in filter(P -> !(P <: Hermite), Ps)
        for i in 2:5

            ps = basis.(P, i+1:-1:i-1)
            dps = derivative.(ps)
            pᵢ,dpᵢ = ps[2], dps[2]
            
            a,b,c,d,e = SP.abcde(P)
            x = variable(P)
            σ = a*x^2 + b*x + c

            
            # x⋅p_n   = [an, bn, cn]    ⋅ [p_{n+1}, p_n, p_{n-1}]     #  Eqn (7)
            as= SP.an(P,i), SP.bn(P,i), SP.cn(P,i)
            @test x*pᵢ ≈ sum(a*p for (a,p) in zip(as, ps))

            
            # σ⋅p'_n  = [αn, βn, γn]    ⋅  [p_{n+1}, p_n, p_{n-1}]    # Eqn (9), n ≥ 1
            αs = SP.αn(P,i), SP.βn(P,i), SP.γn(P,i)
            @test σ*dpᵢ ≈ sum(a*p for (a,p) in zip(αs, ps))
                
            # p_n    = [ân, b̂n, ĉn]    ⋅  [p'_{n+1}, p'_n, p'_{n-1}] # Eqn (19)
            âs = SP.ân(P,i), SP.b̂n(P,i), SP.ĉn(P,i)
            @test pᵢ ≈ sum(a*dp for (a,dp) in zip(âs, dps))
            
            # x⋅p'_n  = [αᴵn, βᴵn, γᴵn] ⋅  [p'_{n+1}, p'_n, p'_{n-1}] # Eqn  (14) with  α^*, β^*,  γ^*
            @test (a=a,b=b,c=c,d=d+2a,e=e+b)  == SP.abcdeᴵ(P)
            αᴵs = SP.αᴵn(P,i), SP.βᴵn(P,i), SP.γᴵn(P,i) 
            @test x * dpᵢ ≈ sum(a*dp for (a,dp) in zip(αᴵs, dps))
        end
    end
end
        

@testset "Conversion" begin

    for P in Ps
        ## convert to polynomial and back
        Q = Polynomial
        for _ in 1:10
            ps = rand(1.0:8, 5)
            p = P(ps)
            @test convert(P, convert(Q, p)) ≈ p
        end

    end
    
    # through evaluation
    for P in  Ps
        for Q in Ps
            p = P([1,2,3])
            x = variable(Q)
            @test p(x) ≈ convert(Q,  p)
        end
    end

  

    x = variable(Polynomial)
    # Connections: groupings have same sigma = a⋅x² + b⋅x + c
    for Ps in ((Chebyshev, ChebyshevU, Gegenbauer{1/2}, Gegenbauer{1},
                Jacobi{1/4, 1/4}, Jacobi{3/4, 3/4}),
               (Laguerre{0},  Laguerre{1/2}, Laguerre{1}),
               (Bessel{1/4}, Bessel{1/2}, Bessel{3/4})#, Bessel{2})
               )

        
        for  P in Ps
            for Q in Ps
                for n in 2:5
                    q = basis(Q, n)
                    # various ways to convert to P
                    
                    @test q(variable(P))(variable()) ≈ q(variable())  # evaluation
                    @test SP._convert_ccop(P, q)(variable()) ≈ q(variable())  # structural  equations
                    @test convert(P, convert(Polynomial, q))(variable()) ≈ q(variable()) # through Standard  basis
                    @test convert(P, q)(x) ≈ q(x) 

                end
            end
        end
    end
    
end


@testset "Arithmetic" begin

    x = variable(Polynomial{Float64})

    for P in Ps
        pNULL = zero(P)
        p0,p1,p2,p3,p4 = ps = basis.(P,0:4)

        as =  [1,2,3,4,5]
        @test sum(a*pᵢ for (a,pᵢ) in zip(as, ps)) ≈ P(as)

        x = variable(Polynomial)
        for p in ps
            for q in ps
                for op in (+, *)
                    @test op(p,q)(x) ≈ op(p(x),q(x))
                end
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
        @test p ≈ q  *  a +  b 
    end

end


@testset "Derivatives and integrals" begin

    for P in Ps
        Q = Polynomial
        x = variable(Q)

        for i in 1:6
            p = basis(P,i)
            @test derivative(p)(x) ≈ derivative(p(x))
        end


        p = P([1.0, 2, 3, 4, 3, 2, 1])
        q = convert(Polynomial, p)
        @test maximum(abs, [integrate(p, a, a+1/2) - integrate(q, a, a+1/2) for a in range(0, stop=1/2, length=10)]) <= sqrt(eps())

        for _ in 1:10
            ps = rand(1:10, 5)
            p = P(ps)
            @test  (derivative ∘ derivative)(p) ≈ derivative(p, 2)
            
            @test (derivative ∘ integrate)(p) ≈ p

            q = (integrate ∘ derivative)(p) - p
            q = chop(q, atol=1e-8)
            @test degree(q) <= 0
        end
    end

end

# @testset "roots" begin

#     for P in Ps
#         for n in 6:2:12
#             rts = roots(basis(P, n))
#             evals = eigvals(SpecialPolynomials.jacobi_matrix(P, n))
#             @test maximum(abs, sort(rts) - sort(evals)) <= 1e-6
#         end
#     end
    
# end


@testset "fitting" begin

    f(x) = exp(-x) * cospi(x)

    # Interpolating
    for P in Ps
        P <: SP.AbstractOrthogonalPolynomial || continue
        dom = domain(P)
        (isinf(first(dom)) || isinf(last(dom))) && continue
        q = fit(P, f, 10)
        @test maximum(abs, q(x) -  f(x) for x in range(0, stop=1/2, length=10)) <= 1e-4
    end

    # least squares
    for P in (Chebyshev,)
        q = fit(Val(:lsq), P, f, 10)
        @test maximum(abs, q(x) -  f(x) for x in range(0, stop=1/2, length=10)) <= 1e-4
    end

    # series
    for  P in  (Chebyshev,)
        q = fit(Val(:series), P, f)
        @test maximum(abs, q(x) -  f(x) for x in range(0, stop=1/2, length=10)) <= sqrt(eps())
    end
        
end

@testset "quadrature" begin

    f(x) = x^7
    n = 4

    for P in Ps
        P <: SP.AbstractOrthogonalPolynomial || continue
        !all(isfinite.(extrema(P))) && continue
        q = sum(f(tau)*w for (tau, w)  in  zip(SP.gauss_nodes_weights(P,n)...))
        p = SP.innerproduct(P, f, one)
        @test abs(p - q)  <= 10sqrt(eps(T))
     end

end


