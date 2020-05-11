

using  Test
function testit()
    Ps  = (Hermite,
           ChebyshevHermite,
           Laguerre{0}, Laguerre{1},
           Gegenbauer{1/2}, Gegenbauer{2},
           Bessel{1/2}, Bessel{1/4},
           Jacobi{1,1/2},Jacobi{1/2,2},
           Jacobi{1/2,1/2},
           Jacobi{-1/2,1/2},
           Jacobi{1/2,-1/2}, Jacobi{-1/2,-1/2}
           )

    x = variable()
    @show :arithmetic
    for P in Ps
        # p(x), -p,p+s,p*s,p/s, p+q, p*q
        as,bs = [1,2,3], [1,2,0]
        asc = [2,2,3]
        s,c,p,q = 1, P(1,:θ), P(as), P(bs)
        p(1/2)
        p(x)
        @test -p == P(-as)
        @test p+s == P(asc)
        @test p*s == P(as.*s)
        @test p/s == P(as ./ s)
        @test p+c == P(asc)
        @test p*c == P(as.*s)
        @test p+q == P(as+bs)
        #p*q
    end

    @show :integrate
    for P in Ps
        p = basis(P,3)
        out = integrate(p)(x) ≈ integrate(p(x))
        !out && @show  P
    end

    @show :conversion_within
    for  (P,Q) in  (
        (Hermite, Hermite),
        (Laguerre{1/2}, Laguerre{3/2}),
        (Gegenbauer{1/4}, Gegenbauer{3/4}),
        (Bessel{1/4},  Bessel{3/4}),
        (Jacobi{1/2,1/2}, Jacobi{3/4, 3/4}),
        (Jacobi{1/2,1/2}, Jacobi{-1/2, -1/2}),
                    )
        @show  P,Q
        p,q = basis(P,5), basis(Q,5)
        @test convert(P,q)(x) ≈ q(x)
        @test convert(Q,p)(x) ≈ p(x)
    end

    
    @show :conversion_Polynomial
    Q = Polynomial
    for P in Ps
        # P ∈ (#Hermite,
        #      #ChebyshevHermite,
        #      #Laguerre{0}, Laguerre{1},
        #      #Gegenbauer{1/2}, Gegenbauer{2},
        #      #Bessel{1}, Bessel{2},
        #      #Jacobi{1,1/2},
        #      #Jacobi{1/2,2},
        #      #Jacobi{1/2,1/2},
        #      #Jacobi{-1/2,1/2},
        #      #Jacobi{1/2,-1/2},
        #      #Jacobi{-1/2,-1/2},
        #      ) && continue
        @show P, Q
        p = P([1,2,3,4,5])
        q = Q([1,2,3,4,5])
        @test convert(Q, p) ≈ p(x) || hasnan(p(x))
        #@test convert(P, q)(x) ≈ q
    end

    @show :Jacobi
    for i in 1:5
        α, β = 2*rand(2) .- 1
        P = Jacobi{α,β}
        Q = Polynomial
        x = Q(:x)
        
        for n in 3:6
            p,q = basis(P,n), x^n
            @test convert(Q, _convert_connection_m(P, q)) ≈ q
            @test convert(P, _convert_connection_m(Q, p)) ≈ p            
            @test convert(Q, p) ≈ p(x)
            @test convert(P, q)(x) ≈ q
        end
    end

    @show :Jacobi_12
    for (α,β) = ((1/2, 1/2),
                 (1/2, -1/2),
                 #(-1/2, 1/2),
                 (-1/2, -1/2)
                 )
        @show α,β

        P = Jacobi{α,β}
        Q = Polynomial
        x = Q(:x)
        
        for n in 3:6
            p,q = basis(P,n),x^n
            @test convert(Q, convert(P, q)) ≈ q
            @test convert(P, convert(Q, p)) ≈ p            
            @test convert(Q, p) ≈ p(x) || hasnan(p(x))
            @test convert(P, q)(x) ≈ q
        end
    end
    
end






