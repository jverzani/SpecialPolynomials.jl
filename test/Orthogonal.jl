## General test for orthogonal polynomials
import SpecialPolynomials: domain
using LinearAlgebra

T = Float64
Ps = (
    Chebyshev,
    OrthonormalChebyshev,
    ChebyshevU,
    Laguerre{0},
    Laguerre{1 / 2},
    OrthonormalLaguerre{1 / 2},
    Hermite,
    ChebyshevHermite,
    Jacobi{1 / 2,1 / 2},
    Jacobi{1 / 4,3 / 4},
    Jacobi{-1 / 2,1 / 2},
    Jacobi{1 / 2,-1 / 2},
    OrthonormalJacobi{1 / 2,1 / 2},
    Legendre,
    OrthonormalLegendre,
    Gegenbauer{1 / 2},
    OrthonormalGegenbauer{1 / 2},
    Bessel{3 / 2}, # Bessel{1} is an issue
    Bessel{1 / 2},
)
DPs = (Charlier{1 / 2}, Meixner{1 / 2,1 / 2}, Krawchouk{1 / 2,10}, Hahn{1 / 4,1 / 2,10})

@testset "Construction" begin
    @testset for P in union(Ps, DPs)
        x = variable(P)
        # basic recurrence
        @testset for n in 1:6
            @test basis(P, n + 1) ≈
                  (SP.An(P, n) * x + SP.Bn(P, n)) * basis(P, n) -
                  SP.Cn(P, n) * basis(P, n - 1)
        end

        # leading term and ratios
        x = variable(Polynomial)
        @testset for n in 0:5
            @test basis(P, n)(x)[end] ≈ SP.kn(P, n)
        end
        @testset for n in 0:5
            @test SP.k1k0(P, n) ≈ SP.kn(P, n + 1) / SP.kn(P, n)
        end
        @testset for n in 1:5
            @test SP.k1k_1(P, n) ≈ SP.kn(P, n + 1) / SP.kn(P, n - 1)
        end

        # basis
        x = variable(P)
        p0, p1, p2, p3, p4, p5, p6 = ps = basis.(P, 0:6)
        p = P([1, 2, 3, 4, 5, 6])
        @testset for x in range(
            max(-1, first(domain(P))),
            stop=min(1, last(domain(P))),
            length=10,
        )
            @test p(x) ≈ 1p0(x) + 2p1(x) + 3p2(x) + 4p3(x) + 5p4(x) + 6p5(x)
        end

        ## Different type
        T = Float64
        @test eltype(P(ones(Int32, 4))) == Int32
        @test eltype(P{T}(ones(Int32, 4))) == T
        @test degree(P(1) - 1) == -1
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
        @test !(variable(P, :x) == variable(P, :y))

        @test !(zero(P, :x) === zero(P, :y))
        @test !(one(P, :x) === one(P, :y))
        @test !(variable(P, :x) === variable(P, :y))

        @test zero(P, :x) ≈ zero(P, :y)
        @test one(P, :x) ≈ one(P, :y)
        @test (variable(P, :x) ≈ variable(P, :x))
        @test_throws ArgumentError variable(P, :x) ≈ variable(P, :y)
    end

    @testset for P in Ps
        @inferred P([1, 2, 3]) # not @inferred P([1,2,3], :x)
        @inferred P{Int}([1, 2, 3]) # not @inferred P{Int}([1,2,3], :x)
        @inferred P{Int,:x}([1, 2, 3])
    end
end

@testset "Structural equations" begin
    @testset for P in Ps
        P <: Hermite && continue #  Hermite  has issues
        @testset for i in 2:5
            ps = basis.(P, (i + 1):-1:(i - 1))
            dps = derivative.(ps)
            pᵢ, dpᵢ = ps[2], dps[2]

            a, b, c, d, e = SP.abcde(P)
            x = variable(P)
            σ = a * x^2 + b * x + c

            # x⋅p_n   = [an, bn, cn]    ⋅ [p_{n+1}, p_n, p_{n-1}]     #  Eqn (7)
            as = SP.an(P, i), SP.bn(P, i), SP.cn(P, i)
            @test x * pᵢ ≈ sum(a * p for (a, p) in zip(as, ps))

            # σ⋅p'_n  = [αn, βn, γn]    ⋅  [p_{n+1}, p_n, p_{n-1}]    # Eqn (9), n ≥ 1
            αs = SP.αn(P, i), SP.βn(P, i), SP.γn(P, i)
            @test σ * dpᵢ ≈ sum(a * p for (a, p) in zip(αs, ps))

            # p_n    = [ân, b̂n, ĉn]    ⋅  [p'_{n+1}, p'_n, p'_{n-1}] # Eqn (19)
            âs = SP.ân(P, i), SP.b̂n(P, i), SP.ĉn(P, i)
            @test pᵢ ≈ sum(a * dp for (a, dp) in zip(âs, dps))

            # x⋅p'_n  = [αᴵn, βᴵn, γᴵn] ⋅  [p'_{n+1}, p'_n, p'_{n-1}] # Eqn  (14) with  α^*, β^*,  γ^*
            @test (a=a, b=b, c=c, d=d + 2a, e=e + b) == SP.abcdeᴵ(P)
            αᴵs = SP.αᴵn(P, i), SP.βᴵn(P, i), SP.γᴵn(P, i)
            @test x * dpᵢ ≈ sum(a * dp for (a, dp) in zip(αᴵs, dps))
        end
    end

    @testset for P in DPs
        @testset for i in 1:5
            ps = basis.(P, (i + 1):-1:(i - 1))
            dps = SP.:∇ₓ.(ps)
            Δps = SP.:Δₓ.(ps)
            pᵢ, dpᵢ, Δpᵢ = ps[2], dps[2], Δps[2]

            a, b, c, d, e = SP.abcde(P)
            x = variable(P)
            σ = a * x^2 + b * x + c

            # x⋅p_n   = [an, bn, cn]    ⋅ [p_{n+1}, p_n, p_{n-1}]     #  Eqn (7)
            as = SP.an(P, i), SP.bn(P, i), SP.cn(P, i)
            @test x * pᵢ ≈ sum(a * p for (a, p) in zip(as, ps))

            αs = SP.αn(P, i), SP.βn(P, i), SP.γn(P, i)
            @test σ * dpᵢ ≈ sum(a * p for (a, p) in zip(αs, ps))

            âs = SP.ân(P, i), SP.b̂n(P, i), SP.ĉn(P, i)  #  Eqn 20
            @test pᵢ ≈ sum(a * dp for (a, dp) in zip(âs, Δps))

            αᴵs = SP.αᴵn(P, i), SP.βᴵn(P, i), SP.γᴵn(P, i)   # Eqn (16)
            @test x * Δpᵢ ≈ sum(a * dp for (a, dp) in zip(αᴵs, Δps))
        end
    end
end

@testset "Conversion" begin
    for (PP, Q) in ((Ps, Polynomial), (DPs, FallingFactorial))
        @testset for P in PP
            for _ in 1:10
                ps = rand(1.0:8, 5)
                p = P(ps)
                @test convert(P, convert(Q, p)) ≈ p
            end
        end
    end
    # through evaluation
    @testset for PP in (Ps, DPs)
        @testset for P in PP
            @testset for Q in PP
                p = P([1, 2, 3])
                x = variable(Q)
                @test p(x) ≈ convert(Q, p)
            end
        end
    end

    x = variable(Polynomial)
    # Connections: groupings have same sigma = a⋅x² + b⋅x + c
    @testset for Ps in (
        (
            Chebyshev,
            ChebyshevU,
            Gegenbauer{1 / 2},
            Gegenbauer{1},
            Jacobi{1 / 4,1 / 4},
            Jacobi{3 / 4,3 / 4},
        ),
        (Laguerre{0}, Laguerre{1 / 2}, Laguerre{1}),
        (Bessel{1 / 4}, Bessel{1 / 2}, Bessel{3 / 4}, Bessel{2}),
    )
        @testset for P in Ps
            @testset for Q in Ps
                @testset for n in 2:5
                    q = basis(Q, n)
                    # various ways to convert to P

                    @test q(variable(P))(variable()) ≈ q(variable())  # evaluation
                    @test SP._convert_cop(P, q)(variable()) ≈ q(variable())  # structural  equations
                    @test convert(P, convert(Polynomial, q))(variable()) ≈ q(variable()) # through Standard  basis
                    @test convert(P, q)(x) ≈ q(x)
                end
            end
        end
    end

    # Connections via FastTransforms
    x = variable()
    @testset for P in Ps
        @testset for Q in Ps
            p = P([1.0, 2.0, 3.0])
            @test convert(Q, p)(x) ≈ p(x)
        end
    end

    # Issue #43 conversion from Polynomials.ChebyshevT
    #  @test_throws ArgumentError convert(Legendre, ChebyshevT([0,1]))
end

@testset "Evaluation" begin
    @testset for P in Ps
        (SP.ismonic(P) || SP.isorthonormal(P)) && continue
        n, x = 4, 1.23
        p = basis(P, n)
        # compares  clenshaw  to  hypergeometric
        @test p(x) ≈ Basis(P, n)(x)
    end
end

@testset "Arithmetic" begin
    x = variable(Polynomial{Float64})

    @testset for P in union(Ps, DPs)
        pNULL = zero(P)
        p0, p1, p2, p3, p4 = ps = basis.(P, 0:4)

        as = [1, 2, 3, 4, 5]
        @test sum(a * pᵢ for (a, pᵢ) in zip(as, ps)) ≈ P(as)

        x = variable(Polynomial)
        @testset for p in ps
            @testset for q in ps
                @testset for op in (+, *)
                    @test op(p, q)(x) ≈ op(p(x), q(x))
                end
            end
        end

        for _ in 1:10
            p, q = P(rand(5)), P(rand(5))
            @testset for op in (+, *)
                @test (op(p, q))(x) ≈ op(p(x), q(x))
            end
        end
    end
end

@testset "Elementwise operations" begin
    @testset for P in union(Ps, DPs)
        xs = [1, 2, 3, 4]
        p = P(xs)
        @test -p == P(-xs)
        ys = copy(xs ./ 1) # might need float
        ys[1] += one(P)[0]
        @test p + 1 ≈ p + P(1)
        @test 2p == P(2xs)
    end
end

# too slow
@testset "Orthogonality" begin
    n = 5
    @testset for P in Ps
        B = Polynomials.basistype(P)
        a, b = extrema(domain(B))
        (isinf(a) || isinf(b)) && continue
        @testset for i in 2:n
            @testset for j in (i + 1):n
                val = SP.innerproduct(B, Polynomials.basis(P, i), Polynomials.basis(P, j))
                @test abs(val) <= 1e-4
            end
        end
    end

    @testset for P in (Krawchouk{1 / 2,n}, Hahn{1 / 4,1 / 2,n})
        @testset for i in 2:n
            @testset for j in (i + 1):n
                B = Polynomials.basistype(P)
                val = SP.innerproduct(B, Polynomials.basis(P, i), Polynomials.basis(P, j))
                @test abs(val) <= 1e-12
            end
        end
    end
end

@testset "divrem" begin
    @testset for P in union(Ps, DPs)
        ps, qs = rand(1:6, 4), rand(1:6, 3)
        p, q = P(ps), P(qs)
        a, b = divrem(p, q)
        @test p ≈ q * a + b
    end
end

@testset "Derivatives and integrals" begin
    @testset for P in Ps
        Q = Polynomial
        x = variable(Q)

        @testset for i in 1:6
            p = basis(P, i)
            @test derivative(p)(x) ≈ derivative(p(x))
        end

        p = P([1.0, 2, 3, 4, 3, 2, 1])
        q = convert(Polynomial, p)
        @test maximum(
            abs,
            [
                integrate(p, a, a + 1 / 2) - integrate(q, a, a + 1 / 2) for
                a in range(0, stop=1 / 2, length=10)
            ],
        ) <= sqrt(eps())

        for _ in 1:10
            ps = rand(1:10, 5)
            p = P(ps)
            @test (derivative ∘ derivative)(p) ≈ derivative(p, 2)

            q = (derivative ∘ integrate)(p) - p
            q = chop(q, atol=1e-8)
            @test degree(q) <= 0

            q = (integrate ∘ derivative)(p) - p
            q = chop(q, atol=1e-8)
            @test degree(q) <= 0
        end
    end

    # special  case 𝐍  cases
    𝐍 = 9
    @testset for P in (Krawchouk{1 / 2,𝐍}, Krawchouk{1 / 4,𝐍}, Hahn{1 / 2,1 / 2,𝐍})
        @testset for i in 0:𝐍
            p = basis(P, i)
            (p(𝐍) - p(0)) <= 1e-12 && continue
            @test sum(SP.∇ₓ(p)(j) for j in 1:𝐍) ≈ p(𝐍) - p(0)
            @test sum(SP.Δₓ(p)(j) for j in 0:(𝐍 - 1)) ≈ p(𝐍) - p(0)
        end
    end
end

@testset "roots" begin
    @testset for P in union(Ps, DPs)
        P <: Bessel && continue
        @testset for n in 6:2:10
            rts = roots(convert(Polynomial, basis(P, n)))
            evals = eigvals(SP.jacobi_matrix(P, n))
            @test maximum(abs, rts - evals) <= 1e-4
        end
    end

    # comrade matrix approach compared to conversion
    @testset for P in (
        Legendre,
        MonicLegendre,
        Chebyshev,
        MonicChebyshev,
        ChebyshevU,
        MonicChebyshevU,
    ) #(Legendre, Chebyshev, ChebyshevU)
        @testset for T in (Float64, Complex{Float64})
            ps = vcat(rand(T, 4), 1)
            p = P(ps)
            q = convert(Polynomial, p)
            λs = sort(norm.(roots(p)))
            γs = sort(norm.(roots(q)))
            @test maximum(abs, λs - γs) < sqrt(eps())
        end
    end

    # comrade pencil
    @testset for P in (Legendre, Chebyshev, ChebyshevU)
        p = P(rand(5))
        λ = rand()
        C₀, C₁ = SP.comrade_pencil(p)
        M = λ * C₁ - C₀
        L, Ũ = SP.comrade_pencil_LU(p)(λ)
        @test det(Ũ) ≈ 1
        Ũ[end, end] *= p(λ)
        @test L * Ũ ≈ M
    end
end

@testset "fitting" begin
    f(x) = exp(-x) * cospi(x)

    # Interpolating
    @testset for P in Ps
        P <: SP.AbstractOrthogonalPolynomial || continue
        dom = domain(P)
        (isinf(first(dom)) || isinf(last(dom))) && continue
        q = fit(P, f, 10)
        @test maximum(abs, q(x) - f(x) for x in range(0, stop=1 / 2, length=10)) <= 1e-4
    end

    # least squares
    @testset for P in (Chebyshev,)
        q = fit(Val(:lsq), P, f, 10)
        @test maximum(abs, q(x) - f(x) for x in range(0, stop=1 / 2, length=10)) <= 1e-4
    end

    # series
    @testset for P in (Chebyshev,)
        q = fit(Val(:series), P, f)
        @test maximum(abs, q(x) - f(x) for x in range(0, stop=1 / 2, length=10)) <=
              sqrt(eps())
    end
end

@testset "quadrature" begin
    f(x) = x^7
    n = 4

    @testset for P in Ps
        P <: SP.AbstractOrthogonalPolynomial || continue
        B = Polynomials.basistype(P)
        !all(isfinite.(extrema(domain(P)))) && continue
        q = sum(f(tau) * w for (tau, w) in zip(SP.gauss_nodes_weights(B, n)...))
        p = SP.innerproduct(B, f, one)
        @test abs(p - q) <= 10sqrt(eps(T))
    end
end

@testset "Showing" begin
    @test sprint(show, Hermite([1, -2, 3])) == "Hermite(1⋅H₀(x) - 2⋅H₁(x) + 3⋅H₂(x))"

    ps = [1.0 + 2.0im, -0.5, 3.7im, 0.1 + 0.01im]     # issue 26 with complex coefficients
    p = Hermite(ps)
    @test sprint(show, p) ==
          "Hermite((1.0 + 2.0im)⋅H₀(x) - 0.5⋅H₁(x) + 3.7im⋅H₂(x) + (0.1 + 0.01im)⋅H₃(x))"
end

# issue 44 allows fully typed array construction, as `N` parameter is removed
@testset "array elements" begin
    @testset for P in Ps
        p, q = basis.(P, (2, 3))
        @test typeof([p]) == typeof([p, q])
    end

    q = basis(LaurentPolynomial, 3) # not Polynomial v4.0
    @testset for P in Ps
        p = basis(P, 2)
        @test typeof([q]) == typeof([p, q])
    end
end

@testset "Shifted" begin
    for (P, Q) in ((Legendre, ShiftedLegendre), (Jacobi{1/2,1/2}, ShiftedJacobi{1/2,1/2}))
        p, q = basis(P, 4), basis(Q, 4)
        x = Polynomial(:x)
        ϕ(x) = 2x - 1

        @test p(ϕ(x)) ≈ q(x)
        @test p(ϕ(0.5)) ≈ q(0.5)
    end
end
