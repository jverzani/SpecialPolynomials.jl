using Polynomials
using SpecialFunctions
using SpecialPolynomials
using SpecialPolynomials: Pochhammer

@testset "Construction" begin
    n,α,β = 5,0,0
    D = DualBernstein{n,α,β}
    p = D([0,0,1,0,0,0])
    @test p(.5) ≈ 9.75
end

@testset "Alternate construction" begin
    # Corollary 2.2
    n = 5
    α, β = 1/2, 1/2
    D = DualBernstein{n,α,β}
    x = 0.2
    for i in 0:n
        D′ = SpecialPolynomials.classical_hypergeometric(D,i,x)
        @test D′ ≈ basis(DualBernstein{n,α,β},i)(x)
    end
end


@testset "Orthogonality" begin
    # check Remark 1.1 of https://arxiv.org/abs/2004.09801
    ϵ = 1e-6

    for (α,β) ∈ ((0,0), (1/2,1/2), (1,1))
        n = 5
        i,j = 2, 3
        D = DualBernstein{n, α, β}
        Di, Dj = basis.(D, (i,j))
        Bi, Bj = basis.(Bernstein{n}, (i,j))
        R = ShiftedJacobi{α, β}
        dij = ip(Bi,Dj, α, β)
        @test dij ≈ 0 atol = ϵ
        dij = ip(Bj,Di, α, β)
        @test dij ≈ 0 atol = ϵ
        dii = SpecialPolynomials.innerproduct(R,Bi,Di)
        @test dii ≈ 1 atol=ϵ
    end
end

#= too unreliable
@testset "O(n)" begin
    for (α, β) ∈ ((0,0), (1/2, 1/2), (1, 1))
        ts = zeros(6)
        for (i, n) ∈ enumerate([10^i for i in 1:6])
            Dᵢ = basis(DualBernstein{n,α,β}, n÷2)
            Dᵢ(rand())
            ts[i] = @elapsed Dᵢ(rand()) # noisier but quicker than @belapsed
        end
        Δs = ts[2:end] ./ ts[1:end-1]
        x̄ = sum(Δs)/length(Δs)
        @test x̄ <= 25 # should be about 10
    end
end
=#

@testset "minimize least square error" begin
    f(x) = sinpi(x)
    n = 5
    α, β = 1/2, 1/2
    D = DualBernstein{n,α,β}
    R = ShiftedJacobi{α, β}
    Iₖ = [ip(f, basis(D,k), α, β) for k in 0:n]


    Jₖ = [ip(f, basis(R,k), α, β) for k in 0:n] # say, just some comparison

    B = Bernstein{n}
    pn = B(Iₖ)
    qn = B(Jₖ)
    Fp, Fq = x -> f(x) - pn(x), x -> f(x) - qn(x)
    @test ip(Fp, Fp, α, β) < ip(Fq, Fq, α, β)
end
