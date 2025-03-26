using SpecialPolynomials: Pochhammer
using Polynomials

# compute <>ᵅᵝ inner product
function ip(f,g,α,β; n=100)
    simpsons  = (f,a,b) -> (c = a/2 + b/2;(1/6) * (f(a) + 4*f(c) + f(b)))
    λ = x -> (1-x)^α*x^β * f(x) * g(x)
    xs = range(0, 1, n+1)
    xs′ = zip(Iterators.take(xs, n), Iterators.drop(xs, 1))
    sum(simpsons(λ, xᵢ₋₁, xᵢ) * (xᵢ-xᵢ₋₁) for (xᵢ₋₁, xᵢ) ∈ xs′)
end

@testset "Construction" begin
    α,β = 1/2, 1/2
    R = ShiftedJacobi{α, β}
    p = R([1,2,3])
    @test p(.2) ≈ 1/40
end

@testset "Other Construction" begin

    # from 1.2 of https://arxiv.org/abs/2004.09801
    # check alternate construcdtion
    for (α,β) ∈ ((0,0), (1/2,1/2), (1,1))
        n = 3
        Polynomials.@variable x
        R = ShiftedJacobi{α, β}
        Rn₁ = SpecialPolynomials.classical_hypergeometric(R,n,x)
        Rn₂ = basis(R,n)(x)
        @test Rn₁ ≈ Rn₂
    end

end

@testset "Orthogonality" begin
    # check Remark 1.1 of https://arxiv.org/abs/2004.09801
    for (α,β) ∈ ((0,0), (1/2,1/2), (1,1))
        σ = α + β + 1

        P = ShiftedJacobi{α,β}
        i,j = 3,4
        Ri,Rj = basis(P,i), basis(P,j)
        a = ip(Ri, Rj, α, β)
        @test a ≈ 0 atol = 1e-5

        a = ip(Ri, Ri, α, β)
        K = gamma(α+1)*gamma(β+1)/gamma(σ+1)
        k = i
        hₖ = K * Pochhammer(α+1,k)*Pochhammer(β+1,k)/factorial(k)/ (2k/σ + 1) /Pochhammer(σ, k)
        @test a ≈ hₖ atol = 1e-3
    end
end
