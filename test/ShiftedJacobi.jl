using SpecialPolynomials: Pochhammer
using Polynomials

# compute <>ᵅᵝ inner product
# also used in DualBernstein tests.
function ip(f,g,α,β; n=100)
    SpecialPolynomials.innerproduct(ShiftedJacobi{α,β}, f, g)
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
        Rn₁ = Basis{SpecialPolynomials.ShiftedJacobiBasis{α, β}}(n)(x)
        Rn₂ = basis(ShiftedJacobi{α, β},n)(x)
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
