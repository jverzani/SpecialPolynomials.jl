## Comrade Matrix
## The comrade matrix plays the role of the companion matrix for static basis polynomials

#
# generate the  matrix for an orhtogonal polynomial family
# should have eigvals(comrade_matrix(p)) ≈ roots(convert(Polynomial, p))
# This isn't exported
function comrade_matrix(p::P) where {P <: AbstractCCOP}
    U,V = comrade_pencil(p)
    U * inv(V)
end

# Generate the [comrade pencil](http://www.cecm.sfu.ca/personal/pborwein/MITACS/papers/LinMatPol.pdf)
# C₀, C₁ with λ C₁ - C₀ the linearization; and C₀ ⋅ C₁⁻¹ the companion matrix
# can be used with square matrix coefficients
function comrade_pencil(p::P, symmetrize=false) where {T, P <: AbstractCCOP{T}}
    comrade_pencil(p, αβγk(p)..., symmetrize)
end

function αβγk(p::P) where {T, P <: AbstractCCOP{T}}
    n = Polynomials.degree(p)
    𝐏 = _comrade_pencil_type(p)

    As = An.(𝐏, 0:n-1)
    Bs = Bn.(𝐏, 0:n-1)
    Cs = Cn.(𝐏, 1:n-1)
    kₙ, kₙ₋₁ = leading_term(𝐏{T}, n), leading_term(𝐏{T},n-1)

    # α = 1/A
    # β = -B/A
    # γ = C/A


    α = inv.(As[1:end-1])
    β = (-Bs ./ As)
    γ = (Cs ./ As[2:end])

    (α,β,γ,kₙ₋₁, kₙ)
end
_comrade_pencil_type(p::P) where {T, P <: AbstractCCOP{T}} =  ⟒(P){T}
_comrade_pencil_type(p::P) where {T, M<: AbstractMatrix{T}, P <: AbstractCCOP{M}} = ⟒(P){T}

# α = [α₀, α₁, ⋯, αₙ₋₂]
# β = [β₀, β₁, ⋯, βₙ₋₂, βₙ₋₁]
# γ = [γ₁, γ₂ ⋯, αₙ₋₁]
# kₙ₋₁/kₙ = aₙ₋₁
function comrade_pencil(p::P, α, β, γ, kₙ₋₁, kₙ, symmetrize=false) where {T, P <: Polynomials.AbstractPolynomial{T}}

    n = degree(p)
    R = eltype(one(T)/1)
    C₀ = diagm(2 => zeros(R, n-2),
               1 => γ,
               0 => β,
               -1 => α
               )

    C₀[end-1,end] *= kₙ * p[end]
    C₀[end,end] *= kₙ * p[end]

    for i ∈ 0:n-1
        C₀[1+i, end] -= kₙ₋₁ * p[i]
    end

    C₁ = diagm(0 => ones(R, n))
    C₁[end,end] = kₙ * p[end]

    C₀, C₁
end


function comrade_pencil(p::P, α, β, γ, kₙ₋₁, kₙ, symmetrize=false) where {T, M <: AbstractMatrix{T}, P <: Polynomials.AbstractPolynomial{M}}

    n = Polynomials.degree(p)
    m,m′ = size(p[0])
    @assert m == m′ # square
    𝐈 = I(m)

    R = eltype(one(T)/1)
    C₀ = zeros(R, n*m, n*m)
    C₁ = diagm(0 => ones(R, n*m))
    Δ = 1:m

    C₀[0*m .+ Δ, 0*m .+ Δ] .= β[1]
    C₀[1*m .+ Δ, 0*m .+ Δ] .= α[1]

    for i ∈ 1:n-2
        C₀[(i-1)*m .+ Δ, i*m .+ Δ] .= γ[i] * 𝐈
        C₀[i*m .+ Δ, i*m .+ Δ] .= β[1+i] * 𝐈
        C₀[(i+1)*m .+ Δ, i*m .+ Δ] .= α[1+i] * 𝐈
    end

    C₀[(n-2)*m .+ Δ, (n-1)*m .+ Δ] .= γ[end] * kₙ * p[end] * 𝐈
    C₀[(n-1)*m .+ Δ, (n-1)*m .+ Δ] .= β[end] * kₙ * p[end] * 𝐈

    for i ∈ 0:n-1
        C₀[i*m .+ Δ, (n-1)*m .+ Δ] .-= kₙ₋₁ * p[i]
    end



    C₁[(n-1)*m .+ Δ, (n-1)*m .+ Δ] .=  kₙ * p[end]

    if symmetrize
        H₀ = zeros(T, n*m, n*m)
        for k ∈ 1:n
            for j ∈ 1:k
                H₀[(k-j)*m .+ Δ ,(j-1)*m .+ Δ] .= p[k]
            end
        end
        return (C₀*inv(C₁)) * H₀, H₀
    else
        return C₀, C₁
    end
end

"""
    comrade_pencil_LU(p)

Return an LU decomposition of `λC₁ - C₀`, where `C₀` and `C₁` are
the pieces of the comrade pencil.

The returned value is a function. We have `L,Ũ = comrade_pencil_LU(p)(λ)`.
The upper-triangular part, `Ũ` is not quite `U`, as `U` is singluar at eigenvalues. Rather,
this holds up to floating point

```
C₀, C₁ = comrade_pencil(p)
L,Ũ = comrade_pencil_LU(p)(λ)
@assert det(Ũ) ≈ 1
Ũ[end,end] *= p(λ)
L*Ũ - (λ*C₁ - C₀) ≈ 0 # if atol=sqrt(eps())
```
"""
comrade_pencil_LU(p::P) where {P <: AbstractCCOP} =  comrade_pencil_LU(p, αβγk(p)...)

function comrade_pencil_LU(p::P, α, β, γ, kₙ₋₁, kₙ) where {T, P <: Polynomials.AbstractPolynomial{T}}
    n = degree(p)
    R = typeof(one(T)/1)
    # return a function of λ
    λ -> begin
        𝐏 = _comrade_pencil_type(p)
        ϕs = [basis(𝐏,i)(λ) for i ∈ 0:n-1]
        L = zeros(R, n, n)
        U = zeros(R, n, n)
        for i ∈ 1:n-1
            L[i,i] = one(R)
            L[i+1,i] = -ϕs[i]/ϕs[i+1]

            U[i,i+1] = -γ[i]
            U[i,i] = α[i] * ϕs[i+1]/ϕs[i]
        end
        L[end,end] = one(R)

        U[1,end] = kₙ₋₁ * p[0]
        for i ∈ 2:(n-2)
            j = 1+i
            U[i,end] = kₙ₋₁ * p[j-2] + ϕs[j-2]/ϕs[j-1] * U[i-1,end]
        end

        U[n-1,end] = kₙ₋₁ * p[end-2] + ϕs[end-2]/ϕs[end-1]*U[n-2,end] - kₙ*γ[end]*p[end]
        U[end,end] = ϕs[1]/ϕs[end] / prod(α) * 1 # p(λ)
        (L=L, U=U)
    end
end


function comrade_pencil_LU(p::P, α, β, γ, kₙ₋₁, kₙ) where {T, M <: AbstractMatrix{T}, P <: Polynomials.AbstractPolynomial{M}}

    m, m′ = size(p[0])
    @assert m == m′
    Iₘ = I(m)
    Δ = 1:m

    n = degree(p)
    R = typeof(one(T)/1)

    # return a function of λ
    λ -> begin
        𝐏 = _comrade_pencil_type(p)
        ϕs = [basis(𝐏,i)(λ) for i ∈ 0:n-1]
        L = zeros(R, n*m, n*m)
        Ũ = zeros(R, n*m, n*m)
        for i ∈ 1:n-1
            i′ = i - 1
            L[i′*m .+ Δ, i′*m .+ Δ] .= Iₘ
            L[(i′+1)*m .+ Δ, i′*m .+ Δ] .= -ϕs[i]/ϕs[i+1] * Iₘ

            Ũ[i′*m .+ Δ,(i′+1)*m .+ Δ] .= -γ[i] * Iₘ
            Ũ[i′*m .+ Δ, i′*m .+ Δ] .= α[i] * ϕs[i+1]/ϕs[i] * Iₘ
        end
        END = (n-1)*m .+ Δ
        L[END, END] .= Iₘ

        Ũ[Δ, END] .= kₙ₋₁ * p[0]
        for i ∈ 2:(n-2)
            j = 1+i
            i′ = i - 1
            Ũ[i′*m .+ Δ, END] .= kₙ₋₁ * p[j-2] + ϕs[j-2]/ϕs[j-1] * Ũ[(i′-1)*m .+ Δ,END]
        end

        Ũ[(n-1-1)*m .+ Δ,END] = kₙ₋₁ * p[end-2] +
            ϕs[end-2]/ϕs[end-1]*Ũ[(n-2-1)*m .+ Δ,END] - kₙ*γ[end]*p[end]


        Ũ[END,END] = ϕs[1]/ϕs[end] / prod(α) * Iₘ # use Ũ, not U, which uses p(λ)
        (L=L, U=Ũ)
    end

end



## -----

## Implement algorithm for eigenvalue of a comrade matrix from:
## Fast computation of eigenvalues of companion, comrade, and related matrices
## Jared L. Aurentz · Raf Vandebril · David S. Watkins

## Overlap with AMRVW.jl, but not really, as this paper uses different transforms:
## * a general 2x2 core transformation
## * a GaussTransform of the type `[1 0; m 1]`

abstract type AbstractTransform{T} end
# Ci CoreTransform [c 1; 1 0]
struct CoreTransform{T} <: AbstractTransform{T}
    g11::T
    g12::T
    g21::T
    g22::T
end
CoreTransform(m::Matrix) = CoreTransform(m[1,1], m[1,2], m[2,1], m[2,2])
LinearAlgebra.norm(c::CoreTransform) = maximum(abs, splat(c))
Base.isnan(c::CoreTransform) = isnan(c.g12) || isnan(c.g22) || isnan(c.g21) || isnan(c.g11)


# Gi GaussTransform [1 0, m 1]
struct GaussTransform{T} <: AbstractTransform{T}
    m::T
end

splat(c::CoreTransform) = (c.g11, c.g12, c.g21, c.g22)
splat(c::GaussTransform{T}) where {T}  = (one(T), zero(T), c.m, one(T))

function Base.inv(c::CoreTransform)
    b11,b12,b21,b22 = splat(c)
    Δ = b11 * b22 - b12 * b21
    CoreTransform(b22/Δ, -b12/Δ, -b21/Δ, b11/Δ)
end

Base.inv(g::GaussTransform) = GaussTransform(-g.m)


# overwrite R
# non-allocating would be better!
function LinearAlgebra.lmul!(c::AbstractTransform, R, i)
    ri = copy(R[i,:])
    c11,c12,c21,c22 = splat(c)
    R[i,:] = R[i,:] * c11 + R[i+1,:] * c12
    R[i+1,:] = ri * c21 + R[i+1,:] * c22
    nothing
end

## ----

# Helpful for developing
function Base.show(io::IO, c::AbstractTransform)
    c11,c12,c21,c22 = splat(c)
    show(IOContext(io, :compact => true), "text/plain", [c11 c12;
                                                         c21 c22])
    println(io)
end


## realize transform as a matrix
function realize(c::AbstractTransform{T}, i, n) where {T}
    A = diagm(0 => ones(T, n))
    a,b,c,d = splat(c)
    A[i,i] = a
    A[i,i+1] = b
    A[i+1,i] = c
    A[i+1,i+1] = d
    A
end

function Base.Matrix(Cs, B)
    n = length(Cs)
    m = size(B)[1]
    for i ∈ n:-1:1
        B = realize(Cs[i], i,m)*B
    end
    B
end


## ----- primary operations with CoreTransforms

# type 1
function turnover(Ci::CoreTransform,Ci1::CoreTransform,Gi::GaussTransform)
    a11,a12,a21,a22 = splat(Ci)
    b22,b23, b32, b33 = splat(Ci1)
    g21 = Gi.m

    â11 = a11 + a12 * (b22*g21)
    â12 = a12
    â21 = a21 + a22 * (b22 * g21)
    â22 = a22

    ĝ32 = b32 * g21 / â21

    b̂22 = b22
    b̂23 = b23
    b̂32 = b32 - (ĝ32 * a22) * b22
    b̂33 = b33 - (ĝ32 * a22) * b23

    GaussTransform(ĝ32), CoreTransform(â11,â12,â21,â22), CoreTransform(b̂22,b̂23,b̂32,b̂33)

end

# type 2
function turnover(Gi::GaussTransform, Gi1::GaussTransform, Ci::CoreTransform)
    a21 = Gi.m
    b32 = Gi1.m
    g11,g12,g21,g22 = Ci.g11,Ci.g12,Ci.g21,Ci.g22

    â11 = g11
    â12 = g12
    â21 = g21 + a21 * g11
    â22 = g22 + a21 * g12

    ĝ32 = b32*g21/â21
    b̂32 = b32*g22 - ĝ32*â22

    GaussTransform(ĝ32), CoreTransform(â11, â12, â21, â22), GaussTransform(b̂32)

end

# type 3
function turnover(Gi::GaussTransform, Gi1::GaussTransform, G′i::GaussTransform)
    a21 = Gi.m
    b32 = Gi1.m
    g21 = G′i.m

    â21 = a21 + g21
    ĝ32 = b32 * g21 / â21
    b̂32 = b32 - ĝ32

    GaussTransform(ĝ32), GaussTransform(â21), GaussTransform(b̂32)
end

function fuse(c::CoreTransform, g::GaussTransform)
    a11,a12,a21,a22 = splat(c)
    m = g.m
    CoreTransform(a11 + a12 * m, a12,
                  a21 + a22 * m, a22)
end

function fuse(g::GaussTransform, c::CoreTransform)
    m = g.m
    a11,a12,a21,a22 = splat(c)

    CoreTransform(a11,         a12,
                  a11*m + a21, a12*m + a22)
end

function fuse(g1::GaussTransform, g2::GaussTransform)
    GaussTransform(g1.m + g2.m)
end

# Rlᵢ -> l′ᵢR′
# modifies R; returns l′ᵢ
function passthrough!(R, l, i)
    m = l.m

    m′ = -m * R[i+1,i+1]/(R[i,i] + R[i,i+1]*m)

    for j ∈ 1:i+1
        R[j,i] += m * R[j,i+1] # R[:,i]   += m * R[:,i+1]  # R * Gᵢ
    end
    for j ∈ i:size(R)[2]
        R[i+1, j]  += m′ * R[i,j] # R[i+1,:] += m′ * R[i,:] # Ĝᵢ * R
    end

    R[i+1,i] = 0 # may accumulate round off

    GaussTransform(-m′)

end

## ---- initial shift

# get eigenvalues of lower 2x2 matrix of A = Cs * B
function lower_eigs(Cs, R::Matrix{T}, n=length(Cs)) where {T}

    cn11, cn12, cn21, cn22 = splat(Cs[n])

    if n >= 2
        cn1_11, cn1_12, cn1_21, cn1_22 = splat(Cs[n-1])

        u11 = cn1_21 * R[n-1,n] + cn1_22 * cn11 * R[n,n]
        u12 = cn1_21*R[n-1,n+1] + cn1_22 * cn11 * R[n,n+1] + cn1_22 * cn12 * R[n+1, n+1]
        u21 = cn21 * R[n,n]
        u22 = cn21 * R[n,n+1] + cn22 * R[n+1,n+1]

    else
        u11 = cn11 * R[1,1]
        u12 = cn11 * R[1,2] + cn12 * R[2,2]
        u21 = zero(T)
        u22 = cn22 * R[2,2]
    end

    ρ1, ρ2 = eigen_values(u11, u12, u21, u22)
    ρ1, ρ2
end


# We create bulge through Â L⁻¹ A L and chase it down through as sequence of
# unitary operations of the type G⁻¹ Â G, G a Gauss transform

L(Cs, R::Matrix{T}, random::Bool=false) where {T <: Complex} = L1(Cs, R, random)
L(Cs, R::Matrix{T}, random::Bool=false) where {T}  = L1L2(Cs, R, random)

function L1(Cs, R::Matrix{T}, random::Bool=false) where {T}
    random && return GaussTransform(rand(T))

    c11,c12,c21,c22 = splat(Cs[1])
    r11 = R[1,1]
    a11, a21 = c11*r11, c21*r11

    ρs = lower_eigs(Cs, R)
    ρ1 = norm(a21 - ρs[1]) < norm(a21 - ρs[2]) ? ρs[1] : ρs[2]

    m = a21 / (a11 - ρ1)

    (isinf(m) || isnan(m)) && return GaussTransform(rand(T))

    GaussTransform(m)
end

function L1L2(Cs,R::Matrix{T}, random::Bool=false) where {T}
    random && return (GaussTransform(rand(T)), GaussTransform(rand(T)))

    ρ1, ρ2 = lower_eigs(Cs, R)

    ρ1ρ2 = real(ρ1*ρ2)
    ρ1_ρ2 = real(ρ1 + ρ2)


    cn11,cn12,cn21,cn22 = splat(Cs[2])
    cn1_11,cn1_12, cn1_21, cn1_22 = splat(Cs[1])

    a11 = cn1_11 * R[1,1]
    a12 = cn1_11 * R[1,2] + cn1_12 * cn11 * R[2,2]
    a21 = cn1_21 * R[1,1]
    a22 = cn1_21 * R[1,2] + cn1_22 * cn11 * R[2,2]
    a32 = cn21 * R[2,2]

    x1 = a11^2 +  a12 * a21 - ρ1_ρ2*a11 + ρ1ρ2
    x2 = a21 * a11 + a22 * a21 - ρ1_ρ2*a21
    x3 = a32 * a21
    m1 = x2/x1
    m2 = x3/x2

    (isinf(m1) || isinf(m2) || isnan(m1) || isnan(m2)) && return (GaussTransform(rand(T)), GaussTransform(rand(T)))

    # return l1, l2 (L = L2 * L1 though)
    GaussTransform(m1), GaussTransform(m2)

end


# single-shift bulge chase
# updates Cs, R
function chase_bulge!(Cs, R, l1::GaussTransform)
    n = length(Cs)
    Cs[1] = fuse(inv(l1), Cs[1])
    for i ∈ 1:n-1
        l1 = passthrough!(R, l1, i)
        l1, Cs[i], Cs[i+1] = turnover(Cs[i],Cs[i+1], l1)
    end
    l1 = passthrough!(R, l1, n)
    Cs[n] = fuse(Cs[n], l1)
    nothing
end

# double-shift bulge chase
function chase_bulge!(Cs, R, l1l2::Tuple)
    n = length(Cs)

    l1, l2 = l1l2
    l1′, l2′ = inv(l1), inv(l2)
    l1′, Cs[1], l2′ = turnover(l1′, l2′, Cs[1])
    Cs[2] = fuse(l2′, Cs[2])

    for i in 1:n-2

        l2 = passthrough!(R, l2, i+1)
        l2, Cs[i+1], Cs[i+2] = turnover(Cs[i+1], Cs[i+2], l2)

        l1 = passthrough!(R, l1, i)
        l1, Cs[i], Cs[i+1] = turnover(Cs[i], Cs[i+1], l1)

        l2, l1, l1′ = turnover(l1′, l2, l1)

    end
    # knit in ends
    l2 = passthrough!(R, l2, n)
    Cs[n] = fuse(Cs[n], l2)


    l1 = passthrough!(R, l1, n-1)
    l1, Cs[n-1], Cs[n] = turnover(Cs[n-1], Cs[n], l1)


    l1 = fuse(l1′, l1)
    l1 = passthrough!(R, l1, n)
    Cs[n] = fuse(Cs[n], l1)

    nothing
end


# use AVW algorithm to find eigen value of matrix
# decomposed into Cs::Vector{GaussTransform}, B::UpperTriangular
function AVW_eigvals(Cs, B)



    max_ctr = 1
    ctr = 0
    while length(Cs) > 1

        # get L # random if no deflation
        if 0 < ctr < 15
            ctr += 1
            l = L(Cs, B)
        else
            l = L(Cs, B, true) # random=true
            ctr = 1
        end

        chase_bulge!(Cs, B, l)
        deflateQ!(Cs, B) && (ctr = 0)

        # check runaway
        max_ctr += 1
        (!isempty(Cs) && isnan(Cs[end])) && break
        max_ctr > 1000 && error("Failure to converge, too many steps")

    end

    !isempty(Cs) && lmul!(Cs[1], B, 1)

    _eigvals(B)

end

# modifies Cs and B on deflation
# return true if deflation
# we deflate when off-diagonal entry < ϵ
# paper would have |aₙ₊₁,ₙ| ≤ eps() * (|aₙ,ₙ| + |aₙ₊₁,ₙ₊₁|)
# but we just use eps()^(2/3) which seems to work better for Float64
function deflateQ!(Cs, B::Matrix{T}) where {T}
    # simple case
    # is aₙ₊₁,ₙ < eps(T)(|ann| + |an+1,n+1|)
    n = length(Cs)

    c_11, c_12, c_21, c_22 = splat(Cs[n])
    c1_11, c1_12, c1_21, c1_22 = n > 1 ? splat(Cs[n-1]) : (zero(T), zero(T), zero(T), zero(T))


    Rn_1n = n > 1 ? B[n-1,n] : zero(T)

    an1_n = c_21 * B[n,n]
    ann = c1_21 * Rn_1n + c1_22 * c_11 * B[n,n]
    an1n1 = c_21 * B[n,n+1] + c_22 * B[n+1,n+1]

    ϵ = cbrt(eps(real(T))^2)

    #if norm(an1_n) <=  ϵ * (norm(ann) + norm(an1n1))
    if norm(an1_n) <=  ϵ
        n = length(Cs)
        c = pop!(Cs)
        lmul!(c, B, n)
        return true
    end

    !(T <: Real) && return false
    n == 1 && return false

    c2_11, c2_12, c2_21, c2_22 = n > 2 ? splat(Cs[n-2]) : (zero(T), zero(T), zero(T), zero(T))
    Rn_1n_1 = n > 1 ? B[n-1, n-1] : zero(T)
    Rn_2n_1 = n > 2 ? B[n-2, n-1] : zero(T)
    ann_1 = c1_21 * Rn_1n_1
    an_1n_1 = c1_11 * c2_22 * Rn_1n_1 + c2_21 * Rn_2n_1

    # if norm(ann_1) <=  eps(real(T)) * (norm(an_1n_1) + norm(ann))
    if norm(ann_1) <=  ϵ
        A = Matrix(Cs, B)
        n = length(Cs)
        c = pop!(Cs)
        lmul!(c, B, n)
        c = pop!(Cs)
        lmul!(c, B, n-1)
        return true
    end

    return false
end

# once algorithm is complete, this grabs eigenvalues from B
function _eigvals(B::Matrix{T}) where {T}
    m = size(B)[1]
    es = Vector{complex(T)}(undef, m)
    i = 1

    #ϵ = 1000 * eps(real(T))
    ϵ = cbrt(eps(real(T))^2)
    while i < m
        #if abs(B[i+1,i]) <= ϵ * (norm(B[i,i]) + norm(B[i+1,i+1]))
        if abs(B[i+1,i]) <= ϵ
            es[i] = B[i,i]
            i += 1
        else
            r1, r2 = eigen_values(B[i,i], B[i,i+1], B[i+1,i], B[i+1,i+1])
            es[i] = r1
            es[i+1] = r2
            i += 2
        end
    end

    i == m && (es[end] = B[end,end])

    LinearAlgebra.sorteig!(es)
    es
end

# one step improves accuracy
function newton_refinement(λs, p::P)  where {P <: AbstractCCOP}
    p′ = derivative(p)
    λs .- p.(λs) ./ p′.(λs)
end


## -----
## from AMRVW
# altrnative to eigvals that works with other types
function  eigen_values(a11::T, a12::T, a21::T, a22::T) where {T <: Real}
    b = (a11 + a22) / 2
    c = a11 * a22 - a12 * a21

    e1r, e1i, e2r, e2i = qdrtc(b,c) #qdrtc(one(b), b, c)
    return  complex(e1r, e1i), complex(e2r, e2i)

end

# from `modified_quadratic.f90`
function  eigen_values(a11::S, a12::S, a21::S, a22::S) where {S <: Complex}

    tr = a11 + a22
    detm = a11 * a22 - a21 * a12
    disc::S = sqrt(tr * tr - 4.0 * detm)

    u::S = abs(tr + disc) > abs(tr - disc) ? tr + disc : tr - disc
    if iszero(u)
        return zero(S), zero(S)
    else
        e1::S = u / 2.0
        e2::S = detm / e1
        return e1, e2
    end

end

## Kahan quadratic equation with fma
##  https://people.eecs.berkeley.edu/~wkahan/Qdrtcs.pdf

## solve ax^2 - 2bx + c
function qdrtc(a::T, b::T, c::T) where {T <: Real}
    # z1, z2 roots of ax^2 - 2bx + c
    d = discr(a,b,c)  # (b^2 - a*c), as 2 removes 4

    if d <= 0
        r = b/a  # real
        s = sqrt(-d)/a #imag
        return (r,s,r,-s)
    else
        r = sqrt(d) * (sign(b) + iszero(b)) + b
        return (r/a, zero(T), c/r, zero(T))
    end
end

## more work could be done here.
function discr(a::T,b::T,c::T) where {T}
    pie = 3.0 # depends on 53 or 64 bit...
    d = b*b - a*c
    e = b*b + a*c

    pie*abs(d) > e && return d

    p = b*b
    dp = muladd(b,b,-p)
    q = a*c
    dq = muladd(a,c,-q)

    (p-q) + (dp - dq)
end

function qdrtc(b::T, c::T) where {T <: Real}
    d = b*b - c
    if d <= 0
        r = b  # real
        s = sqrt(-d) #imag
        return (r,s,r,-s)
    else
        r = sqrt(d) * (sign(b) + iszero(b)) + b
        return (r, zero(T), c/r, zero(T))
    end
end

## -----
# In SpecialPolynomials we have P_{n+1} = (An*x + Bn) * P_n + Cn P_{n-1}
# In paper, we have
# αₖ pₖ = (z  - βₖ) ⋅ pₖ₋₁ - γₖ ⋅ pₖ₋₂
## Decompose p::P into Cs, B where
## * Cs -- a vector of Gauss Transfroms
## * B -- an Upper triangular matrix
## Matrix(Cs, B) is comrade matrix
function comrade_decomposition(P::Type{<:AbstractCCOP}, ps)

    n = length(ps) - 1
    As = An.(P, 0:n-1)
    Bs = Bn.(P, 0:n-1)
    Cs = Cn.(P, 1:n-1)
    kₙ, kₙ₋₁ = leading_term(P, n), leading_term(P,n-1)

    p̂s = copy(ps) ./ 1

    # αs = [1/A0, 1/A1, ..., 1/An-2]
    # βs = [B0/A0, ..., Bn-1/An-1]
    # γs = [C1/A1, ..., Cn-1/An-1]
    # cs [ 0 ... 0 γ[end] β[end]] - ps[1:end-1]/ps[end] / An-1


    aₙ = An(P,n)
    βs = -Bs ./ As
    βn = pop!(βs)
    γs = [Cs[i]/As[i+1] for i ∈ 1:length(Cs)]
    γn = pop!(γs)
    pop!(As)
    αs = 1 ./ As

    p̂n = pop!(p̂s)
    p̂s ./= p̂n
    p̂s *= -kₙ₋₁ / kₙ
    p̂s[end-1] += γn #* kₙ
    p̂s[end] += βn #* kₙ

    comrade_decomposition(p̂s, αs, βs, γs)

end

# The comrade matrix is the following (though the paper flips it for Gaussian elimination)
# [β₁ γ₂ .  .  c₀;
#  α₁ β₂ γ₃ .  c₁;
#  .  α₂ β₃ γ₄ c₂;
#  .  .  α₃ β₄ c₃;
#  .  .  .  α₄ c₄]
# returns Cs, B with eigvals(C1 ⋅ C2 ⋅ ⋯ ⋅ Cn * B) = roots(P(p))
function comrade_decomposition(cs::Vector{T}, αs, βs, γs) where {T}


    n = length(cs)
    R = zeros(T, n, n)

    CTs = Vector{CoreTransform{T}}(undef, n-1)

    ĉ₀, ĉ₁ = cs[end+1-1], cs[end+1-2]

    for i ∈ 1:n-2
        λ = ĉ₀ / αs[end+1-i]
        ĉ₀ = ĉ₁ - λ * βs[end+1-i]
        ĉ₁ = cs[end+1 - (i+2)] - λ * γs[end+1-i]
        CTs[i] = CoreTransform(λ, one(T), one(T), zero(T))
    end
    λ = ĉ₀ / αs[1]
    CTs[end] = CoreTransform(λ, one(T), one(T), zero(T))
    ĉ₀ = ĉ₁ - λ * βs[1]

    for i ∈ 1:n-2
        R[i, i]   = αs[end+1-i]
        R[i, i+1] = βs[end+1-i]
        R[i, i+2] = γs[end+1-i]
    end

    R[n-1, n-1] = αs[1]
    R[n-1, n] = βs[1]
    R[n,n] = ĉ₀

    CTs, R
end


# find roots of P(p) without converting into standard basis
# using comrade_matrix
"""
    roots(p::P) where {P <: AbstractCCOP}

When `P` is a classical orthogonal family, finds the roots of `p` using the comrade matrix; otherwise falls back to finding the roots after conversion to the `Polynomial` type.

The eigenvalue of the comrade matrix are identified using an algorithm from

*Fast computation of eigenvalues of companion, comrade, and related matrices*
by Jared L. Aurentz, Raf Vandebril, and David S. Watkins
DOI [10.1007/s10543-013-0449-x](https://doi.org/10.1007/s10543-013-0449-x)
BIT Numer Math (2014) 54:7–30

The paper presents both a single- and double-shift version. The paper
has this caveat: "The double-shift code was less successful. It too is
very fast, but it is too unstable in its current form. We cannot
recommend it for accurate computation of roots of high-degree
polynomials. It could find use as an extremely fast method to get a
rough estimate of the spectrum"

"""
function Polynomials.roots(p::P) where {P <: AbstractCCOP}
    Cs, B = comrade_decomposition(P, coeffs(p))
    λs = AVW_eigvals(Cs, B)

    # hacky check on whether roots are good enough
    ns = p.(λs) ./ derivative(p).(λs)
    δ = maximum(abs,ns)
    if δ  >= sqrt(eps())
        rts = eigvals(Polynomials.companion(p))
        if maximum(abs, p.(λs)) > maximum(abs, p.(rts))
            return complex.(rts)
        end
    end
    return newton_refinement(λs, p)
end


## -----

## test root quality
## could presumably be much more efficient
## This test backwards stability
function a_posteriori_check(λs, p::P) where {P <: AbstractCCOP}
    ps = coeffs(p)
    n = length(ps) - 1
    v(λ) = [basis(P, i)(λ) for i ∈ n-1:-1:0]
    αₙ = 1/An(P,  n - 1)
    cₙ = ps[end-1]/ps[end] / αₙ
    p = P(ps)
    A = comrade_matrix(p)
    Aₒₒ = norm(A, Inf)
    [norm(αₙ * p(λ) / cₙ) /  (Aₒₒ * norm(v(λ), Inf)) for λ ∈ λs]
end
