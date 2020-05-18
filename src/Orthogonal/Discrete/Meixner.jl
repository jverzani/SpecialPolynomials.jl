@register2  Meixner AbstractCDOP2
export  Meixner

"""
    Meixner{γ,μ}


"""
Meixner

basis_symbol(::Type{<:Meixner{γ,μ}}) where {γ,μ} = "m⁽ᵞᵝ⁾" 
Polynomials.domain(::Type{<:Meixner{γ, μ}}) where {γ, μ} = Polynomials.Interval(-Inf, Inf)
abcde(::Type{<:Meixner{γ,μ}})  where {γ,μ} = NamedTuple{(:a,:b,:c,:d,:e)}((0,1,0,μ-1,γ*μ))

function kn(P::Type{<:Meixner{γ,μ}}, n::Int) where {γ,μ}
    (1 - 1/μ)^n
end

function k1k0(P::Type{<:Meixner{γ,μ}}, n::Int) where {γ,μ}
    (1 - 1/μ)
end
function k1k_1(P::Type{<:Meixner{γ,μ}}, n::Int) where {γ,μ}
    (1 - 1/μ)^(n>0 ? 2 : 1)
end

function  classical_hypergeometric(P::Type{<:Meixner{γ,μ}}, n::Int,   x) where {γ,μ}
    Pochhammer(γ,n) * pFq((-n, -x), γ, 1 - 1/μ)
end

