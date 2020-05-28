@registerN  Charlier AbstractCDOP1 μ
export  Charlier


"""
    Charlier

References: [Koekoek and Swarttouw §1.12](https://arxiv.org/pdf/math/9602214.pdf)
"""
Charlier

basis_symbol(::Type{<:Charlier{μ}}) where {μ} = "cᵘ"
Polynomials.domain(::Type{<:Charlier{μ}}) where {μ} = Polynomials.Interval(0, Inf)
abcde(::Type{<:Charlier{μ}})  where {μ} = NamedTuple{(:a,:b,:c,:d,:e)}((0,1,0,-1,μ))

function kn(P::Type{<:Charlier{μ}}, n::Int)  where {μ}
    (-one(eltype(P))/μ)^n
end

function k1k0(P::Type{<:Charlier{μ}}, n::Int)  where {μ}
    -one(eltype(P))/μ
end
function k1k_1(P::Type{<:Charlier{μ}}, n::Int)  where {μ}
    m = n ==  0 ?  1 : 2
    (-one(eltype(P))/μ)^m
end

function classical_hypergeometric(P::Type{<:Charlier{μ}}, n::Int, x) where {μ}
    pFq((-n,-x),(),-1/μ)
end

    
    
