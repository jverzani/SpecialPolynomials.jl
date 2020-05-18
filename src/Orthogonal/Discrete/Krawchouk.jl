@register2  Krawchouk AbstractCDOP2
export  Krawchouk


"""
"""
Krawchouk

basis_symbol(::Type{<:Krawchouk{p,𝐍}}) where {p,𝐍} = "kᵖ" * "₍" * sprint(io -> unicode_subscript(io, 𝐍)) * "₎"
Polynomials.domain(::Type{<:Krawchouk{p}}) where {p} = Polynomials.Interval(-Inf, Inf)
abcde(::Type{<:Krawchouk{p,𝐍}})  where {p,𝐍} = NamedTuple{(:a,:b,:c,:d,:e)}((0, p-1, 0, 1, -p*𝐍))

function kn(P::Type{<:Krawchouk}, n::Int) 
    one(eltype(P))/gamma(1+n)
end

function k1k0(P::Type{<:Krawchouk}, n::Int) 
    one(eltype(P))/(n+1)
end
function k1k_1(P::Type{<:Krawchouk}, n::Int)
    one(eltype(P))/ (n*(n+1))
end

function classical_hypergeometric(P::Type{<:Krawchouk{p,N}}, n::Int, x) where {p,N}
    monic = classical_hypergeometric(MonicKrawchouk{p,N},n,x)
    monic * kn(P,n)
end

@register2 MonicKrawchouk AbstractCDOP2
export MonicKrawchouk
abcde(::Type{<:MonicKrawchouk{p,𝐍}})  where {p,𝐍} = NamedTuple{(:a,:b,:c,:d,:e)}((0, p-1, 0, 1, -p*𝐍))

function classical_hypergeometric(P::Type{<:MonicKrawchouk{p,N}}, n::Int, x) where {p,N}
    Pochhammer( -N,n) * p^n * pFq((-n, -x), -N, 1/p)
end
