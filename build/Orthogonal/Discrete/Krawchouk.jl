@registerN  Krawchouk AbstractCDOP2 p 𝐍
export  Krawchouk


"""
     Krawchouk{p,𝐍}

Also spelled  Krawtchouk,  Kravhcuk,….

References: [Koekoek and Swarttouw §1.10](https://arxiv.org/pdf/math/9602214.pdf);  see  also  [Coleman](https://arxiv.org/pdf/1101.1798.pdf) for a different  parameterization.
"""
Krawchouk

abcde(::Type{<:Krawchouk{p,𝐍}})  where {p,𝐍} = NamedTuple{(:a,:b,:c,:d,:e)}((0, p-1, 0, 1, -p*𝐍))

basis_symbol(::Type{<:Krawchouk{p,𝐍}}) where {p,𝐍} = "kᵖ" * "₍" * sprint(io -> unicode_subscript(io, 𝐍)) * "₎"
Polynomials.domain(::Type{<:Krawchouk{p,𝐍}}) where {p,𝐍} = Polynomials.Interval(0, 𝐍)
weight_function(::Type{<:Krawchouk{p,𝐍}})  where {p,𝐍} = x  -> generalized_binomial(𝐍,x)  * p^x * (1-p)^(𝐍-x)


function k1k0(P::Type{<:Krawchouk}, n::Int) 
    one(eltype(P))/(n+1)
end


function classical_hypergeometric(P::Type{<:Krawchouk{p,N}}, n::Int, x) where {p,N}
    monic = classical_hypergeometric(MonicKrawchouk{p,N},n,x)
    monic * kn(P,n)
end

##
## --------------------------------------------------
##

@registerN MonicKrawchouk AbstractCDOP2 p 𝐍
export MonicKrawchouk
ϟ(::Type{<:MonicKrawchouk{p,𝐍}}) where {p,𝐍} = Krawchouk{p,𝐍}
ϟ(::Type{<:MonicKrawchouk{p,𝐍,T}}) where {p,𝐍,T} = Krawchouk{p,𝐍,T}
@register_monic(MonicKrawchouk)

function classical_hypergeometric(P::Type{<:MonicKrawchouk{p,N}}, n::Int, x) where {p,N}
    Pochhammer( -N,n) * p^n * pFq((-n, -x), -N, 1/p)
end
