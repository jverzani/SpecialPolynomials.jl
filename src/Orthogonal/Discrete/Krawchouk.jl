@register2  Krawchouk AbstractCDOP2
export  Krawchouk


"""
"""
Krawchouk

basis_symbol(::Type{<:Krawchouk{p,𝐍}}) where {p,𝐍} = "kᵖ" * "₍" * sprint(io -> unicode_subscript(io, 𝐍)) * "₎"
Polynomials.domain(::Type{<:Krawchouk{p}}) where {p} = Polynomials.Interval(-Inf, Inf)
abcde(::Type{<:Krawchouk{p,𝐍}})  where {p,𝐍} = NamedTuple{(:a,:b,:c,:d,:e)}((0, p-1, 0, 1, -p*𝐍))

function kn(P::Type{<:Krawchouk}, n::Int) 
    1/gamma(1+n)
end

function k1k0(P::Type{<:Krawchouk}, n::Int) 
    n+1
end
function k1k_1(P::Type{<:Krawchouk}, n::Int)
    n ==  0 ?  n+1 : (n+1)*(n+2)
end

