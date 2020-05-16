@register3 DiscreteChebyshev AbstractCDOP3
export DiscreteChebyshev

"""
"""
DiscreteChebyshev

basis_symbol(::Type{<:DiscreteChebyshev{α,β,𝐍}}) where {α,β,𝐍} = "K⁽ᵅᵝ⁾"
Polynomials.domain(::Type{<:DiscreteChebyshev{α, β}}) where {α, β} = Polynomials.Interval(-Inf, Inf)
abcde(::Type{<:DiscreteChebyshev{α,β,𝐍}})  where {α,β,𝐍} = NamedTuple{(:a,:b,:c,:d,:e)}((0,0,1,α,β))

function kn(P::Type{<:DiscreteChebyshev{α,β}}, n::Int) where {α,β}
    α^n
end

function k1k0(P::Type{<:DiscreteChebyshev{α,β}}, n::Int) where {α,β}
    α
end
function k1k_1(P::Type{<:DiscreteChebyshev{α,β}}, n::Int) where {α,β}
    m = n==0 ? 1 : 2
    α^2
end


