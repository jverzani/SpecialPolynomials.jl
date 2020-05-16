@register3 Hahn AbstractCDOP3
export Hahn

"""
"""
Hahn

basis_symbol(::Type{<:Hahn{α,β,𝐍}}) where {α,β,𝐍} = "h⁽ᵅᵝ⁾" * "₍" * sprint(io -> unicode_subscript(io, 𝐍)) * "₎"
Polynomials.domain(::Type{<:Hahn{α, β}}) where {α, β} = Polynomials.Interval(-Inf, Inf)
abcde(::Type{<:Hahn{α,β,𝐍}})  where {α,β,𝐍} = NamedTuple{(:a,:b,:c,:d,:e)}((1, -(β+𝐍+1), 0,  α+β+2, -𝐍*(α+1)))

function kn(P::Type{<:Hahn{α,β}}, n::Int) where {α,β}
    generalized_binomial(α  +  β +  2n,   n)
end

#function k1k0(P::Type{<:Hahn{α,β}}, n::Int) where {α,β}
#end
#function k1k_1(P::Type{<:Hahn{α,β}}, n::Int) where {α,β}
#end


##################################################

@register3 HahnQ AbstractCDOP3
export HahnQ

"""
"""
HahnQ

basis_symbol(::Type{<:HahnQ{α,β}}) where {α,β} = "Q⁽ᵅᵝ⁾"
Polynomials.domain(::Type{<:HahnQ{α, β}}) where {α, β} = Polynomials.Interval(-Inf, Inf)
abcde(::Type{<:HahnQ{α,β,𝐍}})  where {α,β,𝐍} = NamedTuple{(:a,:b,:c,:d,:e)}((1, -(β+𝐍+1), 0,  α+β+2, -𝐍*(α+1)))
#  p22
function kn(P::Type{<:HahnQ{α,β,𝐍}}, n::Int) where {α,β,𝐍}

    as = [(α+β+1),-𝐍, (α+1)]
    
    val =one(P)  * as[1]/as[2]/as[3]
    for  i in  1:n
        as .+= 1
        val *  as[1]/as[2]/as[3]
    end
    val
end

#function k1k0(P::Type{<:HahnQ{α,β}}, n::Int) where {α,β} end function
#k1k_1(P::Type{<:HahnQ{α,β}}, n::Int) where {α,β} end
