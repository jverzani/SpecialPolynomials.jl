@registerN Hahn AbstractCDOP3 α β 𝐍
export Hahn

"""
    Hahn

References: [Koekoek and Swarttouw §1.5](https://arxiv.org/pdf/math/9602214.pdf)

!!!  Note
      This is not  correct; tests are broken
"""
Hahn

abcde(::Type{<:Hahn{α,β,𝐍}})  where {α,β,𝐍} = NamedTuple{(:a,:b,:c,:d,:e)}((1, -(β+𝐍+1), 0,  (α+β+2), -𝐍*(α+1)))

basis_symbol(::Type{<:Hahn{α,β,𝐍}}) where {α,β,𝐍} = "h⁽ᵅᵝ⁾" * "₍" * sprint(io -> unicode_subscript(io, 𝐍)) * "₎"
Polynomials.domain(::Type{<:Hahn{α, β, 𝐍}}) where {α, β, 𝐍} = Polynomials.Interval(0,𝐍)

#  p22 (α+β+2n,  n)
function kn(P::Type{<:Hahn{α,β,𝐍}}, n::Int) where {α,β,𝐍}
    generalized_binomial(one(eltype(P))*(α  +  β +  2n),   n)
end

#function k1k0(P::Type{<:Hahn{α,β}}, n::Int) where {α,β}
#end
#function k1k_1(P::Type{<:Hahn{α,β}}, n::Int) where {α,β}
#end


function classical_hypergeometric(::Type{<:Hahn{α,β,𝐍}}, n::Int, x)  where {α,β,𝐍}
    (-1)^n / gamma(1+n) * Pochhammer(β+1,n) * Pochhammer(𝐍-n,n)  *  pFq((-n, -x, n+1+α+β),  (β+1,  1-𝐍),  1)
end

##################################################

@registerN HahnQ AbstractCDOP3  α β 𝐍
export HahnQ

"""
"""
HahnQ

abcde(::Type{<:HahnQ{α,β,𝐍}})  where {α,β,𝐍} = NamedTuple{(:a,:b,:c,:d,:e)}((1, -(β+𝐍+1), 0,  α+β+2, -𝐍*(α+1)))

basis_symbol(::Type{<:HahnQ{α,β}}) where {α,β} = "Q⁽ᵅᵝ⁾"
Polynomials.domain(::Type{<:HahnQ{α, β}}) where {α, β} = Polynomials.Interval(-Inf, Inf)


#  p22 table  (α+β+n+1)_n/((-N)_n * (α+1)_n )
function kn(P::Type{<:HahnQ{α,β,𝐍}}, n::Int) where {α,β,𝐍}

    as = [(α+β+n+1),-𝐍, (α+1)]
    
    val = one(eltype(P))
    val *= as[1]/as[2]/as[3]
    for  i in  1:n-1
        as .+= 1
        val *=  as[1]/as[2]/as[3]
    end
    
    val
end

#function k1k0(P::Type{<:HahnQ{α,β}}, n::Int) where {α,β} end function
#k1k_1(P::Type{<:HahnQ{α,β}}, n::Int) where {α,β} end

function classical_hypergeometric(::Type{<:HahnQ{α,β,𝐍}}, n::Int, x)  where {α,β,𝐍}
    pFq((-n, -x, n+1+α+β),  (α+1, -𝐍),  1)
end
