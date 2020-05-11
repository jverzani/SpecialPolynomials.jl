## Bessel
@register1 Bessel


"""
    Bessel{α}

Implements the [Bessel](https://dlmf.nist.gov/18.34) polynomials, introduced by [Krall and Frink](https://www.ams.org/journals/tran/1949-065-01/S0002-9947-1949-0028473-1/S0002-9947-1949-0028473-1.pdf) (with `b=2`). The  case `a=2` corresponds to the 
[Bessel](https://en.wikipedia.org/wiki/Bessel_polynomials) polynomials of Wikipedia. The Bessel  polynomials are not orthogonal over  a domain of the real  line, rather over an arbitray curve in the complex plane enclosing the  origin.  The weight  function is `ρ(x)=(2πi)^(-1)∑Γ(α)/Γ(α+n-1)(-β/x)^n`,   where `β=2`.

```jldoctest
julia> using Polynomials, SpecialPolynomials

julia> 𝐐 = Rational{Int}
Rational{Int64}

julia> x = variable(Polynomial{𝐐})
Polynomial(x)

julia> [basis(Bessel{1//2, 𝐐}, i)(x) for i in 0:5]

6-element Array{Polynomial{Rational{Int64}},1}:
 Polynomial(1//1)
 Polynomial(1//1 + 1//4*x)
 Polynomial(1//1 + 3//2*x + 15//16*x^2)
 Polynomial(1//1 + 15//4*x + 105//16*x^2 + 315//64*x^3)
 Polynomial(1//1 + 7//1*x + 189//8*x^2 + 693//16*x^3 + 9009//256*x^4)
 Polynomial(1//1 + 45//4*x + 495//8*x^2 + 6435//32*x^3 + 96525//256*x^4 + 328185//1024*x^5)
```

"""
Bessel
export Bessel
basis_symbol(::Type{<:Bessel{α}}) where {α} = "C^($α)"
abcde(::Type{<:Bessel{α}})  where {α} = NamedTuple{(:a,:b,:c,:d,:e)}((1,0,0,α, 2))

# From https://www.ams.org/journals/tran/1949-065-01/S0002-9947-1949-0028473-1/S0002-9947-1949-0028473-1.pdf we use
# kn = (n + α - 1)_n / 2^n not (n+α+1)_n/2^n
# Koepf suggests kn = (n + α + 1)_n/2^n
function kn(::Type{<:Bessel{α}}, n::Int) where {α}
    one(α)*Pochhammer(n+α-1,n)/2^n
end
function k1k0(::Type{<:Bessel{α}}, k, ::Type{S}=Float64) where {S,α}
    iszero(k) && return α/2
    val = one(S)
    val *=  (2k+α)*(2k+α-1)
    val /= (k+α-1)*2
    val
end
function k1k_1(P::Type{<:Bessel{α}}, k, ::Type{S}=Float64) where {S,α}
    @assert k > 0
    val = one(S)
    if k == 1
        val *= (α + 1)*(α + 2) # cancels  factor (α-1)(α-1)
        val /= 4
    else
        val *= (2k+α-3) * (2k+α-2) * (2k+α-1) * (2k+α)
        val /= 4 * (k+α-2) * (k+α-1)
    end
    return val
end

## Overrides
Bn(::Type{<:Bessel{2}}, ::Val{0}, ::Type{S}) where  {S}  =  zero(S)
Cn(::Type{<:Bessel{α}}, ::Val{1}, ::Type{S}) where {α,S} = -α^2*(1+α)/4
