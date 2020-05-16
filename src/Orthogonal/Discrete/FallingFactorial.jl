# represent the   basis x^(mbar) = x⋅(x-1)⋯(x-m+1) ∈ Πm

struct FallingFactorial{T <: Number} <: AbstractSpecialPolynomial{T}
    coeffs::Vector{T}
    var::Symbol
    function FallingFactorial{T}(coeffs::Vector{T}, var::Symbol) where {T <: Number}
        length(coeffs) == 0 && return new{T}(zeros(T, 1), var)
        last_nz = findlast(!iszero, coeffs)
        last = max(1, last_nz === nothing ? 0 : last_nz)
        return new{T}(coeffs[1:last], var)
    end
end

Polynomials.@register FallingFactorial
export  FallingFactorial

# show as 1⋅(x)₀ + 2⋅(x)₁ + 3⋅(x)₂
function Polynomials.showterm(io::IO, ::Type{P}, pj::T, var, j, first::Bool, mimetype) where {N, T, P <: FallingFactorial}
    iszero(pj) && return false
    !first &&  print(io, " ")
    print(io, Polynomials.hasneg(T)  && Polynomials.isneg(pj) ? "- " :  (!first ? "+ " : ""))
    print(io, "$(abs(pj))⋅")
    print(io,"($(var))")
    unicode_subscript(io, j)
    return true
end


function  (p::FallingFactorial)(x)
    d = degree(p)
    d <= 0 &&  return  p[0]*one(x)

    xⁱ̲ =  one(x)
    tot =  p[0]*xⁱ̲
    for i in  1:d
        xⁱ̲ *= x-(i-1)
        tot = muladd(xⁱ̲, p[i], tot)
    end
    
    tot
end

Base.:*(p::FallingFactorial, q::FallingFactorial) = convert(FallingFactorial,  convert(Polynomial,p)*convert(Polynomial,q))


## Connect FallingFactorial with AbstractCDOP
function Base.convert(::Type{P}, q::Q) where {
    P <: FallingFactorial,
    Q <: AbstractCDOP}
    _convert_ccop(P,q)
end
function Base.convert(::Type{P}, q::Q) where {
    P <: AbstractCDOP,
    Q <: FallingFactorial}
    _convert_ccop(P,q)
end

##  Koepf and Schmersa Thm 6 connection for P, Q=x^n̲
function connection_m(::Type{P},::Type{Q},m,n) where {
    P <: AbstractCDOP,
    Q <: FallingFactorial}

    a,b,c,d,e = abcde(P)

    c₀ = (a*n + a*m - a + d) * (n - m)

    c₁ = (m+1) * (a*n^2 - 2a*m^2 - a*n - a*m  + n*d - 2d*m  - b*m - d - e)

    c₂ = -(m+1) * (m+2) * (a*m^2 +  2a*m + d*m + b*m + a + d + b + c + e)

    (c₀,c₁,c₂)
end


##  Koepf and Schmersa Thm 7, p25 connection for P=x^n̲, Q
function connection_m(::Type{P},::Type{Q},m,n) where {
    P <: FallingFactorial,
    Q <: AbstractCDOP}

    ā,b̄,c̄,d̄,ē = abcde(Q)

    c₀ =  (2m*ā + ā + d̄)  * (2m*ā +  3ā + d̄)*(2m*ā + 2ā + d̄)^2 * (2m*ā + d̄) * (n-m)

    c₁ =  (2m*ā + ā + d̄)  * (2m*ā +  3ā + d̄)*(2m*ā + 2ā + d̄)  * (m+1)
    c₁a = 2m^2*n*ā^2  - 2m^2*ā^2
    c₁a += m^2*ā*d̄ + 2m^2*ā*b̄ +2m*n*ā^2  +  2m*n*ā*d̄ - 2m*ā^2 - m*ā*d̄ + 2m*ā*b̄ + m*d̄^2
    c₁a += 2m*d̄*b̄ + n*ā*d̄ +  2n*ā*ē - n*d̄*b̄ -  ā*d̄ + d̄*b̄  + d̄*ē
    c₁ *= c₁a

    c₂ = (m+1) * (2m*ā + d̄)
    c₂a = m^4*ā^3 + 4m^3*ā^3 + 2m^3*ā^2*d̄ + 6m^2*ā^3 + 6m*ā^2*d̄
    c₂a += 8m*ā^2*c̄ + 4m*ā^2*ē + 2m*ā*d̄^2 - 2m*ā*d̄*b̄ + 4m*ā*d̄*c̄ + 2m*ā*d̄*ē
    c₂a += -2m*ā*b̄^2 - m*d̄^2*b̄ - m*d̄*b̄^2 + ā^3 + 2ā^2*d̄ +  4ā^2*c̄ + 2ā^2*ē +  ā*d̄^2
    c₂a += -ā*d̄*b̄  + 4ā*d̄*c̄ + 2*ā*d̄*ē - ā*b̄^2 + ā*ē^2 - d̄^2*b̄ + d̄^2*c̄
    c₂a += -d̄*b̄^2 - d̄*b̄*ē
    c₂ *= c₂a * (m+2) *  (m*ā + n*ā + ā + d̄)

    (c₀,c₁,c₂)
end



function Base.convert(::Type{P}, q::Q) where {
    P <: Polynomials.StandardBasisPolynomial,
    Q <: FallingFactorial}
    q(variable(P))
end
function Base.convert(::Type{P}, q::Q) where {
    P <: FallingFactorial,
    Q <: Polynomials.StandardBasisPolynomial}
    connection(P,q)
end


# stirling of 2nd kind
@memoize function sterling²(n::Int, k::Int)
    (iszero(n) && iszero(k)) &&  return 1
    (iszero(n) || iszero(k)) && return 0

    k*sterling²(n-1,k) + sterling²(n-1, k-1)
end

# xʲ = ∑ᵢ sterling²(j,i) (x)ᵢ [wikipedia](https://en.wikipedia.org/wiki/Falling_and_rising_factorials)
function Base.iterate(o::Connection{P, Q}, state=nothing) where {P<:FallingFactorial, Q<:Polynomials.StandardBasisPolynomial}
    n,k = o.n, o.k
    if state  == nothing
        j = k
    else
        j = state
        j += 1
    end
    j > n && return nothing
    val = sterling²(j,k)
    (j,val), j
end


    


