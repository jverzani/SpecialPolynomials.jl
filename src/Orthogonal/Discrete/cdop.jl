# generic classical discrete orthogonal polynomial, CDOP
abstract type AbstractCDOP{T,N} <: AbstractContinuousOrthogonalPolynomial{T} end

## Structural Equations
function  An(P::Type{<:Abstra    
    num = (2b*n*(a*n+d-a)-e*(-d+2a))
    den = (d+2a*n) * (d-2a+2a*n)

    iszero(den) && return Bn(P, Val(n))  
    
    val = one(S) * num / den
    
    val
end

function Cn(P::Type{<:AbstractCDOP}, n::Int)
    a,b,c,d,e = abcde(P)
    _Cn(P, a,b,c,d,e,n) *  k1k_1(P, n)
end

function _Cn(P::Type{<:AbstractCDOP}, a,b,c,d,e, n::Int)

    S = eltype(P)
    
    numa = (a*n+d-2a) * n * (4c*a-b^2) + 4a^2*c -a*b^2 + a*e^2 - 4*a*c*d
    numa += d*b^2 - b*e*d + d^2*c
    num = -numa * (a*n + d - 2a) * n
    den = (d - 2a + 2a*n)^2 * (2a*n - 3a + d) * (2a*n - a + d)

    iszero(den) && return Cn(P, Val(n))
    
    val = one(S) * num  / den
    
        # oops, this is the discrete case
#        val *= -((n-1)*(d+a*n-a)*(a*n*d-d*b-a*d+a^2*n^2-2a^2*n+4c*a+a^2+2e*a-b^2)-d*b*e+d^2*c+a*e^2)
#        val *= (a*n+d-2a)*n
#        val /=  (d-a+2a*n)*(d+2a*n-3a)*(2a*n-2a+d)^2
#        val *= k1k_1(P, n, S)

    val
end

