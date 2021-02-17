## Connection and Linearization

## A connection between polynomial systems `P` and `Q` is a function `C(m,n)` satisfying
## p_n = ∑C(m,n) q_m or  pᵢ = ∑ Cʲᵢqⱼ (if unicode is helpful)
##
## From this a polynomial in P = ∑ aᵢ⋅pᵢ can be represented in a Q through
##
##  ∑ᵢ aᵢ⋅pᵢ = ∑ᵢ aᵢ  ∑ⱼ Cʲᵢqⱼ  = ∑ᵢ (∑ⱼ (aᵢ⋅Cʲᵢ)) qⱼ or ∑ⱼ (∑ᵢ aᵢCʲᵢ) qⱼ
##
## The choice of which depends on how the connection coefficients are presented.
##
## A connection between P and Q defines a `convert` method `Base.convert(Q,  p::P)`
## Conversion can also be achieved through polynomial evaluation. For example
## if `x=variable(Q)` then `Base.convert(Q,p) = p(x)` will also convert `p` into
## the basis Q. This is helpful when Q is a standard basis polynomial.
##
## However, if Q is not a standard basis polynomial, it isn't helpful unless
## one can multiply monomomials within Q. For standard basis polynomials this
## is known, as  xⁱ⋅xʲ = xⁱ⁺ʲ.
##
## For general Q though, one needs a linearization: a function α with
##
##  qᵢqⱼ = ∑_k α(k,i,j) q_k
##
## with a given α, multiplication can be defined. The linearization can be defined
## through composition: If P is the standard basis, connect Q to P, multiply ini P, then connect
## to Q.
##
## In short, if there is a linearization then multiplication can be defined and conversion defined through
## polynomial evaluation
## If there is a connection defined then conversion is defined and multiplication can be defined through composition
##


## For `AbstractCCOP` polynomials, Koepf and Schmersau present recursive formulas for the connection coefficients
## between Q and the standard basis, the standard and Q, and between P and Q when σ and σ̂ are the same. (This is the case when a polynomial system is parameterized, like Laguerre, Gegenbauer, Bessel, and  Jacobi.)

ConvertibleTypes = Union{AbstractCOP, Polynomials.StandardBasisPolynomial,  FallingFactorial}
kn(::Type{P},  n) where{P <: Union{Polynomials.StandardBasisPolynomial,  FallingFactorial}} = one(eltype(one(P)))
k1k0(::Type{P},  n) where{P <: Union{Polynomials.StandardBasisPolynomial,  FallingFactorial}} = one(eltype(one(P)))
knk_1(::Type{P},  n) where{P <: Union{Polynomials.StandardBasisPolynomial,  FallingFactorial}} = one(eltype(one(P))) 
                         


## If  the connection  coefficients for  (P,Q) satisfy c₀⋅ C̃ⁱⱼ + c₁⋅C̃ⁱ⁺¹ⱼ  + c₂⋅*C̃ⁱ⁺²ⱼ = 0 for some computable
##  values c₀,c₁,c₂ then this  will implement  `convert(Q,p::P)`
function _convert_cop(::Type{Q}, p::P) where {P <: ConvertibleTypes,
                                               Q <: ConvertibleTypes}

    d = degree(p)
    T,S = eltype(Q), eltype(p)  #
    R = typeof(one(promote_type(T,S))/1)

    as = zeros(R, 1+d)

    λⱼⱼ = one(R)  * k0(P) / k0(Q)
    
    for j in  0:d

        λ = λⱼⱼ 

        λⱼⱼ *= k1k0(P,j) / k1k0(Q,j)  

        pⱼ = p[j]
        iszero(pⱼ) && continue

        C̃ʲⱼ = one(R)  

        as[1+j] += pⱼ * (λ * C̃ʲⱼ)

        C̃ⁱ⁺²ⱼ, C̃ⁱ⁺¹ⱼ = zero(R), C̃ʲⱼ  # for recursion starting  with C̃ʲ⁺¹ⱼ, C̃ʲⱼ
        for i in j-1:-1:0

            λ *= k1k0(Q,i)

            # c₀⋅ C̃ⁱⱼ + c₁⋅C̃ⁱ⁺¹ⱼ  + c₂⋅*C̃ⁱ⁺²ⱼ) = 0
            c₀,c₁,c₂ = connection_m(P, Q, i, j)
            C̃ⁱⱼ = iszero(c₀) ? zero(R) : -one(R)*(c₁*C̃ⁱ⁺¹ⱼ + c₂*C̃ⁱ⁺²ⱼ) / c₀
            as[1+i] += pⱼ * (λ * C̃ⁱⱼ)

            C̃ⁱ⁺²ⱼ, C̃ⁱ⁺¹ⱼ = C̃ⁱ⁺¹ⱼ, C̃ⁱⱼ
        end
    end
    X = Polynomials.indeterminate(Q,p)
    ⟒(Q){R,X}(as)
end

## Use FastTransform for conversion, as  possible
## FastTransforms.kind2string.(0:15)
# Legendre<--Chebyshev
function _convert(::Type{Q}, p::P) where {
    T <: AbstractFloat,
    Q <: Union{Legendre, OrthonormalLegendre},
    P <: Union{Chebyshev{T}, OrthonormalChebyshev{T}}}
    ps = coeffs(p)
    Q(cheb2leg(ps, normcheb=isorthonormal(P), normleg=isorthonormal(Q)))
end

# Chebyshev<--Legendre
function _convert(::Type{Q}, p::P) where {
    T <: AbstractFloat,
    Q <: Union{Chebyshev, OrthonormalChebyshev},
    P <: Union{Legendre{T}, OrthonormalLegendre{T}}
}
    ps =  coeffs(p)
    Q(leg2cheb(ps, normleg=isorthonormal(P), normcheb=isorthonormal(Q)) )
end

# ultraspherical<--ultraspherical
function _convert(::Type{Q}, p::P) where {
    β, α,  T <: AbstractFloat,
    Q <: Union{Gegenbauer{β}, OrthonormalGegenbauer{β}},
    P <: Union{Gegenbauer{α, T}, OrthonormalGegenbauer{α, T}}
}
    
    ps =  coeffs(p)
    Q( ultra2ultra(ps, α, β, norm1=isorthonormal(P), norm2=isorthonormal(Q)) )
       
end
       
# Jacobi<--Jacobi
function _convert(::Type{Q}, p::P) where {
    γ, δ, α, β,  T <: AbstractFloat,
    Q <: Union{Jacobi{γ, δ}, OrthonormalJacobi{γ, δ}},
    P <: Union{Jacobi{α, β, T}, OrthonormalJacobi{α, β, T}}
}
    
    ps =  coeffs(p)
    Q( jac2jac(ps, α, β, γ, δ,  norm1=isorthonormal(P), norm2=isorthonormal(Q)) )
    
end

# Laguerre<--Laguerre
function _convert(::Type{Q}, p::P) where {
    α, β,  T <: AbstractFloat,
    Q <: Union{Laguerre{β}, OrthonormalLaguerre{β}},
    P <: Union{Laguerre{α, T}, OrthonormalLaguerre{α, T}}
}

    ps =  coeffs(p)
    Q( lag2lag(ps, α, β,  norm1=isorthonormal(P), norm2=isorthonormal(Q)) )
    
end


# Jacobi<--ultraspherical
function _convert(::Type{Q}, p::P) where {
    γ, δ, α,  T <: AbstractFloat,
    Q <: Union{Jacobi{γ, δ}, OrthonormalJacobi{γ, δ}},
    P <: Union{Gegenbauer{α, T}, OrthonormalGegenbauer{α, T}}
}
    
    ps =  coeffs(p)
    Q( ultra2jac(ps, α, γ, δ,  normultra=isorthonormal(P), normjac=isorthonormal(Q)) )
    
end

# ultraspherical<--Jacobi
function _convert(::Type{Q}, p::P) where {
    α, γ, δ,  T <: AbstractFloat,
    Q <: Union{Gegenbauer{α, T}, OrthonormalGegenbauer{α, T}},
    P <: Union{Jacobi{γ, δ}, OrthonormalJacobi{γ, δ}}
}
    
    ps =  coeffs(p)
    Q( jac2ultra(ps,γ, δ, α, normjac=isorthonormal(P), normultra=isorthonormal(Q)) )
    
end

# Jacobi<--Chebyshev
function _convert(::Type{Q}, p::P) where {
    γ, δ,  T <: AbstractFloat,
    Q <: Union{Jacobi{γ, δ}, OrthonormalJacobi{γ, δ}},
    P <: Union{Chebyshev{T}, OrthonormalChebyshev{T}}
}
    
    ps =  coeffs(p)
    Q( cheb2jac(ps,γ, δ, normcheb=isorthonormal(P), normjac=isorthonormal(Q)) )
    
end

# Chebyshev<--Jacobi
function _convert(::Type{Q}, p::P) where {
    γ, δ,  T <: AbstractFloat,
    Q <: Union{Chebyshev, OrthonormalChebyshev},
    P <: Union{Jacobi{γ, δ,T}, OrthonormalJacobi{γ, δ,T}}
}
    
    ps =  coeffs(p)
    Q( jac2cheb(ps,γ, δ, normjac = isorthonormal(P), normcheb=isorthonormal(Q)) )
    
end

# ultraspherical<--Chebyshev
function _convert(::Type{Q}, p::P) where {
    α,  T <: AbstractFloat,
    Q <: Union{Gegenbauer{α}, OrthonormalGegenbauer{α}},
    P <: Union{Chebyshev{T}, OrthonormalChebyshev{T}}
}

    ps =  coeffs(p)
    Q( cheb2ultra(ps,α, normcheb=isorthonormal(P), normultra=isorthonormal(Q)) )
    
end

# Chebyshev<--ultraspherical
function _convert(::Type{Q}, p::P) where {
    α,  T <: AbstractFloat,
    Q <: Union{Chebyshev, OrthonormalChebyshev},
    P <: Union{Gegenbauer{α,T}, OrthonormalGegenbauer{α,T}}
}
    ps =  coeffs(p)
    Q( ultra2cheb(ps,α, normultra=isorthonormal(P), normcheb=isorthonormal(Q)) )
    
end


##  Koepf and Schmersa Thm 2, p11 connection for P to Q when σ = σ̂
function connection_m(::Type{P},::Type{Q},m,n) where {
    P <: AbstractCCOP,
    Q <: AbstractCCOP}

    a,b,c,d,e = abcde(P)
    ā,b̄,c̄,d̄,ē = abcde(Q)

    (a == ā && b == b̄ && c == c̄) || throw(ArgumentError("This connection assume σ = σ̄"))
    a²,b²,c²,d²,d̄²,e²,ē²,m²,n² = a^2, b^2, c^2, d^2, d̄^2, e^2, ē^2, m^2, n^2
    
    c0 = -(m-n) ⋅ (a⋅m + d - a + a⋅n) ⋅ (d̄ + 2a⋅m) ⋅ (d̄ + a + 2a⋅m) ⋅ (d̄ + 3a + 2a⋅m)
    c0 *= (d̄ + 2a⋅m  + 2a)^2

    c1 = -d⋅b⋅n⋅d̄ + 2d⋅a⋅m²⋅b + d⋅b⋅d̄ + 2d⋅a⋅m⋅b + 2d⋅ē⋅n⋅a
    c1 += d⋅d̄⋅ē + 2d⋅d̄⋅b⋅m - m⋅b⋅d̄² - e⋅d̄² - 4a²⋅m²⋅e - m²⋅a⋅b⋅d̄ + b⋅n⋅d̄⋅a - 2e⋅d̄⋅a
    c1 += -4a²⋅m⋅e - 4e⋅d̄⋅a⋅m + 2m²⋅a²⋅ē + 2ē⋅a²⋅n² - 2ē⋅a²⋅n - m⋅a⋅b⋅d̄ + 2m⋅d̄⋅ē⋅a
    c1 += 2m⋅ē⋅a² - b⋅n²⋅d̄⋅a
    c1 *= (d̄ + 2a⋅m + 2a) ⋅ (m + 1) ⋅  (d̄ + a + 2a⋅m) ⋅ (d̄ + 3a + 2a⋅m)

    c2 =  -(d̄ + 2a⋅m) ⋅ (m + 1) ⋅ (-a⋅m - 2a + a⋅n - d̄ + d) ⋅ (a⋅m + a⋅n + a + d̄)
    c2a = a⋅b²⋅m² - 4a²⋅m²⋅c - 8a²⋅m⋅c + 2a⋅m⋅b² - 4a⋅d̄⋅m⋅c + m⋅b²⋅d̄ - 4a⋅d̄⋅c - a⋅ē² + a⋅b² - c⋅d̄²
    c2a += b⋅ē⋅d̄ - 4⋅a²⋅c + b²⋅d̄
    c2 *= c2a ⋅ (m + 2)

    
    return (c0, c1, c2)

end

##  Koepf and Schmersa Thm 4, p14 connection for P, Q=x^n
function connection_m(::Type{P},::Type{Q},m,n) where {
    P <: AbstractCCOP,
    Q <: Polynomials.StandardBasisPolynomial}

    a,b,c,d,e = abcde(P)

    c0 = (m-n) ⋅ (a⋅n + d - a + a⋅m)

    c1 = (m + 1) ⋅ (b⋅m + e)

    c2 = c ⋅ (m + 1) ⋅ (m + 2)
                  
    return (c0, c1, c2)

end

##  Koepf and Schmersa Thm 5, p16 connection for P=x^n, Q
function connection_m(::Type{P},::Type{Q},m,n) where {
    P <: Polynomials.StandardBasisPolynomial,
    Q <: AbstractCCOP}

    ā,b̄,c̄,d̄,ē = abcde(Q)
    ā²,b̄²,c̄²,d̄²,ē²,m²,n² = ā^2, b̄^2, c̄^2, d̄^2, ē^2, m^2, n^2
    
    c0 =  (n - m) ⋅ (d̄ + 2ā⋅m) ⋅ (d̄ + 3ā + 2ā⋅m) ⋅ (d̄ + ā + 2ā⋅m) ⋅ (d̄ + 2ā⋅m + 2ā)^2

    c1 =  (d̄⋅ē + b̄⋅d̄ +2d̄⋅b̄⋅m + 2ā⋅m²⋅b̄ + 2ā⋅m⋅b̄ +2ē⋅ā⋅n - d̄⋅b̄⋅n) ⋅ (d̄ + 2ā⋅m + 2ā)
    c1 *= (m+1) ⋅ (d̄ + 3ā + 2ā⋅m) ⋅  (d̄ + ā + 2ā⋅m)

    c2 = -(m+2)
    c2a = -4ā²⋅c̄⋅m²
    c2a +=  ā⋅b̄²⋅m² + 2ā⋅b̄²⋅m - 4ā⋅c̄⋅m⋅d̄ - 8ā²⋅c̄⋅m + m⋅b̄²⋅d̄ - ā⋅ē^2 - d̄²⋅c̄ + b̄⋅ē⋅d̄ -4ā²⋅c̄
    c2a += -4ā⋅c̄⋅d̄ + ā⋅b̄² +b̄²⋅d̄
    c2 *= c2a ⋅ (ā⋅m + ā⋅n + ā + d̄) ⋅ (m + 1) ⋅ (d̄ + 2ā⋅m)

    return (c0, c1, c2)

end

## With an explicit connectioni
function connection(::Type{P}, q::Q) where
    {P <: Polynomials.AbstractPolynomial,
     Q<:Polynomials.AbstractPolynomial}
    
    n = degree(q)
    T = eltype(q)
    T = eltype(one(eltype(one(P)))*one(q)/1)
    cs = zeros(T, n+1)

    for k in 0:n
        for (i,val) in Connection{P,Q}(n, k)
            cs[1+k] = muladd(q[i], val, cs[1+k])
        end
    end
    X = Polynomials.indeterminate(P,q)
    ⟒(P){T,X}(cs)
    
end

