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
    T,S = eltype(one(Q)), eltype(p)  #
    R = typeof(one(promote_type(T,S))/1)

    as = zeros(R, 1+d)

    λⱼⱼ = one(R) #  kn(P,0,R)/kn(Q,0,R) = 1
    
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

    ⟒(Q){R}(as, p.var)
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

    ⟒(P)(cs, q.var)
    
end


## For this implementation, we suppose there is an is an iterator
## C_{n, k} which iterates
##
## (k, α(k,k)), (k+1, α(k+1, k)), ..., (n, α(n, k))
##
## defined by iterate(Connection{ϕ, ψ}(n, k)
##
## When ψ_n  = ∑_k^n C_{n,k} ϕ_k
##
## the iterator for (n,k) returns
## (i, val) ∈  { (k, α_{k,k}), (k+1, α_{k+1,k}),  ..., (n, α_{n,k}) }
##
## A sample would be to convert polys in Q to polys in P:
# function Base.iterate(o::Connection{P, Q}, state=nothing) where
#     {P <: Polynomial, Q <: Polynomial}
#
#     n,k = o.n, o.k
#     if state  == nothing
#         j = k
#     else
#         j = state
#         j += 1  # likely j += 2 if parity condition
#     end
#     j > n && return nothing
#     val = connection_α(P,Q,j,k)
#     (j,val), j
# end



##  Linearization...

#Linearization{P}(l,m,n) where {P} = Linearization{P, Val{:diagonal}}(l,m,n)

# function linearization_product_pq(p::P,  q::Q) where {P <: AbstractOrthogonalPolynomial, Q <: AbstractOrthogonalPolynomial}

#     ⟒(P) == ⟒(Q) || throw(ArgumentError("Base polynomial type must match"))
#     PP = promote_type(P,Q)

#     n,m = length(p)-1, length(q)-1 # n,m = degree(p), degree(q)
#     T,S = eltype(p), eltype(q)
#     R = eltype(one(T)*one(S)/1)
#     cs = zeros(R, n+m+1)

#     ## pn*pm = ∑ a(l,m,n) H_l
#     for n in eachindex(p)
#         p_n = p[n]
#         iszero(p_n) && continue
#         for m in eachindex(q)
#             c_nm = p_n * q[m]
#             for (k, val) in Linearization{P, Val{:pq}}(0, n,m)
#                 cs[1+k] = muladd(c_nm, val, cs[1+k])
#             end
#         end
#     end
    
#     ⟒(P)(cs, p.var)
# end

# Sparser product, implented by fixing k and finding combinations of (m,n) that produce a
# non-zero α(k,n,m) for P_k through the linearization formula P_n⋅P_m = ∑ α(k,n,m) P_k
#
# compute product using a linearization formula, as implemented in an `iterate` method for the type
#
# @inline function linearization_product(pn::P,  pm::Q) where {P <: AbstractOrthogonalPolynomial, Q <: AbstractOrthogonalPolynomial}


#     Polynomials.isconstant(pn) && return pm * pn[0]
#     Polynomials.isconstant(pm) && return pn * pm[0]
#     pn.var == pm.var || throw(ArgumentError("bases must match"))
    
#     ⟒(P) == ⟒(Q) || throw(ArgumentError("Base polynomial type must match"))
#     PP = promote_type(P,Q)

#     n,m = length(pn)-1, length(pm)-1 # n,m = degree(pn), degree(pm)
#     T,S = eltype(pn), eltype(pm)
#     R = eltype(one(T)*one(S)/1)
#     cs = zeros(R, n+m+1)
    
#     # pn*pm = ∑ a(l,m,n) H_l
#     # so ∑_l ∑_p ∑_q a(l,p,q) H_l
#     # can  use iterator to keep as simple sum
#     # `sum(pn[p] * pm[q] * val for (p,q,val) in Linearization{PP}(l,n,m))`
#     for l in 0:n+m
#         for (p,q,val) in Linearization{PP, R}(l,n,m)
#             cs[1+l] += pn[p] * pm[q] * val
#         end
#     end
#     ⟒(P)(cs, pn.var)
# end


