## General weight function
## Wheeler algorithm, https://authors.library.caltech.edu/86868/1/1.4822929.pdf
## algorithm not considered, just text

struct Weight{P, F, T <: Number} <: OrthogonalPolynomial{T}
    pis::P
    W::F
    coeffs::Vector{T}
    var::Symbol
#    function Weight(pis::Type{P}, w::F, coeffs::AbstractVector{T}, var::Symbol=:x) where {P<:OrthogonalPolynomial, F, T <: Number}
#        length(coeffs) == 0 && return new{T}(zeros(T, 1), var)
#        last_nz = findlast(!iszero, coeffs)
#        last = max(1, last_nz === nothing ? 0 : last_nz)
#        return new{P,F,T}(pis, w, coeffs[1:last], var)
#    end
end


export Weight

## conversion functions ....

basis_symbol(::Type{<:Weight}) = "W"

Polynomials.domain(w::Weight) = domain(w.pis)
weight_function(W::Type{<:Weight}) = W.W

## v_j = \int pi_j W(x) dx
## should cache thesex
function modified_moment(w::Weight,  j)
    dom = domain(w.pis)
    p = Polynomials.basis(w.pis, j)
    @show convert(Polynomial, p)
    fn = x -> p(x)  * weight_function(w.pis)(x)
    a,err =  quadgk(fn, first(dom), last(dom))
    a
end

Polynomials.variable(::Type{P}, var::Polynomials.SymbolLike=:x) where {P <: Weight} = P([0, 1], var)

## If the moment, v_k =\int_a^b x^k W dx are known, then Golub and
## Welsh
## (https://web.stanford.edu/class/cme335/S0025-5718-69-99647-1.pdf)
## give an algorithm tot find An, Bn, Cn.
## However, the moments arer numerically unstable.
## The Wheeler approach, detailed in Press and Tuloksky (https://authors.library.caltech.edu/86868/1/1.4822929.pdf)
## uses a different formula.

##  Use Wheeler here

#@memoize
alpha(p::P, n) where {P <: OrthogonalPolynomial} = -Bn(p,n)/ An(P,n)
alpha(P::Type{<:OrthogonalPolynomial}, n) = -Bn(P,n)/An(P,n)

#@memoize
function beta(p::P, n) where {P <: OrthogonalPolynomial}
    iszero(n) &&  return 0/1
    -Cn(P,n)/An(P,n)/An(P, n-1)
end

function beta(P::Type{<:OrthogonalPolynomial}, n)
    iszero(n) &&  return 0/1
    -Cn(P,n)/An(P,n)/An(P, n-1)
end

# <p_k, pi_l>
@memoize function sigma(w::W, k, l) where {W <: Weight}
    P = w.pis

    if k == -1
        return 0
    elseif  iszero(k)
        return modified_moment(w, l)
    elseif iszero(l)
        # pi_l = 1
        dom = domain(w.pis)
        zs = zeros(Int, k+1); zs[end]=1
        pk = Weight(w.pis, w.W, zs, w.var)
        fn = x -> pk(x) * weight_function(w.pis)(x)
        a,err =  quadgk(fn, first(dom), last(dom))
        return a
    end
    sigma(w, k-1, l-1) - (An(w, k-1) - alpha(P, l)) * sigma(w, k-1, l) - Bn(w,k-1)*sigma(w, k-2, l) + beta(P,l)*sigma(w, k-1, l-1)
end


##
An(::Type{Weight{P,F,T}}, n) where {P,F,T} = 1.0

@memoize function Bn(w::Weight, k)
    P = w.pis
    if k == 0
        v0, v1 =  modified_moment(w, 0), modified_moment(w,1)
        alpha(P,0) + v1 / v0
    else
        alpha(P, k) - sigma(w, k-1,k)/sigma(w, k-1, k-1) + sigma(w, k, k+1)/sigma(w, k,k)
    end
end

@memoize function Cn(w::Weight, k)
    iszero(k) && return 0/1
    sigma(w, k, k)/sigma(w, k-1, k-1)
end


(ch::Weight)(x::S) where {S} = orthogonal_polyval(ch, x)
