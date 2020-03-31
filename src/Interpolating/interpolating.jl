abstract type AbstractInterpolatingPolynomial{T} <: AbstractSpecialPolynomial{T} end


function Polynomials.showterm(io::IO, ::Type{P}, pj::T, var, j, first::Bool, mimetype) where {N, T, P <: AbstractInterpolatingPolynomial}
    iszero(pj) && return false
    !first &&  print(io, " ")
    print(io, Polynomials.hasneg(T)  && Polynomials.isneg(pj) ? "- " :  (!first ? "+ " : ""))
    print(io, "$(abs(pj))⋅$(basis_symbol(P))_$(j)($var)")
    return true
end

basis_symbol(::Type{<:AbstractInterpolatingPolynomial}) = "P"

function Polynomials.variable(::Type{P}, var::Polynomials.SymbolLike=:x) where {P <: AbstractInterpolatingPolynomial}
    @warn "No `variable` function is defined for a type, just an instance"
    throw(MethodError())
end
