##
## --------------------------------------------------
##
# Macros to register POLY{T,X},  POLY{αs..., T,X}
#
# We use `Vector{T}` with N =  d+1 to store  a degree d  polynomial - no trailing zeros permitted.
##
macro register0(name, parent)
    poly = esc(name)
    parent_type = esc(parent)
    quote
        struct $poly{T,X} <: $parent_type{T,X}
            coeffs::Vector{T}
            function $poly{T,X}(::Val{false}, coeffs::Vector{T}) where {T,X}
                isempty(coeffs) && return new{T,X}(T[])
                if iszero(coeffs[end])
                    N = findlast(!iszero, coeffs)
                    N == nothing && return new{T,X}(T[])
                    resize!(coeffs, N)
                end
                new{T,X}(coeffs)
            end
            function $poly{T,X}(coeffs::AbstractVector) where {T,X}
                if Base.has_offset_axes(coeffs)
                    @warn "Ignoring offset indices of the coefficient vector"
                end
                $poly{T,X}(Val(false), collect(T, coeffs))
            end
            function $poly{T,X}(coeffs::Tuple) where {T,X}
                $poly{T,X}(Val(false), collect(T, coeffs))
            end
            function $poly{T}(coeffs::Vector, var::Polynomials.SymbolLike = Polynomials.Var(:x)) where {T}
                $poly{T, Symbol(var)}(coeffs)
            end
            function $poly{T}(coeffs::Tuple, var::Polynomials.SymbolLike = Polynomials.Var(:x)) where {T}
                $poly{T,Symbol(var)}(collect(T, coeffs))
            end
            function $poly(coeffs::Vector{T}, var::Polynomials.SymbolLike = Polynomials.Var(:x)) where {T}
                $poly{T}(coeffs, var)
            end
            function $poly(coeffs::Tuple{Vararg{T}}, var::Polynomials.SymbolLike = Polynomials.Var(:x)) where {T}
                $poly{T}(coeffs, var)
            end
        end

        Base.length(p::$poly{T,X}) where {T,X} = length(p.coeffs)
        Polynomials.evalpoly(x, p::$poly) = eval_cop(typeof(p), p.coeffs, x)
        Polynomials.@register $poly

        # work around N parameter in pomote(p,q) usage in defaults
        Base.:+(p::P, q::Q) where {T,X,P<:$poly{T,X},S,Q<:$poly{S,X}} = ⊕(P, p, q)
        Base.:*(p::P, q::Q) where {T,X,P<:$poly{T,X},S,Q<:$poly{S,X}} = ⊗(P, p, q)
        Base.divrem(p::$poly{T}, q::$poly{T}) where {T} = _divrem(p, q)
    end
end

# register macro for polynomial families with parameters
# would be nice to consolidate, but ...
macro registerN(name, parent, params...)
    poly = esc(name)
    parent_type = esc(parent)
    αs = tuple(esc.(params)...)
    quote
        struct $poly{$(αs...),T,X} <: $parent_type{$(αs...),T,X}
            coeffs::Vector{T}
            function $poly{$(αs...),T,X}(::Val{false}, coeffs::Vector{S}) where {$(αs...),T,X,S}
                N = findlast(!iszero, coeffs)
                N == nothing && return new{$(αs...),T,X}(T[])
                new{$(αs...),T,X}(T[coeffs[i] for i in firstindex(coeffs):N])
            end
            function $poly{$(αs...),T,X}(coeffs::AbstractVector) where {$(αs...),T,X}
                if Base.has_offset_axes(coeffs)
                    @warn "Ignoring offset indices of the coefficient vector"
                end
                $poly{$(αs...),T,X}(Val(false), collect(T, coeffs))
            end
            function $poly{$(αs...),T,X}(coeffs::Tuple) where {$(αs...),T,X}
                $poly{$(αs...),T,X}(Val(false), collect(T, coeffs))
            end
            function $poly{$(αs...),T}(coeffs::Vector,
                    var::Polynomials.SymbolLike = Polynomials.Var(:x)) where {$(αs...),T}
                $poly{$(αs...), T, Symbol(var)}(coeffs)
            end

            function $poly{$(αs...),T}(coeffs::Tuple,
                    var::Polynomials.SymbolLike = Polynomials.Var(:x)) where {$(αs...),T}
                $poly{$(αs...),T,Symbol(var)}(coeffs)
            end

            function $poly{$(αs...)}(coeffs::Vector{T},
                    var::Polynomials.SymbolLike = Polynomials.Var(:x)) where {$(αs...),T}
                $poly{$(αs...),T}(coeffs, var)
            end

            function $poly{$(αs...)}(coeffs::Tuple{Vararg{T}},
                    var::Polynomials.SymbolLike = Polynomials.Var(:x)) where {$(αs...),T}
                $poly{$(αs...),T}(coeffs, var)
            end
        end

        Base.length(ch::$poly{$(αs...),T,X}) where {$(αs...),T,X} = length(ch.coeffs)
        Polynomials.evalpoly(x, ch::$poly) = eval_cop(typeof(ch), ch.coeffs, x)
        (ch::$poly)(x) = evalpoly(x, ch)

        Base.convert(
            ::Type{P},
            q::Q,
        ) where {$(αs...),T,P<:$poly{$(αs...),T},Q<:$poly{$(αs...),T}} = q
        Base.convert(
            ::Type{$poly{$(αs...)}},
            q::Q,
        ) where {$(αs...),T,Q<:$poly{$(αs...),T}} = q
        Base.promote(
            p1::P,
            p2::Q,
        ) where {$(αs...),T,P<:$poly{$(αs...),T},Q<:$poly{$(αs...),T}} = p1, p2
        Base.promote_rule(
            ::Type{<:$poly{$(αs...),T}},
            ::Type{<:$poly{$(αs...),S}},
        ) where {$(αs...),T,S} = $poly{$(αs...),promote_type(T, S)}
        Base.promote_rule(
            ::Type{<:$poly{$(αs...),T}},
            ::Type{S},
        ) where {$(αs...),T,S<:Number} = $poly{$(αs...),promote_type(T, S)}

        $poly{$(αs...),T,X}(n::Number) where {$(αs...),T,X} = n * one($poly{$(αs...),T,X})
        $poly{$(αs...),T}(n::Number, var::Polynomials.SymbolLike=:x) where {$(αs...),T} =
            n * one($poly{$(αs...),T}, var)
        $poly{$(αs...)}(n::S, var::Polynomials.SymbolLike=:x) where {$(αs...),S<:Number} =
            $poly{$(αs...),S}(n, var)
        $poly{$(αs...),T}(var::Polynomials.SymbolLike=:x) where {$(αs...),T} =
            variable($poly{$(αs...),T}, var)
        $poly{$(αs...)}(var::Polynomials.SymbolLike=:x) where {$(αs...)} =
            variable($poly{$(αs...)}, var)

        Base.:+(
            p::P,
            q::Q,
        ) where {$(αs...),T,X,P<:$poly{$(αs...),T,X},S,Q<:$poly{$(αs...),S,X}} =
            ⊕(P, p, q)
        Base.:*(
            p::P,
            q::Q,
        ) where {$(αs...),T,X,P<:$poly{$(αs...),T,X},S,Q<:$poly{$(αs...),S,X}} =
            ⊗(P, p, q)
        Base.divrem(p1::$poly{$(αs...),T}, p2::$poly{$(αs...),T}) where {$(αs...),T} =
            _divrem(p1, p2)

        Polynomials._indeterminate(::Type{P}) where {$(αs...),T,X,P<:$poly{$(αs...),T,X}} =
            X
    end
end

## Delegation macros
"""
    ϟ(P) = Q; ϟ(P{T}) where {T} =   Q{T}

Used internally to delegate methods of `P` to `Q`.  (`[slash]upkoppa[tab]`)
"""
ϟ(P::Type{<:Polynomials.AbstractPolynomial}) =
    throw(ArgumentError("No default for  delegation"))

## Make monic version
## delegate to monic versions
macro register_monic(B)
    monic = esc(B)

    quote
        Polynomials.basis_symbol(::Type{P}) where {B<:$monic,P<:AbstractUnivariatePolynomial{B}} =
            Polynomials.basis_symbol(ImmutableDensePolynomial{ϟ(B)}) * "̃"
        SpecialPolynomials.ismonic(::Type{<:$monic}) = true
        SpecialPolynomials.abcde(B::Type{<:$monic}) = abcde(ϟ(B))
        SpecialPolynomials.weight_function(B::Type{<:$monic}) = weight_function(ϟ(B))
        Polynomials.domain(B::Type{<:$monic}) = domain(ϟ(B))
        SpecialPolynomials.B̃n(B::Type{<:$monic}, n::Int) = B̃n(ϟ(B), n)
        SpecialPolynomials.B̃n(B::Type{<:$monic}, v::Val{N}) where {N} = B̃n(ϟ(B), v)
        SpecialPolynomials.C̃n(B::Type{<:$monic}, n::Int) = C̃n(ϟ(B), n)
        SpecialPolynomials.C̃n(B::Type{<:$monic}, v::Val{N}) where {N} = C̃n(ϟ(B), v)
        SpecialPolynomials.ẫn(B::Type{<:$monic}, v::Val{N}) where {N} = ẫn(ϟ(B), v)
        SpecialPolynomials.b̂̃n(B::Type{<:$monic}, v::Val{N}) where {N} = b̂̃n(ϟ(B), v)
        SpecialPolynomials.ĉ̃n(B::Type{<:$monic}, v::Val{N}) where {N} = ĉ̃n(ϟ(B), v)
    end
end

## Make ortho
macro register_orthonormal(B)
    orthonormal = esc(B)

    quote
        Polynomials.basis_symbol(::Type{P}) where {B<:$orthonormal,P<:AbstractUnivariatePolynomial{B}} =
            Polynomials.basis_symbol(ImmutableDensePolynomial{ϟ(B)}) *  "̃"
        SpecialPolynomials.isorthonormal(::Type{<:$orthonormal}) = true
        SpecialPolynomials.abcde(B::Type{<:$orthonormal}) = abcde(ϟ(B))
        SpecialPolynomials.weight_function(B::Type{<:$orthonormal}) = weight_function(ϟ(B))
        Polynomials.domain(B::Type{<:$orthonormal}) = domain(ϟ(B))
        SpecialPolynomials.k0(B::Type{P}) where {P<:$orthonormal} =
            k0(ϟ(B)) / sqrt(norm2(ϟ(B), 0))
        SpecialPolynomials.k1k0(B::Type{P}, n::Int) where {P<:$orthonormal} =
            k1k0(ϟ(B), n) / ω₁₀(ϟ(B), n)
        SpecialPolynomials.ω₁₀(B::Type{<:$orthonormal}, n::Int) = one(eltype(P))

        SpecialPolynomials.B̃n(B::Type{<:$orthonormal}, n::Int) = B̃n(ϟ(B), n)
        SpecialPolynomials.B̃n(B::Type{<:$orthonormal}, v::Val{N}) where {N} = B̃n(ϟ(B), v)
        SpecialPolynomials.C̃n(B::Type{<:$orthonormal}, n::Int) = C̃n(ϟ(B), n)
        SpecialPolynomials.C̃n(B::Type{<:$orthonormal}, v::Val{N}) where {N} = C̃n(ϟ(B), v)
        SpecialPolynomials.ẫn(B::Type{<:$orthonormal}, v::Val{N}) where {N} =
            ẫn(ϟ(B), v)
        SpecialPolynomials.b̂̃n(B::Type{<:$orthonormal}, v::Val{N}) where {N} =
            b̂̃n(ϟ(B), v)
        SpecialPolynomials.ĉ̃n(B::Type{<:$orthonormal}, v::Val{N}) where {N} =
            ĉ̃n(ϟ(B), v)
    end
end

# P̃ = P(αx+β), so P̃'' = α^2P''(α x+ β); P̃' = αP'(αx+β)
macro register_shifted(B,  alpha, beta)
    shifted = esc(B)

    quote
        Polynomials.basis_symbol(::Type{P}) where {B<:$shifted,P<:AbstractUnivariatePolynomial{B}} =
            Polynomials.basis_symbol(ImmutableDensePolynomial{ϟ(B)}) *  "̃"

        function SpecialPolynomials.abcde(B::Type{<:$shifted})
            α, β = $alpha, $beta
            oP = 1
            a, b, c, d, e = abcde(ϟ(B))
            ã = oP * a
            b̃ = (oP * 2a * β) / α + (oP * b) / α
            c̃ = (oP * (a * β^2 + b * β + c)) / α^2
            d̃ = oP * d
            ẽ = (oP * (d * β + e)) / α
            (a=ã, b=b̃, c=c̃, d=d̃, e=ẽ)
        end
        function Polynomials.domain(B::Type{<:$shifted})
            α, β = $alpha, $beta
            a, b = extrema(ϟ(B))
            Polynomials.Interval((a - β) / α, (b - β) / α)
        end
        function SpecialPolynomials.k1k0(B::Type{<:$shifted}, n::Int)
            α, β = $alpha, $beta
            α * k1k0(ϟ(B), n)
        end
    end
end


# macro register_monic(name)
#     monic = esc(name)

#     quote
#         SpecialPolynomials.ismonic(::Type{<:$monic}) = true
#         SpecialPolynomials.abcde(P::Type{<:$monic}) = abcde(ϟ(P))
#         SpecialPolynomials.basis_symbol(P::Type{<:$monic}) = basis_symbol(ϟ(P)) * "̃"
#         SpecialPolynomials.weight_function(P::Type{<:$monic}) = weight_function(ϟ(P))
#         Polynomials.domain(P::Type{<:$monic}) = domain(ϟ(P))
#         SpecialPolynomials.B̃n(P::Type{<:$monic}, n::Int) = B̃n(ϟ(P), n)
#         SpecialPolynomials.B̃n(P::Type{<:$monic}, v::Val{N}) where {N} = B̃n(ϟ(P), v)
#         SpecialPolynomials.C̃n(P::Type{<:$monic}, n::Int) = C̃n(ϟ(P), n)
#         SpecialPolynomials.C̃n(P::Type{<:$monic}, v::Val{N}) where {N} = C̃n(ϟ(P), v)
#         SpecialPolynomials.ẫn(P::Type{<:$monic}, v::Val{N}) where {N} = ẫn(ϟ(P), v)
#         SpecialPolynomials.b̂̃n(P::Type{<:$monic}, v::Val{N}) where {N} = b̂̃n(ϟ(P), v)
#         SpecialPolynomials.ĉ̃n(P::Type{<:$monic}, v::Val{N}) where {N} = ĉ̃n(ϟ(P), v)
#     end
# end

# ## Make ortho
# macro register_orthonormal(name)
#     orthonormal = esc(name)

#     quote
#         SpecialPolynomials.isorthonormal(::Type{<:$orthonormal}) = true
#         SpecialPolynomials.abcde(P::Type{<:$orthonormal}) = abcde(ϟ(P))
#         SpecialPolynomials.basis_symbol(P::Type{<:$orthonormal}) = basis_symbol(ϟ(P)) * "̃"
#         SpecialPolynomials.weight_function(P::Type{<:$orthonormal}) = weight_function(ϟ(P))
#         Polynomials.domain(P::Type{<:$orthonormal}) = domain(ϟ(P))
#         SpecialPolynomials.k0(::Type{P}) where {P<:$orthonormal} =
#             k0(ϟ(P)) / sqrt(norm2(ϟ(P), 0))
#         SpecialPolynomials.k1k0(::Type{P}, n::Int) where {P<:$orthonormal} =
#             k1k0(ϟ(P), n) / ω₁₀(ϟ(P), n)
#         SpecialPolynomials.ω₁₀(P::Type{<:$orthonormal}, n::Int) = one(eltype(P))

#         SpecialPolynomials.B̃n(P::Type{<:$orthonormal}, n::Int) = B̃n(ϟ(P), n)
#         SpecialPolynomials.B̃n(P::Type{<:$orthonormal}, v::Val{N}) where {N} = B̃n(ϟ(P), v)
#         SpecialPolynomials.C̃n(P::Type{<:$orthonormal}, n::Int) = C̃n(ϟ(P), n)
#         SpecialPolynomials.C̃n(P::Type{<:$orthonormal}, v::Val{N}) where {N} = C̃n(ϟ(P), v)
#         SpecialPolynomials.ẫn(P::Type{<:$orthonormal}, v::Val{N}) where {N} =
#             ẫn(ϟ(P), v)
#         SpecialPolynomials.b̂̃n(P::Type{<:$orthonormal}, v::Val{N}) where {N} =
#             b̂̃n(ϟ(P), v)
#         SpecialPolynomials.ĉ̃n(P::Type{<:$orthonormal}, v::Val{N}) where {N} =
#             ĉ̃n(ϟ(P), v)
#     end
# end

# # P̃ = P(αx+β), so P̃'' = α^2P''(α x+ β); P̃' = αP'(αx+β)
# macro register_shifted(name, alpha, beta)
#     shifted = esc(name)
#     quote
#         function SpecialPolynomials.abcde(P::Type{<:$shifted})
#             α, β = $alpha, $beta
#             oP = one(eltype(P))
#             a, b, c, d, e = abcde(ϟ(P))
#             ã = oP * a
#             b̃ = (oP * 2a * β) / α + (oP * b) / α
#             c̃ = (oP * (a * β^2 + b * β + c)) / α^2
#             d̃ = oP * d
#             ẽ = (oP * (d * β + e)) / α
#             (a=ã, b=b̃, c=c̃, d=d̃, e=ẽ)
#         end

#         SpecialPolynomials.basis_symbol(P::Type{<:$shifted}) = basis_symbol(ϟ(P)) * "ˢ"
#         function Polynomials.domain(P::Type{<:$shifted})
#             α, β = $alpha, $beta
#             a, b = extrema(ϟ(P))
#             Polynomials.Interval((a - β) / α, (b - β) / α)
#         end
#         function SpecialPolynomials.k1k0(P::Type{<:$shifted}, n::Int)
#             α, β = $alpha, $beta
#             α * k1k0(ϟ(P), n)
#         end
#     end
# end

# Weight Functions
macro register_weight_function(W, P, ω)
    W′ = esc(W)
    P′ = esc(P)
    ω′ = esc(ω)

    quote
        SpecialPolynomials.:ϟ(::Type{<:$W′}) = $P′
        SpecialPolynomials.weight_function(::Type{<:$W′}) = $ω′
        (ch::$W′)(x) = clenshaw_eval(ch, x)
    end
end

macro register_discrete_weight_function(W, xs, ws)
    W′ = esc(W)
    xs′ = esc(xs)
    ws′ = esc(ws)

    quote
        SpecialPolynomials.xs_ws(::Type{<:$W′}) = ($xs′, $ws′)
        SpecialPolynomials.k0(::Type{<:$W′}) = 1
        (p::$W′)(x) = clenshaw_eval(typeof(p), coeffs(p), x)
    end
end
