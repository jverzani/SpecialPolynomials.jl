##
## --------------------------------------------------
##
# Macros to register POLY{T},  POLY{αs..., T}
#
# We use `Vector{T}` with N =  d+1 to store  a degree d  polynomial - no trailing zeros permitted.
##
macro register0(name, parent)
    poly = esc(name)
    parent_type = esc(parent)
    quote
        struct $poly{T,N} <: $parent_type{T,N}
            coeffs::Vector{T}
            var::Symbol
            function $poly{T,N}(coeffs::Vector{T}, var::Polynomials.SymbolLike=:x) where {T, N}
                M = length(coeffs)
                (M  != N || (N > 0 &&iszero(coeffs[end]))) &&  throw(ArgumentError("wrong  size")) 
                new{T,N}(coeffs, Symbol(var))

            end
            function $poly{T,N}(coeffs::NTuple{N,T}, var::Polynomials.SymbolLike=:x) where {T, N}
                M = length(coeffs)
                (M  != N || (N > 0 &&iszero(coeffs[end]))) &&  throw(ArgumentError("wrong  size")) 
                new{T,N}(collect(coeffs), Symbol(var))
            end

            function $poly{T}(coeffs::Vector{S},  var::Polynomials.SymbolLike=:x) where {T,S}
                N = findlast(!iszero, coeffs)
                if N ==  nothing
                    new{T,0}(T[], Symbol(var))
                else
                    cs = T.(coeffs[1:N])
                    new{T,N}(cs,  Symbol(var))
                end
            end
            
            function $poly(coeffs::Vector{S},  var::Polynomials.SymbolLike=:x) where {S}
                $poly{S}(coeffs, var)
            end

        end

        Base.length(p::$poly{T,N}) where {T,N} = N
        (p::$poly)(x::S) where  {S} = eval_ccop(typeof(p), p.coeffs, x)
        
        Polynomials.@register $poly

        # work around N parameter in promote(p,q) usage in defaults
        Base.:*(p::$poly{T}, q::$poly{T}) where {T}  = ⊗(p,q)
        Base.:+(p::$poly{T}, q::$poly{T}) where {T}  = ⊕(p,q)
        Base.divrem(p::$poly{T}, q::$poly{T}) where {T}  = _divrem(p,q)
    end
end

# register macro for polynomial families with parameters
# would be nice to consolidate, but ...
macro registerN(name,  parent, params...)
    poly = esc(name)
    parent_type = esc(parent)
    αs = tuple(esc.(params)...)
    quote
        struct $poly{$(αs...), T,N} <: $parent_type{$(αs...), T,N}
            coeffs::Vector{T}
            var::Symbol
            function $poly{$(αs...), T,N}(coeffs::Vector{T}, var::Polynomials.SymbolLike=:x) where {$(αs...), T, N}
                M = length(coeffs)
                (M  != N || (N > 0 &&iszero(coeffs[end]))) &&  throw(ArgumentError("wrong  size")) 
                new{$(αs...), T,N}(coeffs, Symbol(var))
            end

            function $poly{$(αs...),T,N}(coeffs::NTuple{N,T}, var::Polynomials.SymbolLike=:x) where {$(αs...), T, N}
                M = length(coeffs)
                (M  != N || (N > 0 &&iszero(coeffs[end]))) &&  throw(ArgumentError("wrong  size")) 
                new{$(αs...),T,N}(collect(coeffs), Symbol(var))
            end

            function $poly{$(αs...),T}(coeffs::Vector{S},  var::Polynomials.SymbolLike=:x) where {$(αs...),T,S}
                N = findlast(!iszero, coeffs)
                if N ==  nothing
                    new{$(αs...),T,0}(T[], Symbol(var))
                else
                    new{$(αs...),T,N}(T.(coeffs[1:N]), Symbol(var))
                end
            end

            function $poly{$(αs...)}(coeffs::Vector{T},  var::Polynomials.SymbolLike=:x) where {$(αs...), T}
                $poly{$(αs...),T}(coeffs, var)
            end

        end

        Base.length(p::$poly{$(αs...),T,N}) where {$(αs...),T,N} = N        
        (p::$poly)(x::S) where  {S} = eval_ccop(typeof(p),  p.coeffs, x)
    
        Base.convert(::Type{P}, q::Q) where {$(αs...),T, P<:$poly{$(αs...),T}, Q <: $poly{$(αs...),T}} = q
        Base.convert(::Type{$poly{$(αs...)}}, q::Q) where {$(αs...),T, Q <: $poly{$(αs...),T}} = q        
        Base.promote(p::P, q::Q) where {$(αs...),T, P<:$poly{$(αs...),T}, Q <: $poly{$(αs...),T}} = p,q
        Base.promote_rule(::Type{<:$poly{$(αs...),T}}, ::Type{<:$poly{$(αs...),S}}) where {$(αs...),T,S} =
            $poly{$(αs...),promote_type(T, S)}
        Base.promote_rule(::Type{<:$poly{$(αs...),T}}, ::Type{S}) where {$(αs...),T,S<:Number} = 
            $poly{$(αs...),promote_type(T,S)}

        $poly{$(αs...),T}(n::Number, var::Polynomials.SymbolLike = :x) where {$(αs...),T} =
            $poly{$(αs...)}(T[n], var)
        $poly{$(αs...)}(n::S, var::Polynomials.SymbolLike = :x)        where {$(αs...), S<:Number} =
            $poly{$(αs...),S}(n,var)
        $poly{$(αs...),T}(var::Polynomials.SymbolLike=:x)              where {$(αs...), T} =
            variable($poly{$(αs...),T}, var)
        $poly{$(αs...)}(var::Polynomials.SymbolLike=:x)                where {$(αs...)} =
            variable($poly{$(αs...)}, var)

        # work around N parameter in promote
        Base.:*(p::$poly{$(αs...),T}, q::$poly{$(αs...),T}) where {$(αs...),T}  = ⊗(p,q)
        Base.:+(p::$poly{$(αs...),T}, q::$poly{$(αs...),T}) where {$(αs...),T}  = ⊕(p,q)
        Base.divrem(p::$poly{$(αs...),T}, q::$poly{$(αs...),T}) where {$(αs...),T}  = _divrem(p,q)

        
    end
end

## Delegation macros
"""
    ϟ(P) = Q; ϟ(P{T}) where {T} =   Q{T}

Used to delegate methods of P to Q.  (upkoppa[tab])
"""
ϟ(P::Type{<:Polynomials.AbstractPolynomial}) = throw(ArgumentError("No default for  delegation"))

## Make monic version
macro register_monic(name)

    monic=esc(name)
    
    quote
        SpecialPolynomials.ismonic(::Type{<:$monic}) = true
        SpecialPolynomials.abcde(P::Type{<:$monic})  = abcde(ϟ(P))
        SpecialPolynomials.basis_symbol(P::Type{<:$monic})  = basis_symbol(ϟ(P)) * "̃"
        SpecialPolynomials.weight_function(P::Type{<:$monic}) = weight_function(ϟ(P))
        Polynomials.domain(P::Type{<:$monic}) = domain(ϟ(P))
        SpecialPolynomials.B̃n(P::Type{<:$monic}, n::Int) =  B̃n(ϟ(P),n)
        SpecialPolynomials.B̃n(P::Type{<:$monic}, v::Val{N}) where{N} =  B̃n(ϟ(P),v)
        SpecialPolynomials.C̃n(P::Type{<:$monic}, n::Int) =  C̃n(ϟ(P),n)
        SpecialPolynomials.C̃n(P::Type{<:$monic}, v::Val{N}) where{N} =  C̃n(ϟ(P),v)
        SpecialPolynomials.ẫn(P::Type{<:$monic}, v::Val{N}) where{N} =  ẫn(ϟ(P),v)
        SpecialPolynomials.b̂̃n(P::Type{<:$monic}, v::Val{N}) where{N} =  b̂̃n(ϟ(P),v)
        SpecialPolynomials.ĉ̃n(P::Type{<:$monic}, v::Val{N}) where{N} =  ĉ̃n(ϟ(P),v)
    end
end

# P̃ = P(αx+β), so P̃'' = α^2P''(α x+ β); P̃' = αP'(αx+β)
macro register_shifted(name, alpha, beta)
    shifted = esc(name)
    quote
        function SpecialPolynomials.abcde(P::Type{<:$shifted}) 
            α,β = $alpha, $beta
            oP = one(eltype(P))
            a,b,c,d,e = abcde(ϟ(P))
            ã = oP*a
            b̃ = (oP * 2a*β)/α + (oP*b)/α
            c̃ = (oP * (a*β^2 +b*β + c))/α^2
            d̃ = oP *  d
            ẽ = (oP * (d*β + e))/α
            NamedTuple{(:a,:b,:c,:d,:e)}((ã, b̃, c̃, d̃, ẽ))
        end
        
        SpecialPolynomials.basis_symbol(P::Type{<:$shifted})  = basis_symbol(ϟ(P)) * "ˢ"
        function  Polynomials.domain(P::Type{<:$shifted})
            α,β = $alpha, $beta
            a, b = extrema(ϟ(P))
            Polynomials.Interval((a- β)/α, (b - β)/α)
        end
        function SpecialPolynomials.kn(P::Type{<:$shifted}, n::Int)
            α,β = $alpha, $beta
            α^n * kn(ϟ(P))
        end
        function SpecialPolynomials.k1k0(P::Type{<:$shifted}, n::Int)
            α,β = $alpha, $beta
            α * k1k0(ϟ(P), n)
        end
        function  SpecialPolynomials.k1k_1(P::Type{<:$shifted}, n::Int)
            α,β = $alpha, $beta
            α^2 * k1k_1(ϟ(P),n)
        end
        
    end
end

# Weight Functions
macro register_weight_function(W, P, ω)
    W′ = esc(W)
    P′ = esc(P)
    ω′ = esc(ω)
    
    quote
        SpecialPolynomials.:ϟ(::Type{<:$W′}) =  $P′
        SpecialPolynomials.weight_function(::Type{<:$W′}) =  $ω′
        (ch::$W′)(x)  = clenshaw_eval(ch, x)
    end
end

macro register_discrete_weight_function(W, xs, ws)
    W′ = esc(W)
    xs′ = esc(xs)
    ws′ = esc(ws)
    
    quote
        SpecialPolynomials.xs_ws(::Type{<:$W′}) = ($xs′, $ws′)
        (p::$W′)(x)  = clenshaw_eval(p, x)
    end
end
