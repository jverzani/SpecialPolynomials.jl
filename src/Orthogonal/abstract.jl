##
## --------------------------------------------------
##
# Macros to register POLY{T},  POLY{α, T} and POLY{α, β, T}
#
# We use `Vector{T}` with N =  d+1 to store  a degree d  polynomial - no trailing zeros permitted.
#
macro register0(name, parent)
    poly = esc(name)
    parent_type = esc(parent)
    quote
        struct $poly{T,N} <: $parent_type{T,N}
            coeffs::Vector{T}
            var::Symbol
            function $poly{T,N}(coeffs::Vector{T}, var::Polynomials.SymbolLike=:x) where {T, N}
                M = length(coeffs)
                (M != N  || iszero(coeffs[end])) && throw(ArgumentError("wrong  size"))
                new{T,N}(coeffs, Symbol(var))
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
    
        Base.convert(::Type{P}, q::Q) where {T, P <:$poly{T}, Q <: $poly{T}} = q
        Base.promote(p::P, q::Q) where {T, P <:$poly{T}, Q <: $poly{T}} = p,q
        Base.promote_rule(::Type{<:$poly{T}}, ::Type{<:$poly{S}}) where {T,S} =
            $poly{promote_type(T, S)}
        Base.promote_rule(::Type{<:$poly{T}}, ::Type{S}) where {T,S<:Number} = 
            $poly{promote_type(T,S)}
        
        $poly{T}(n::Number, var::Polynomials.SymbolLike = :x) where {T} = $poly(T[n], var)
        $poly(n::S, var::Polynomials.SymbolLike = :x) where {S<:Number} = $poly{S}(n,var)
        $poly{T}(var::Polynomials.SymbolLike=:x) where {T} = variable($poly{T}, var)
        $poly(var::Polynomials.SymbolLike=:x) = variable($poly, var)

        # work around N parameter in promote(p,q) usage in defaults
        Base.:*(p::$poly{T}, q::$poly{T}) where {T}  = ⊗(p,q)
        Base.:+(p::$poly{T}, q::$poly{T}) where {T}  = ⊕(p,q)
        Base.divrem(p::$poly{T}, q::$poly{T}) where {T}  = _divrem(p,q)
    end
end

macro register1(name,parent)
    poly = esc(name)
    parent_type = esc(parent)
    quote
        struct $poly{α,T,N} <: $parent_type{α,T,N}
            coeffs::Vector{T}
            var::Symbol
            
            function $poly{α,T,N}(coeffs::Vector{T},  var::Polynomials.SymbolLike=:x) where {α,T,N}
                M =  length(coeffs)
                (M  != N || iszero(coeffs[end])) &&  throw(ArgumentError("wrong  size")) 
                new{α,T,N}(coeffs, Symbol(var))
            end
            
            function $poly{α,T}(coeffs::Vector{S},  var::Polynomials.SymbolLike=:x) where {α,T,S}
                N = findlast(!iszero, coeffs)
                if N ==  nothing
                    new{α,T,0}(T[], Symbol(var))
                else
                    new{α,T,N}(T.(coeffs[1:N]), Symbol(var))
                end
            end
            
            function $poly{α}(coeffs::Vector{T},  var::Polynomials.SymbolLike=:x) where {α, T}
                $poly{α,T}(coeffs, var)
            end
        end

        Base.length(p::$poly{α,T,N}) where {α,T,N} = N        
        (p::$poly)(x::S) where  {S} = eval_ccop(typeof(p),  p.coeffs, x)
    
        Base.convert(::Type{P}, q::Q) where {α,T, P<:$poly{α,T}, Q <: $poly{α,T}} = q
        Base.convert(::Type{$poly{α}}, q::Q) where {α,T, Q <: $poly{α,T}} = q        
        Base.promote(p::P, q::Q) where {α,T, P<:$poly{α,T}, Q <: $poly{α,T}} = p,q
        Base.promote_rule(::Type{<:$poly{α,T}}, ::Type{<:$poly{α,S}}) where {α,T,S} =
            $poly{α,promote_type(T, S)}
        Base.promote_rule(::Type{<:$poly{α,T}}, ::Type{S}) where {α,T,S<:Number} = 
            $poly{α,promote_type(T,S)}

        $poly{α,T}(n::Number, var::Polynomials.SymbolLike = :x) where {α,T} = $poly{α}(T[n], var)
        $poly{α}(n::S, var::Polynomials.SymbolLike = :x) where {α,S<:Number} = $poly{α,S}(n,var)
        $poly{α,T}(var::Polynomials.SymbolLike=:x) where {α, T} = variable($poly{α,T}, var)
        $poly{α}(var::Polynomials.SymbolLike=:x) where {α} = variable($poly{α}, var)

        # work around N parameter in promote
        Base.:*(p::$poly{α,T}, q::$poly{α,T}) where {α,T}  = ⊗(p,q)
        Base.:+(p::$poly{α,T}, q::$poly{α,T}) where {α,T}  = ⊕(p,q)
        Base.divrem(p::$poly{α,T}, q::$poly{α,T}) where {α,T}  = _divrem(p,q)

    end
end

macro register2(name, parent)
    poly = esc(name)
    parent_type = esc(parent)
    quote
        struct $poly{α,β,T,N} <: $parent_type{α,β,T,N}
            coeffs::Vector{T}
            var::Symbol
            
            function $poly{α,β,T,N}(coeffs::Vector{T},  var::Polynomials.SymbolLike=:x) where {α,β,T,N}
                M = length(coeffs)
                (M != N  || iszero(coeffs[end])) && throw(ArgumentError("wrong  size"))
                new{α,β,T,N}(coeffs, Symbol(var))
            end
            
            function $poly{α,β,T}(coeffs::AbstractVector{S},  var::Polynomials.SymbolLike=:x) where {α,β,T,S,M}
                N = findlast(!iszero, coeffs)
                if N ==  nothing
                    new{α,β,T,0}(T[], Symbol(var))
                else
                    cs = T[c  for c in coeffs[1:N]]
                    new{α,β,T,N}(cs,  Symbol(var))
                end
            end
            
            function $poly{α,β}(coeffs::AbstractVector{T},  var::Polynomials.SymbolLike=:x) where {α,β, T}
                $poly{α,β,T}(coeffs,  var)
            end
        end

        Base.length(p::$poly{α,β,T,N}) where {α,β,T,N} = N
        (p::$poly)(x::S) where  {S} = eval_ccop(typeof(p), p.coeffs, x)
    
        Base.convert(::Type{P}, q::Q) where {α,β,T,P<:$poly{α,β,T},Q<:$poly{α,β,T}} = q
        Base.convert(::Type{$poly{α,β}}, q::Q) where {α,β,T,Q<:$poly{α,β,T}} = q
        Base.promote(p::P, q::Q) where {α,β,T,P<:$poly{α,β,T},Q<:$poly{α,β,T}} = (p,q)
        Base.promote_rule(::Type{<:$poly{α,β,T}}, ::Type{<:$poly{α,β,S}}) where {α,β,T,S} =
            $poly{α,β,promote_type(T, S)}
        Base.promote_rule(::Type{<:$poly{α,β,T}}, ::Type{S}) where {α,β,T,S<:Number} = 
            $poly{α,β,promote_type(T,S)}

        $poly{α,β,T}(n::Number, var::Polynomials.SymbolLike = :x) where {α,β,T} = $poly{α,β}(T[n], var)
        $poly{α,β}(n::S, var::Polynomials.SymbolLike = :x) where {α,β,S<:Number} = $poly{α,β,S}(n,var)
        $poly{α,β,T}(var::Polynomials.SymbolLike=:x) where {α,β, T} = variable($poly{α,β,T}, var)
        $poly{α,β}(var::Polynomials.SymbolLike=:x) where {α,β} = variable($poly{α,β}, var)

        Base.:*(p::$poly{α,β,T}, q::$poly{α,β,T}) where {α,β,T}  = ⊗(p,q)
        Base.:+(p::$poly{α,β,T}, q::$poly{α,β,T}) where {α,β,T}  = ⊕(p,q)        
        Base.divrem(p::$poly{α,β,T}, q::$poly{α,β,T}) where {α,β,T}  = _divrem(p,q)

    end
end

macro register3(name, parent)
    poly = esc(name)
    parent_type = esc(parent)
    quote
        struct $poly{α,β,γ,T,N} <: $parent_type{α,β,γ,T,N}
            coeffs::Vector{T}
            var::Symbol
            
            function $poly{α,β,γ,T,N}(coeffs::Vector{T},  var::Polynomials.SymbolLike=:x) where {α,β,γ,T,N}
                M = length(coeffs)
                (M != N  || iszero(coeffs[end])) && throw(ArgumentError("wrong  size"))
                new{α,β,γ,T,N}(coeffs, Symbol(var))
            end
            
            function $poly{α,β,γ,T}(coeffs::AbstractVector{S},  var::Polynomials.SymbolLike=:x) where {α,β,γ,T,S,M}
                N = findlast(!iszero, coeffs)
                if N ==  nothing
                    new{α,β,γ,T,0}(T[], Symbol(var))
                else
                    cs = T[c  for c in coeffs[1:N]]
                    new{α,β,γ,T,N}(cs,  Symbol(var))
                end
            end
            
            function $poly{α,β,γ}(coeffs::AbstractVector{T},  var::Polynomials.SymbolLike=:x) where {α,β,γ, T}
                $poly{α,β,γ,T}(coeffs,  var)
            end
        end

        Base.length(p::$poly{α,β,γ,T,N}) where {α,β,γ,T,N} = N
        (p::$poly)(x::S) where  {S} = eval_ccop(typeof(p), p.coeffs, x)
    
        Base.convert(::Type{P}, q::Q) where {α,β,γ,T,P<:$poly{α,β,γ,T},Q<:$poly{α,β,γ,T}} = q
        Base.convert(::Type{$poly{α,β,γ}}, q::Q) where {α,β,γ,T,Q<:$poly{α,β,γ,T}} = q
        Base.promote(p::P, q::Q) where {α,β,γ,T,P<:$poly{α,β,γ,T},Q<:$poly{α,β,γ,T}} = (p,q)
        Base.promote_rule(::Type{<:$poly{α,β,γ,T}}, ::Type{<:$poly{α,β,γ,S}}) where {α,β,γ,T,S} =
            $poly{α,β,γ,promote_type(T, S)}
        Base.promote_rule(::Type{<:$poly{α,β,γ,T}}, ::Type{S}) where {α,β,γ,T,S<:Number} = 
            $poly{α,β,γ,promote_type(T,S)}

        $poly{α,β,γ,T}(n::Number, var::Polynomials.SymbolLike = :x) where {α,β,γ,T} = $poly{α,βγ,}(T[n], var)
        $poly{α,β,γ}(n::S, var::Polynomials.SymbolLike = :x) where {α,β,γ,S<:Number} = $poly{α,β,γ,S}(n,var)
        $poly{α,β,γ,T}(var::Polynomials.SymbolLike=:x) where {α,β,γ, T} = variable($poly{α,β,γ,T}, var)
        $poly{α,β,γ}(var::Polynomials.SymbolLike=:x) where {α,β,γ} = variable($poly{α,β,γ}, var)

        Base.:*(p::$poly{α,β,γ,T}, q::$poly{α,β,γ,T}) where {α,β,γ,T}  = ⊗(p,q)
        Base.:+(p::$poly{α,β,γ,T}, q::$poly{α,β,γ,T}) where {α,β,γ,T}  = ⊕(p,q)        
        Base.divrem(p::$poly{α,β,γ,T}, q::$poly{α,β,γ,T}) where {α,β,γ,T}  = _divrem(p,q)

    end
end
