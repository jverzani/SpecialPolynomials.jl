@testset "Continous Weight" begin

    T = Float64

    ## known set of modified moments
    P = ShiftedLegendre{T}
    f(x) = -log(x)
    q0 = ContinuousWeight(P,f,[1], :x)
    WW = typeof(q0)
    function SpecialPolynomials.modified_moment(w::W, j) where {W <: WW}
        iszero(j) && return 1.0
        (-1)^j * gamma(j+1)^2/(j*(j+1)*gamma(2j+1))
    end

    q0 = ContinuousWeight(P, f, [1], :x)
    q1 = ContinuousWeight(P, f, [0,1], :x)
    q2 = ContinuousWeight(P, f, [0,0,1], :x)
    q3 = ContinuousWeight(P, f, [0,0,0,1], :x)
    q4 = ContinuousWeight(P, f, [0,0,0,0,1], :x)
    q5 = ContinuousWeight(P, f, [0,0,0,0,0,1], :x)
    qs = convert.(Polynomial, (q0,q1,q2,q3,q4,q5))
    a, b = first(domain(q0)), last(domain(q0))
    for i in eachindex(qs)
        for j in eachindex(qs)
            j <= i  && continue
            @test abs(SpecialPolynomials._quadgk(x -> qs[i](x) * qs[j](x) * f(x),  first(domain(q0)), last(domain(q0)))) <= 100*sqrt(eps(T))
        end
    end


    ## non-known set
    T = Float64
    P = ChebyshevTT{T}
    f(x) = SpecialPolynomials.weight_function(ChebyshevU{T})(x)
    q0 = ContinuousWeight(P, f, [1], :x)
    q1 = ContinuousWeight(P, f, [0,1], :x)
    q2 = ContinuousWeight(P, f, [0,0,1], :x)
    q3 = ContinuousWeight(P, f, [0,0,0,1], :x)
    q4 = ContinuousWeight(P, f, [0,0,0,0,1], :x)
    q5 = ContinuousWeight(P, f, [0,0,0,0,0,1], :x)
    qs = convert.(Polynomial, (q0,q1,q2,q3,q4,q5))
    a, b = first(domain(q0)), last(domain(q0))
    for qi in qs
        for qj in qs
            j <= i && continue
            @test abs(SpecialPolynomials._quadgk(x -> qi(x) * qj(x) * f(x), a, b)) <= 100*sqrt(eps(T))
        end
    end



end
