

function innerproduct(P::Type{<:DiscreteOrthogonalPolynomial}, f, g)
    dom = domain(P)
    fn = x -> f(x) * g(x) * weight_function(P)(x)
    a, b = first(dom), last(dom)
    if !isinf(a) && !isinf(b)
        return sum(fn(x)  for x in first(dom):last(dom))
    else
        ## what to do if infinite
    end
end

