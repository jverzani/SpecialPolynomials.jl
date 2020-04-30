# SpecialPolynomials.jl

Documentation for [SpecialPolynomials.jl](https://github.com/jverzani/SpecialPolynomials.jl).



```@meta
CurrentModule = SpecialPolynomials
```


```@meta
DocTestSetup = quote
  using Polynomials, SpecialPolynomials
end
```

## Overview

This package provides a number of different polynomial families to
represent polynomials, extending the `Polynomials` package.

```@docs
SpecialPolynomials.AbstractSpecialPolynomial
```

## Orthogonal Polynomials

```@docs
SpecialPolynomials.AbstractOrthogonalPolynomial
```

There are  several families of orthogonal polynomials available.

```@docs
Legendre
Chebyshev
ChebyshevU
Laguerre
Hermite
ChebyshevHermite
Gegenbauer
Jacobi
GeneralizedLaguerre
ShiftedLegendre
```

```@docs
DiscreteChebyshev
Krawtchouk
```

```@docs
WeightFunction
DiscreteWeightFunction
```

Some  non-exported methods define each:

```@docs
SpecialPolynomials.weight_function
SpecialPolynomials.generating_function
SpecialPolynomials.An
SpecialPolynomials.Bn
SpecialPolynomials.Cn
SpecialPolynomials.alpha
SpecialPolynomials.beta
SpecialPolynomials.jacobi_matrix
SpecialPolynomials.gauss_nodes_weights
```


## Interpolating Polynomials

```@docs
SpecialPolynomials.AbstractInterpolatingPolynomial
```

```@docs
Lagrange
Newton
```

## Other Polynomials

```@docs
Bernstein
```

Example of a [Bezier](https://pomax.github.io/bezierinfo/) curve  (parameterized by `r(t) = ∑₀ᴺ bᵢBᵢ(t)`:


```@example
using Plots, Polynomials, SpecialPolynomials
bs =[[220, 260], [220, 40], [35, 100],  [120, 140]]
N = length(bs)-1
ρ = Bernstein(bs) 
ts = range(0, stop=1, length=500)
rs= ρ.(ts);
p =  plot([r[1] for r in rs], [r[2] for r in rs], legend=false)
scatter!(p, [b[1] for b in bs], [b[2] for b in bs])
savefig("bezier.svg"); nothing # hide
```

![](bezier.svg)


