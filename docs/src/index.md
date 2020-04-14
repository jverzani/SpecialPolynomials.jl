# SpecialPolynomials.jl

Documentation for SpecialPolynomials.jl.



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
ChebyshevTT
ChebyshevU
Gegenbauer
GeneralizedLaguerre
Hermite
ChebyshevHermite
Jacobi
Laguerre
Legendre
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

Example of a Bezier curve  (parameterized by `r(t) = ∑₀ᴺ bᵢBᵢ(t)`:


```@example
using Plots, Polynomials, SpecialPolynomials
bs =[[220,40], [220,260], [35,200],  [120,160]]
N = length(bs)-1
r(t) = sum(b.*B(t) for (b,B) in zip(bs, Polynomials.basis.(Bernstein{N},0:N)))
ts = range(0, stop=1, length=500); rs=r.(ts)
plot([r[1] for r in rs], [r[2] for r in rs], legend=false)
savefig("bezier.svg"); nothing # hide
```

![](bezier.svg)


