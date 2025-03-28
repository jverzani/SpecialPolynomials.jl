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

This package provides a number of different polynomial types to
represent polynomials, extending the `Polynomials` package.

```@docs
SpecialPolynomials.AbstractSpecialPolynomial
```

## Orthogonal polynomials

```@docs
SpecialPolynomials.AbstractOrthogonalPolynomial
SpecialPolynomials.AbstractCCOPBasis
SpecialPolynomials.AbstractCDOPBasis
```

## Implemented polynomial  types

### Classical continuous orthogonal polynomials

There are  several classical continuous  orthogonal polynomials available:

```@docs
Legendre
Chebyshev
ChebyshevU
Laguerre
Hermite
ChebyshevHermite
Gegenbauer
Jacobi
ShiftedJacobi
Bessel
ShiftedLegendre
```

### Classical discrete orthogonal polynomials

There are  several classical discrete  orthogonal polynomials available:

```@docs
Charlier
Krawchouk
Meixner
Hahn
DiscreteChebyshev
FallingFactorial
```


----

Some non-exported methods are available or define each of  the classical orthogonal polynomials:

```@docs
SpecialPolynomials.weight_function
SpecialPolynomials.generating_function
SpecialPolynomials.abcde
SpecialPolynomials.An
SpecialPolynomials.Bn
SpecialPolynomials.Cn
SpecialPolynomials.jacobi_matrix
SpecialPolynomials.gauss_nodes_weights
SpecialPolynomials.lagrange_barycentric_nodes_weights
```

### Defining new types

A new polynomial system  of classical type can  be specified fairly  succinctly,  provided the 5 constants  for the  `abcde`  method are known.

Polynomial systems  can also be generated  through  an associated weight function, but the code base currently needs some tending, so this feature is not enabled.

### Interpolating polynomials

```@docs
SpecialPolynomials.AbstractInterpolatingPolynomial
```

```@docs
Lagrange
Newton
```

### Other polynomials

```@docs
Bernstein
DualBernstein
```

#### Example of a [Bezier](https://pomax.github.io/bezierinfo/) curve (parameterized by `r(t) = ∑₀ᴺ bᵢBᵢ(t)`):

```@example
using Plots, Polynomials, SpecialPolynomials;  # hide
bs = [[220, 260], [220, 40], [35, 100],  [120, 140]]

p = Bernstein(bs)
ts = range(0, stop=1, length=50)
ps = p.(ts)
xs, ys = [[pᵢ[1] for pᵢ ∈ ps], [pᵢ[2] for pᵢ ∈ ps]]
p = plot(xs, ys, legend=false)
scatter!(p, [b[1] for b in bs], [b[2] for b in bs])
savefig("bezier.svg"); nothing  # hide
```

![](bezier.svg)
