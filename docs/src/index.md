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
SpecialPolynomials.AbstractCCOP
SpecialPolynomials.AbstractCDOP
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
```

### Defining new types

A new polynomial system  of classical type can  be specified fairly  succinctly,  provided the 5 constants  for the  `abcde`  method are known.

Polynomial systems  can also be generated  through  an associated weight function.

```@docs
WeightFunction
DiscreteWeightFunction
```


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
```

Example of a [Bezier](https://pomax.github.io/bezierinfo/) curve  (parameterized by `r(t) = ∑₀ᴺ bᵢBᵢ(t)`:


```@example
using Plots, Polynomials, SpecialPolynomials
bs =[[220, 260], [220, 40], [35, 100],  [120, 140]]
N = length(bs)-1

ρ = sum(bᵢ.*basis(Bernstein{N},i-1) for (i,bᵢ)  ∈ enumerate(bs))
ts = range(0, stop=1, length=500)
p =  plot(ρ[1].(ts), ρ[2].(ts), legend=false)
scatter!(p, [b[1] for b in bs], [b[2] for b in bs])
savefig("bezier.svg"); nothing # hide
```

![](bezier.svg)


