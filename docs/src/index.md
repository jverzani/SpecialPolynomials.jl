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

## Orthogonal Polynomials

```@docs
SpecialPolynomials.AbstractOrthogonalPolynomial
SpecialPolynomials.AbstractCCOP
SpecialPolynomials.AbstractCDOP
```

There are  several classic continuous  orthogonal polynomials available:

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

There are  several classic discrete  orthogonal polynomials available:

```@docs
Charlier
Krawchouk
Meixner
Hahn
DiscreteChebyshev
FallingFactorial
```

Polynomial systems  can also be generated  through  an associated weight function.

```@docs
WeightFunction
DiscreteWeightFunction
```

Some  non-exported methods define each:

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

ρ = sum(bᵢ.*basis(Bernstein{N},i-1) for (i,bᵢ)  ∈ enumerate(bs))
ts = range(0, stop=1, length=500)
p =  plot(ρ[1].(ts), ρ[2].(ts), legend=false)
scatter!(p, [b[1] for b in bs], [b[2] for b in bs])
savefig("bezier.svg"); nothing # hide
```

![](bezier.svg)


