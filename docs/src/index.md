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
Jacobi
Laguerre
Legendre
ShiftedLegendre
```

```@docs
DiscreteChebyshev
```

```@docs
WeightFunction
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
