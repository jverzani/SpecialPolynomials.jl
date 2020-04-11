# Examples


```@index
Pages = ["examples.md"]
```

```@meta
DocTestSetup = quote
  using Polynomials, SpecialPolynomials
end
```

## Construction

This package provides alternate families of bases to the standard basis for univariate polynomials. 
For example, the Legendre polynomials are a family of polynomials on `[-1,1]`. The first few may be seen through:

```jldoctest example
julia> using Polynomials, SpecialPolynomials

julia> p0 = Legendre([1])
Legendre(1⋅L_0(x))

julia> p1 = Legendre([0,1])
Legendre(1⋅L_1(x))

julia> p2 = Legendre([0,0,1])
Legendre(1⋅L_2(x))

julia> p3 = Legendre([0,0,0,1])
Legendre(1⋅L_3(x))
```

The coefficients, e.g., `[0,0,0,1]` indicate a polynomial `0⋅p0 +
0⋅p1 + 0⋅p2 + 1⋅p3`. The `show` method expresses these polynomials
relative to their bases. More familiar expressions are seen by
conversion to the standard basis. For example:



```jldoctest example
julia> convert.(Polynomial, [p0,p1,p2,p3])
4-element Array{Polynomial,1}:
 Polynomial(1)
 Polynomial(x)
 Polynomial(-1//2 + 3//2*x^2)
 Polynomial(-3//2*x + 5//2*x^3)
```

Polynomial instances are callable. We have, for example, to evaluate a polynomial at a set of points:

```jldoctest example
julia> p3.([1/4, 1/2, 3/4])
3-element Array{Float64,1}:
 -0.3359375
 -0.4375
 -0.0703125
```



Conversion can also be achieved through polynomial evaluation, using a variable `x` in the `Polynomial` basis:

```jldoctest example
julia> x = variable(Polynomial)
Polynomial(x)

julia> p3(x)
Polynomial(-3//2*x + 5//2*x^3)
```

Representation in another basis can be achieved this way:

```jldoctest example
julia> u = variable(ChebyshevU)
ChebyshevU(0.5⋅U_1(x))

julia> p3(u)
ChebyshevU(- 0.125⋅U_1(x) + 0.3125⋅U_3(x))
```

There are not conversion methods defined between each family, but
typically converting to `Polynomial` and then to the alternate family
will work.



For the basis functions, the unexported `Polynomials.basis` function can be used:

```jldoctest example
julia> h0,h1,h2,h3 = Polynomials.basis.(Hermite, 0:3);

julia> x = variable();

julia> h3(x)
Polynomial(-12*x + 8*x^3)
```

If the coefficients are known, they can be directly passed to the constructor:

```jldoctest example
julia> Laguerre([1,2,3])
Laguerre(1⋅L_0(x) + 2⋅L_1(x) + 3⋅L_2(x))
```

Some polynomial families are parameterized. The parameters are passed as in this example:

```jldoctest example
julia> Jacobi{1/2, -1/2}([1,2,3])
Jacobi(1⋅J^(α, β)_0(x) + 2⋅J^(α, β)_1(x) + 3⋅J^(α, β)_2(x))
```

----

The polynomial families specified above are orthogonal, meaning the
inner product of different basis vectors will be 0. For example:

```jldoctest example
julia> using QuadGK

julia> P = Legendre
Legendre

julia> p4,p5 = Polynomials.basis.(P, [4,5])
2-element Array{Legendre{Int64},1}:
 Legendre(1⋅L_4(x))
 Legendre(1⋅L_5(x))
 
julia> wf, dom = SpecialPolynomials.weight_function(P), domain(P);

julia> quadgk(x -> p4(x) * p5(x) *  wf(x), first(dom), last(dom))
(0.0, 0.0)
```

The unexported `innerproduct` will compute this as well, without the need to specifiy the domain or weight function, which can be gleaned from the type.

```jldoctest example
julia> SpecialPolynomials.innerproduct(P, p4, p5)
6.096184133406375e-16
```



## Polynomial methods

For each polynomial family, this package implements as many of the methods for polynomials defined in `Polynomials`, as possible. 


### Arithemtic

For example, basic arithmetic:

```
julia> P = ChebyshevU
ChebyshevU

julia> p,q = P([1,2,3,4]), P([-2,0,1,2])
(ChebyshevU(1⋅U_0(x) + 2⋅U_1(x) + 3⋅U_2(x) + 4⋅U_3(x)), ChebyshevU(- 2⋅U_0(x) + 1⋅U_2(x) + 2⋅U_3(x)))

julia> p + 1
ChebyshevU(2⋅U_0(x) + 2⋅U_1(x) + 3⋅U_2(x) + 4⋅U_3(x))

julia> -p
ChebyshevU(- 1⋅U_0(x) - 2⋅U_1(x) - 3⋅U_2(x) - 4⋅U_3(x))

julia> p + q
ChebyshevU(- 1⋅U_0(x) + 2⋅U_1(x) + 4⋅U_2(x) + 6⋅U_3(x))

julia> p*q
ChebyshevU(9⋅U_0(x) + 8⋅U_1(x) + 10⋅U_2(x) + 6⋅U_3(x) + 15⋅U_4(x) + 10⋅U_5(x) + 8⋅U_6(x))

julia> p^2
ChebyshevU(30⋅U_0(x) + 40⋅U_1(x) + 51⋅U_2(x) + 44⋅U_3(x) + 41⋅U_4(x) + 24⋅U_5(x) + 16⋅U_6(x))
```

Multiplication formulas may not be defined for each family, and a fall back may be used where the multiplication is done with respect to the standard basis and the answer re-represented:

```jldoctest example
julia> P = Jacobi{1/2, -1/2}
Jacobi{0.5,-0.5,T} where T<:Number

julia> p,q = P([1,2]), P([-2,1])
(Jacobi(1⋅J^(α, β)_0(x) + 2⋅J^(α, β)_1(x)), Jacobi(- 2⋅J^(α, β)_0(x) + 1⋅J^(α, β)_1(x)))

julia> p * q
Jacobi(- 1.5⋅J^(α, β)_0(x) - 2.0⋅J^(α, β)_1(x) + 1.3333333333333333⋅J^(α, β)_2(x))
```

### Derivatives and integrals

Many families have integral and derivative formulas that allow a direct computation:

```jldoctest example
julia> P = ChebyshevU{Float64}
ChebyshevU{Float64}

julia> p = P([1,2,3])
ChebyshevU(1.0⋅U_0(x) + 2.0⋅U_1(x) + 3.0⋅U_2(x))

julia> dp = derivative(p)
ChebyshevU(4.0⋅U_0(x) + 12.0⋅U_1(x))

julia> convert.(Polynomial, (p, dp))
(Polynomial(-2.0 + 4.0*x + 12.0*x^2), Polynomial(4.0 + 24.0*x))
```

There is a default  method defined through conversion to  the standard basis. It is  used behind the  scenes  in the next example.

```jldoctest example
julia> P = Jacobi{1//2, -1//2}
Jacobi{1//2,-1//2,T} where T<:Number

julia> p,q = P([1,2]), P([-2,1])
(Jacobi(1⋅J^(α, β)_0(x) + 2⋅J^(α, β)_1(x)), Jacobi(- 2⋅J^(α, β)_0(x) + 1⋅J^(α, β)_1(x)))

julia> p * q # as above, only with rationals for paramters
Jacobi(- 3//2⋅J^(α, β)_0(x) - 2//1⋅J^(α, β)_1(x) + 4//3⋅J^(α, β)_2(x))

julia> P = Jacobi{1//2, -1//2}
Jacobi{1//2,-1//2,T} where T<:Number

julia> p = P([1,2,3])
Jacobi(1⋅J^(α, β)_0(x) + 2⋅J^(α, β)_1(x) + 3⋅J^(α, β)_2(x))

julia> dp = derivative(p)
Jacobi(- 1//4⋅J^(α, β)_0(x) + 9//1⋅J^(α, β)_1(x))

julia> integrate(p)
Jacobi(1//16⋅J^(α, β)_0(x) + 15//16⋅J^(α, β)_1(x) + 11//12⋅J^(α, β)_2(x) + 3//5⋅J^(α, β)_3(x))

julia> integrate(p, 0, 1)
9//2
```


### Roots

The `roots` function finds the roots of a polynomial.

```jldoctest example
julia> p = Legendre([1,2,2,1])
Legendre(1⋅L_0(x) + 2⋅L_1(x) + 2⋅L_2(x) + 1⋅L_3(x))

julia> rts = roots(p)
3-element Array{Float64,1}:
 -1.0
 -0.2
  0.0

julia> p.(rts)
3-element Array{Float64,1}:
 0.0
 8.326672684688674e-17
 0.0
```
 
 
Here we see `fromroots` and `roots` are related, provided a monic polynomial is used:
 
```jldoctest example
julia> P = Jacobi{1/2,-1/2}
Jacobi{0.5,-0.5,T} where T<:Number


julia> p = P([1,1,2,3])
Jacobi(1⋅J^(α, β)_0(x) + 1⋅J^(α, β)_1(x) + 2⋅J^(α, β)_2(x) + 3⋅J^(α, β)_3(x))

julia> q = p / (convert(Polynomial,p)[end])   # monic
Jacobi(0.13333333333333333⋅J^(α, β)_0(x) + 0.13333333333333333⋅J^(α, β)_1(x) + 0.26666666666666666⋅J^(α, β)_2(x) + 0.4⋅J^(α, β)_3(x))

julia> fromroots(P, roots(q)) - q |> u -> round(u, digits=14)
Jacobi(0.0)
```
 
The roots are found from the eigenvalues of the companion matrix,
which may be directly produced  by `companion`; the default is  to find  it after
conversion to the standard basis.
 
For orthogonal polynomials, the roots of the basis vectors are important for quadrature. For larger values of `n`, the eigenvalues of the unexported `jacobi_matrix` also identify these roots, but the algorithm is more stable.
 
```jldoctest example
julia> using LinearAlgebra

julia> p5 = Polynomials.basis(Legendre, 5)
Legendre(1⋅L_5(x))

julia> roots(p5)
5-element Array{Float64,1}:
 -0.9061798459386644
 -0.5384693101056832
  0.5384693101056831
  0.9061798459386636
  0.0

julia> eigvals(SpecialPolynomials.jacobi_matrix(Legendre, 5))
5-element Array{Float64,1}:
 -0.9061798459386641
 -0.5384693101056831
  2.6689322308431607e-17
  0.5384693101056834
  0.9061798459386642
```

For higher degree, the difference comes out:

```jldoctest example
julia> p50 = Polynomials.basis(Legendre{Float64}, 50); sum(isreal.(roots(p50)))
38

julia> eigvals(SpecialPolynomials.jacobi_matrix(Legendre, 50 ))  .|> isreal |> sum
50
```


The unexported `gauss_nodes_weights` function returns the nodes and weights. For some families (`Legendre`, `Hermite`, `Laguerre`) it uses an ``O(n)` algorithm of Glaser, Liu, and Rokhlin, for others the `O(n²)` algorithm through the Jacobi matrix.

```jldoctest example
julia> xs, ys = SpecialPolynomials.gauss_nodes_weights(Legendre{Float64}, 10)
([-0.9739065285171701, -0.8650633666889836, -0.6794095682990245, -0.43339539412924644, -0.14887433898163127, 0.14887433898163216, 0.43339539412924777, 0.6794095682990245, 0.8650633666889843, 0.9739065285171717], [0.0666713443086869, 0.149451349150583, 0.2190863625159837, 0.26926671930999535, 0.2955242247147498, 0.29552422471475265, 0.26926671930999535, 0.21908636251598357, 0.14945134915058111, 0.06667134430868783])
```

!!! note
    For a broader and faster implementation, the [FastGaussQuadrature](https://github.com/JuliaApproximation/FastGaussQuadrature.jl) package provides `O(n)` algorithms for more families. 

## Fitting


### Interpolation

For any set of points `(x0,y0), (x1,y1), ..., (xn, yn)` with unique `x` values,  there is a unique polynomial of degree `n` or less that *interpolates* these points, that is  `p(x_i) = y_i`. The  `fit` function will perform polynomial interpolation:

```jldoctest example
julia> xs, ys = [0, 1/4,  1/2,  3/4], [1,2,2,3]
([0.0, 0.25, 0.5, 0.75], [1, 2, 2, 3])

julia> p1 = fit(Polynomial,  xs, ys)
Polynomial(1.0000000000000002 + 8.66666666666669*x - 24.000000000000085*x^2 + 21.333333333333385*x^3)

julia> p2 = fit(Lagrange, xs, ys)
Lagrange(1⋅ℓ^3_0(x) + 2⋅ℓ^3_1(x) + 2⋅ℓ^3_2(x) + 3⋅ℓ^3_3(x))

julia> p3 = fit(Newton, xs, ys)
Newton(1.0⋅p_0(x) + 4.0⋅p_1(x) - 8.0⋅p_2(x) + 21.333333333333332⋅p_3(x))
```

The `Lagrange` and `Newton` families represent the polynomial in
convenient bases based on the nodes (`xs`). These all represent the
same interpolating polynomial:

```jldoctest example
julia> [p1.(xs)-ys  p2.(xs)-ys p3.(xs)-ys]
4×3 Array{Float64,2}:
  2.22045e-16  0.0  0.0
  1.33227e-15  0.0  0.0
 -3.33067e-15  0.0  0.0
 -8.88178e-15  0.0  0.0
```


The `Lagrange` and `Newton`  methods allow a function to be specified  in place  of a set of `y` values:

```jldoctest example
julia> p = fit(Newton, [1,2,3], x->x^2)
Newton(1.0⋅p_0(x) + 3.0⋅p_1(x) + 1.0⋅p_2(x))

julia> convert(Polynomial, p)
Polynomial(1.0*x^2)
```

Polynomial interpolation can demonstrate the Runge phenomenon if the
nodes are evenly spaced. For higher degree fitting, the choice of
nodes can greatly effect the approximation of the interpolating
polynomial to the function generating the `y` values. Over the
interval `[-1,1]`, the zeros of the Chebyshev polynomial of the first
kind (`ChebyshevTT`) are well chosen nodes. They have move values near
the edge of the interval (they are the `x` values of points evenly
spaced on the upper half circle). These nodes are related to the zeros
of a trignoometric function, so can readily be generated without need
to call a root finding algorithm. The unexported `lagrange_barycentric_nodes_weights` method
will pass back these nodes (and corresponding weights used for quadrature):

```jldoctest example
julia> xs, ws = SpecialPolynomials.lagrange_barycentric_nodes_weights(ChebyshevTT, 50);

julia> f(x) = exp(-x)*sinpi(x)
f (generic function with 1 method)

julia> p = fit(Lagrange, xs, f);

julia> maximum(norm(p(x)-f(x) for x in range(-1,1,length=500))) <= sqrt(eps())  # ≈ 1e-14
true
```

Fitting with `xs = collect(range(-1, 1, length=50))` will have a maximum, as computed above, of `≈ 1e-4`.


### Polynomial approximation

If there are a lot of points, it is common  to  fit with a  lower  degree polynomial. This won't  be an interpolating polynomial, in general. The criteria  used to select the polynomial is  typically least squares (weighted least squares is also available). Adding a value for the degree, illustrates:

```jldoctest example
julia> xs, ys =  [1,2,3,4], [2.0,3,1,4]
([1, 2, 3, 4], [2.0, 3.0, 1.0, 4.0])

julia> p1 =  fit(Polynomial, xs,  ys, 1)  # degree 1  or less
Polynomial(1.5 + 0.4*x)

julia> p1 =  fit(Polynomial, xs,  ys, 2)  # degree 2 or less
Polynomial(4.000000000000018 - 2.100000000000017*x + 0.5000000000000026*x^2)

julia> p1 =  fit(Polynomial, xs,  ys)     # degree 3 or less (length(xs) - 1)
Polynomial(-10.000000000000302 + 20.16666666666734*x - 9.500000000000297*x^2 + 1.3333333333333721*x^3)
```

For the orthogonal families, fitting a polynomial to a function using least squares can be solved using the polynomial
`a0⋅p0 + a1⋅p1 + ⋅⋅⋅ + an⋅pn` where `ai=∫f⋅pi⋅w⋅dx / ∫pi^2⋅w⋅dx`. There is no need to specify values for `x`:

```jldoctest example
julia> f(x) = exp(-x) * sinpi(x)
f (generic function with 1 method)

julia> p = fit(ChebyshevTT{Float64}, f, 50);

julia> maximum(norm(p(x)-f(x) for x in range(-1,1,length=500))) <= sqrt(eps())
true
```

This wavy example is from Trefethen:

```jldoctest example
julia> f(x) = sin(6x) + sin(60*exp(x))
f (generic function with 1 method)

julia> p50 = fit(ChebyshevTT{Float64}, f, 50);

julia> maximum(norm(p50(x)-f(x) for x in range(-1,1,length=500))) <= sqrt(eps()) # cf. graph below
false
```

(With 50 points, the approximation misses badly over `[-1,1]`. There are 45 local extrema on  this interval.)

However, with more points we have a good fit:

```jldoctest example
julia> p196 = fit(ChebyshevTT{Float64}, f, 196);

julia> maximum(norm(p196(x)-f(x) for x in range(-1,1,length=500))) <= sqrt(eps())  # ≈ 1e-13
true
```

```@example
using Plots, Polynomials, SpecialPolynomials
f(x) = sin(6x) + sin(60*exp(x))
p50 = fit(ChebyshevTT{Float64}, f, 50);
p196 = fit(ChebyshevTT{Float64}, f, 196);
plot(f, -1, 1, legend=false, color=:black)
xs = range(-1, stop=1, length=500) # more points than recipe
plot!(xs, p50.(xs), color=:blue)
plot!(xs, p196.(xs), color=:red)
savefig("wavy.svg"); nothing # hide
```

![](wavy.svg)


!!! note
    The [ApproxFun](https://github.com/JuliaApproximation/ApproxFun.jl) package provides a framework to quickly and accuratately approximate functions using certain polynomial families. The choice of order and methods for most of Julia's built-in functions is conveniently provided.



## Plotting

The `plot` recipe from the `Polynomials` package works as expected for
the polynomial families in this package. The domain to be plotted over
matches that given by `domain`, unless this is infinite. A plot of the first few
Chebyshev Polynomials of the second kind can be produced as follows:

```@example
using Plots, Polynomials, SpecialPolynomials
# U1, U2, U3, and U4:
chebs = [
  ChebyshevU([0, 1]),
  ChebyshevU([0, 0, 1]),
  ChebyshevU([0, 0, 0, 1]),
  ChebyshevU([0, 0, 0, 0, 1]),
]
colors = ["#4063D8", "#389826", "#CB3C33", "#9558B2"]
itr = zip(chebs, colors)
(cheb,col), state = iterate(itr)
p = plot(cheb, c=col,  lw=5, legend=false, label="")
for (cheb, col) in Base.Iterators.rest(itr, state)
  plot!(cheb, c=col, lw=5)
end
savefig("chebs.svg"); nothing # hide
```

![](chebs.svg)
