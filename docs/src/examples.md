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

This package provides several  types to represent  polynomials relative  to
different  bases from the standard polynomial  basis, `1`,`x`,`x²`, `x³` etc.

For example, the Legendre polynomials are a collection of polynomials on `[-1,1]`. The first few may be seen through:

```jldoctest example
julia> using Polynomials, SpecialPolynomials

julia> p0 = Legendre([1])
Legendre(1⋅P₀(x))

julia> p1 = Legendre([0,1])
Legendre(1⋅P₁(x))

julia> p2 = Legendre([0,0,1])
Legendre(1⋅P₂(x))

julia> p3 = Legendre([0,0,0,1])
Legendre(1⋅P₃(x))
```

The coefficients, e.g., `[0,0,0,1]` indicate a polynomial `0⋅p0 +
0⋅p1 + 0⋅p2 + 1⋅p3`. The `show` method expresses these polynomials
relative to their bases. More familiar expressions are seen by
conversion to the standard basis. For example:


```jldoctest example
julia> convert.(Polynomial, [p0,p1,p2,p3])
4-element Array{Polynomial,1}:
 Polynomial(1)
 Polynomial(1.0*x)
 Polynomial(-0.5 + 1.5*x^2)
 Polynomial(-1.5*x + 2.5*x^3)
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
Polynomial(-1.5*x + 2.5*x^3)
```

Representation in another basis can be achieved this way:

```jldoctest example
julia> u = variable(ChebyshevU)
ChebyshevU(0.5⋅U₁(x))

julia> p3(u)
ChebyshevU(- 0.125⋅U₁(x) + 0.3125⋅U₃(x))
```

For most of the orthogonal polynomials, a conversion from the standard basis is provided, and a conversion between different parameter values  for the  same polynomial type are provded. Conversion methods between other polynomial types are not provided, but either evaluation, as above, or conversion through the `Polynomial` type is possible.


For the basis functions, the `basis` function can be used:

```jldoctest example
julia> h0,h1,h2,h3 = basis.(Hermite, 0:3);

julia> x = variable();

julia> h3(x)
Polynomial(-12.0*x + 8.0*x^3)
```

If the coefficients are known, they can be directly passed to the constructor:

```jldoctest example
julia> Laguerre{0}([1,2,3])
Laguerre{0}(1⋅L₀(x) + 2⋅L₁(x) + 3⋅L₂(x))
```

Some polynomial types are parameterized. The parameters are passed as in this example:

```jldoctest example
julia> Jacobi{1/2, -1/2}([1,2,3])
Jacobi{0.5,-0.5}(1⋅Jᵅᵝ₀(x) + 2⋅Jᵅᵝ₁(x) + 3⋅Jᵅᵝ₂(x))
```

----

The polynomial types specified above are orthogonal, meaning the
inner product of different basis vectors will be 0. For example:

```jldoctest example
julia> using QuadGK

julia> P = Legendre
Legendre

julia> p4,p5 = basis.(P, [4,5])
2-element Array{Legendre{Float64,N} where N,1}:
 Legendre(1.0⋅P₄(x))
 Legendre(1.0⋅P₅(x))
 
julia> wf, dom = SpecialPolynomials.weight_function(P), domain(P);

julia> quadgk(x -> p4(x) * p5(x) *  wf(x), first(dom), last(dom))
(-1.3877787807814457e-17, 0.0)
```

The unexported `innerproduct` will compute this as well, without the need to specifiy the domain or weight function, which can be gleaned from the type.

```jldoctest example
julia> SpecialPolynomials.innerproduct(P, p4, p5)
-1.3877787807814457e-17
```



## Polynomial methods

For each polynomial ttype, this package implements as many of the methods for polynomials defined in `Polynomials`, as possible. 


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

Multiplication formulas may not be defined for each type, and a fall back may be used where the multiplication is done with respect to the standard basis and the answer re-represented:

```jldoctest example
julia> P = Jacobi{1/2, -1/2}
Jacobi{0.5,-0.5,T,N} where N where T

julia> p,q = P([1,2]), P([-2,1])
(Jacobi{0.5,-0.5}(1⋅Jᵅᵝ₀(x) + 2⋅Jᵅᵝ₁(x)), Jacobi{0.5,-0.5}(- 2⋅Jᵅᵝ₀(x) + 1⋅Jᵅᵝ₁(x)))

julia> p * q
Jacobi{0.5,-0.5}(- 1.5⋅Jᵅᵝ₀(x) - 2.0⋅Jᵅᵝ₁(x) + 1.3333333333333333⋅Jᵅᵝ₂(x))
```

### Derivatives and integrals

Classic orthogonal polynomials  have known integral and derivative formulas. When implemented, these allow a direct computation, as is donee here:

```jldoctest example
julia> P = ChebyshevU{Float64}
ChebyshevU{Float64,N} where N

julia> p = P([1,2,3])
ChebyshevU(1.0⋅U₀(x) + 2.0⋅U₁(x) + 3.0⋅U₂(x))

julia> dp = derivative(p)
ChebyshevU(4.0⋅U₀(x) + 12.0⋅U₁(x))

julia> convert.(Polynomial, (p, dp))
(Polynomial(-2.0 + 4.0*x + 12.0*x^2), Polynomial(4.0 + 24.0*x))
```

There is a default  method defined through conversion to  the standard basis. It is  used behind the  scenes  in the next example.

```jldoctest example
julia> P = Jacobi{1//2, -1//2}
Jacobi{1//2,-1//2,T,N} where N where T

julia> p,q = P([1,2]), P([-2,1])
(Jacobi{1//2,-1//2}(1⋅Jᵅᵝ₀(x) + 2⋅Jᵅᵝ₁(x)), Jacobi{1//2,-1//2}(- 2⋅Jᵅᵝ₀(x) + 1⋅Jᵅᵝ₁(x)))

julia> p * q # as above, only with rationals for paramters
Jacobi{1//2,-1//2}(- 1.5⋅Jᵅᵝ₀(x) - 2.0⋅Jᵅᵝ₁(x) + 1.3333333333333333⋅Jᵅᵝ₂(x))

julia> P = Jacobi{1//2, 1//2}
Jacobi{1//2,1//2,T,N} where N where T

julia> p = P([1,2,3])
Jacobi{1//2,1//2}(1⋅Jᵅᵝ₀(x) + 2⋅Jᵅᵝ₁(x) + 3⋅Jᵅᵝ₂(x))

julia> dp = derivative(p)
Jacobi{1//2,1//2}(3.0⋅Jᵅᵝ₀(x) + 10.0⋅Jᵅᵝ₁(x))

julia> integrate(p)
Jacobi{1//2,1//2}(0.375⋅Jᵅᵝ₀(x) + 0.24999999999999994⋅Jᵅᵝ₁(x) + 0.6⋅Jᵅᵝ₂(x) + 0.5714285714285714⋅Jᵅᵝ₃(x))

julia> integrate(p, 0, 1)
25//8
```


### Roots

The `roots` function finds the roots of a polynomial.

```jldoctest example
julia> p = Legendre([1,2,2,1])
Legendre(1⋅P₀(x) + 2⋅P₁(x) + 2⋅P₂(x) + 1⋅P₃(x))

julia> rts = roots(p)
3-element Array{Float64,1}:
 -1.0
 -0.20000000000000007
  0.0

julia> p.(rts)
3-element Array{Float64,1}:
 -2.220446049250313e-16
  0.0
  0.0
```
 
 
Here we see `fromroots` and `roots` are related, provided a monic polynomial is used:
 
```jldoctest example
julia> using Polynomials, SpecialPolynomials; const SP=SpecialPolynomials
SpecialPolynomials

julia> P = Jacobi{1/2,-1/2}
Jacobi{0.5,-0.5,T,N} where N where T


julia> p = P([1,1,2,3])
Jacobi{0.5,-0.5}(1⋅Jᵅᵝ₀(x) + 1⋅Jᵅᵝ₁(x) + 2⋅Jᵅᵝ₂(x) + 3⋅Jᵅᵝ₃(x))

julia> q = SP.monic(p) # monic is not exported
Jacobi{0.5,-0.5}(0.13333333333333333⋅Jᵅᵝ₀(x) + 0.13333333333333333⋅Jᵅᵝ₁(x) + 0.26666666666666666⋅Jᵅᵝ₂(x) + 0.4⋅Jᵅᵝ₃(x))

julia> fromroots(P, roots(q)) - q |> u -> truncate(u, atol=sqrt(eps())) 
Jacobi{0.5,-0.5}(0.0)
```
 
The roots are found from the eigenvalues of the companion matrix,
which may be directly produced  by `companion`; the default is  to find  it after
conversion to the standard basis.
 
For orthogonal polynomials, the roots of the basis vectors are important for quadrature. For larger values of `n`, the eigenvalues of the unexported `jacobi_matrix` also identify these roots, but the algorithm is more stable.
 
```jldoctest example
julia> using LinearAlgebra

julia> p5 = basis(Legendre, 5)
Legendre(1.0⋅P₅(x))

julia> roots(p5)
5-element Array{Float64,1}:
 -0.9061798459386636
 -0.5384693101056829
  0.538469310105683
  0.9061798459386635
  0.0

julia> eigvals(SpecialPolynomials.jacobi_matrix(Legendre, 5))
5-element Array{Float64,1}:
 -0.906179845938664
 -0.5384693101056831
 -8.042985227174392e-17
  0.5384693101056834
  0.9061798459386639
```

For higher degree, the difference comes out:

```jldoctest example
julia> p50 = basis(Legendre{Float64}, 50); sum(isreal.(roots(p50)))
34

julia> eigvals(SpecialPolynomials.jacobi_matrix(Legendre, 50 ))  .|> isreal |> sum
50
```


The unexported `gauss_nodes_weights` function returns the nodes and weights. For some types (`Legendre`, `Hermite`, `Laguerre`) it uses an `O(n)` algorithm of Glaser, Liu, and Rokhlin, for others the `O(n²)` algorithm through the Jacobi matrix.

```jldoctest example
julia> xs, ys = SpecialPolynomials.gauss_nodes_weights(Legendre{Float64}, 10)
([-0.9739065285171717, -0.8650633666889845, -0.6794095682990244, -0.4333953941292472, -0.1488743389816312, 0.1488743389816312, 0.4333953941292472, 0.6794095682990244, 0.8650633666889845, 0.9739065285171717], [0.06667134430868804, 0.1494513491505808, 0.21908636251598188, 0.26926671930999635, 0.295524224714753, 0.295524224714753, 0.26926671930999635, 0.21908636251598188, 0.1494513491505808, 0.06667134430868804])
```

!!! note
    For a broader and more robust implementation, the [FastGaussQuadrature](https://github.com/JuliaApproximation/FastGaussQuadrature.jl) package provides `O(n)` algorithms for many classic orthogonal polynomial types.

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

The `Lagrange` and `Newton` types represent the polynomial in
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
polynomial to the function generating the `y` values.

For an orthogonal polynomial type, the zeros of the basis
polynomial `p_{n+1}`, labeled `x_0, x_1, ..., x_n` are often used as
nodes, especially for the Chebyshev nodes (of the first kind).  
[Gil, Segura, and Temme](https://archive.siam.org/books/ot99/OT99SampleChapter.pdf) say
"Interpolation with Chebyshev nodes is not as good as the best
approximation ..., but usually it is the best practical possibility
for interpolation and certainly much better than equispaced
interpolation"

For the orthogonal polynomial types, the default for `fit` for degree `n` will use the zeros of `P_{n+1}` to interpolate. We can see that some interpolation points lead to better fits than others, in this graphic:

```example
f(x) = exp(-x)*sinpi(x)
plot(f, -1, 1, legend=false, color=:black, linewidth=3)
p=fit(Val(:interpolating), Chebyshev, f, 3);  plot!(p, color=:blue)
p=fit(Val(:interpolating), ChebyshevU, f, 3); plot!(p, color=:red)
fit(Val(:interpolating), Legendre, f, 3);     plot!(p, color=:green)
xs = [-0.5, 0.0, 0.5]
p=fit(Newton, xs, f);
ts = range(-1, 1, length=100);                plot!(ts, p.(ts), color=:brown)
savefig("fitting.svg"); nothing # hide
```

![](fitting.svg)




### Polynomial approximation

There are other criteria for fitting that can be used.

least squares...
series...

XXXX


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

For the orthogonal polynomial types, fitting a polynomial to a function using least squares can be solved using the polynomial
`a0⋅p0 + a1⋅p1 + ⋅⋅⋅ + an⋅pn` where `ai=∫f⋅pi⋅w⋅dx / ∫pi^2⋅w⋅dx`. There is no need to specify values for `x`:

```jldoctest example
julia> f(x) = exp(-x) * sinpi(x)
f (generic function with 1 method)

julia> p = fit(Val(:lsq), Chebyshev{Float64}, f, 50);

julia> maximum(norm(p(x)-f(x) for x in range(-1,1,length=500))) <= sqrt(eps())
true
```





This wavy example is from Trefethen:

```jldoctest example
julia> f(x) = sin(6x) + sin(60*exp(x))
f (generic function with 1 method)

julia> p50 = fit(Val(:lsq), Chebyshev{Float64}, f, 50);

julia> maximum(norm(p50(x)-f(x) for x in range(-1,1,length=500))) <= sqrt(eps()) # cf. graph below
false
```

(With 50 points, the approximation misses badly over `[-1,1]`. There are 45 local extrema on  this interval.)

However, with more points we have a good fit:

```jldoctest example
julia> p196 = fit(Chebyshev{Float64}, f, 196);

julia> maximum(norm(p196(x)-f(x) for x in range(-1,1,length=500))) <= sqrt(eps())  # ≈ 1e-13
true
```

```@example
using Plots, Polynomials, SpecialPolynomials
f(x) = sin(6x) + sin(60*exp(x))
p50 = fit(Chebyshev{Float64}, f, 50);
p196 = fit(Chebyshev{Float64}, f, 196);
plot(f, -1, 1, legend=false, color=:black)
xs = range(-1, stop=1, length=500) # more points than recipe
plot!(xs, p50.(xs), color=:blue)
plot!(xs, p196.(xs), color=:red)
savefig("wavy.svg"); nothing # hide
```

![](wavy.svg)


!!! note
    The [ApproxFun](https://github.com/JuliaApproximation/ApproxFun.jl) package provides a framework to quickly and accuratately approximate functions using certain polynomial types. The choice of order and methods for most of Julia's built-in functions are conveniently provided.



## Plotting

The `plot` recipe from the `Polynomials` package works as expected for
the polynomial types in this package. The domain to be plotted over
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
