using SpecialPolynomials
using Plots
using FileIO

ts = range(0, 2pi, 100);
xys = sincos.(ts);
xs,ys=first.(xys), last.(xys);

p = plot(; legend=false,
         aspect_ratio=:equal,
         showaxis = false,
         background_color = :transparent)

julia_colors = (:brown3, :mediumorchid, :forestgreen, :royalblue)
pts = ((0,0), (3, 0), (3/2,2))

for (col, (a,b)) ∈ zip(julia_colors, pts)
    plot!(a .+xs, b .+ys; seriestype=:shape, color=col, alpha=0.25)
end

xs = range(-1,1, 100)
P = Legendre
degs = (3,4,5)
ps = basis.(P, degs)
for (p, (a,b)) ∈  zip(ps, pts)
    ys = p.(xs)
    ys = [yᵢ^2 <= 1 - xᵢ^2 ? yᵢ : NaN for (xᵢ,yᵢ) ∈ zip(xs, ys)]
    plot!([a-1,a+1], [b,b], linewidth=3, color=:white, alpha=0.25)
    plot!(a .+ xs, b .+ ys; linewidth=3, color=:black)
end


for (d, (a,b)) ∈ zip(degs, pts)
    rts, ws = SpecialPolynomials.gauss_nodes_weights(P, d)
    zs = zero.(rts)
    scatter!(a .+ rts, b .+ zs, markersize=5, color=last(julia_colors))
end

p

save("src/assets/logo.png", p)
