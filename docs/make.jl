using Documenter
using SpecialPolynomials

ENV["PLOTS_TEST"] = "true"
ENV["GKSwstype"] = "100"

makedocs(
    sitename = "SpecialPolynomials",
    format = Documenter.HTML(ansicolor=true),
    modules = [SpecialPolynomials],
    checkdocs=:none
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
    )=#

deploydocs(
    repo = "github.com/jverzani/SpecialPolynomials.jl.git",
)
