using Documenter
using SpecialPolynomials

makedocs(
    sitename = "SpecialPolynomials",
    format = Documenter.HTML(),
    modules = [SpecialPolynomials]
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
