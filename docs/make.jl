using Documenter
using GSplines

dir = "md"


makedocs(
    sitename = "GSplines",
    authors = "M. Marsala, B. Mourrain",
    modules = [GSplines],
    build = "GSplines.jl/docs",
    source = "md",
    pages = Any[
        "Home" => "index.md",
        "Functions & types" => [
            "geometric_model.md",
            "converters.md",
            "iga.md"
        ],
        "Using G+Smo" => "gismo.md",
        "About the package"  => "package.md",
    ],
    repo = Remotes.GitHub("AlgebraicGeometricModeling", "GSplines.jl"),
    doctest = false
)

#deploydocs(
#    repo = "github.com/AlgebraicGeometricModeling/GSplines.jl.git",
#    target = "site"
#)
