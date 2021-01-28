using Documenter
using Phylo

makedocs(modules = [Phylo, Phylo.API],
         sitename = "Phylo.jl")

deploydocs(deps = Deps.pip("pygments",
                           "mkdocs",
                           "mkdocs-material",
                           "python-markdown-math"),
           repo = "github.com/EcoJulia/Phylo.jl.git")
