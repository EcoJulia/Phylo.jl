using Documenter
using Phylo

makedocs(modules = [Phylo, Phylo.API],
         sitename = "Phylo.jl", 
         pages = [
             "Home" => "index.md",
             "Manual" => Any[
                 "Plotting" => "man/plotting.md"
             ],
             "API" => "api.md"
         ])

deploydocs(repo = "github.com/EcoJulia/Phylo.jl.git",
           devbranch = "dev",
           deps = Deps.pip("pygments",
                           "mkdocs",
                           "mkdocs-material",
                           "python-markdown-math"))
