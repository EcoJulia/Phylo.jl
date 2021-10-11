using Documenter
using Phylo

makedocs(modules = [Phylo, Phylo.API],
         sitename = "Phylo.jl", 
         pages = [
             "Home" => "index.md",
             "Tutorial / Quick start" => "tutorial.md"
             "Manual" => Any[
                 "Creating phylogenies" => "man/input.md",
                 "Manipulating and building phylogenies" => "man/manipulating.md",
                 "Traversal and iterators" => "man/traversal.md",
                 "Getting phylogeny attributes" => "man/attributes.md",
                 "Evolving traits" => "man/evolving.md",
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
