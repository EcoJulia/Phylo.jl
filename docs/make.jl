using Documenter
using Phylo

makedocs(modules = [Phylo, Phylo.API],
         sitename = "Phylo.jl", 
         pages = ["Home" => "index.md",
                  "Tutorial / Quick start" => "tutorial.md",
                  "Manual" => Any[
                    "Phylogeny data types" => "man/treetypes.md",
                    "Creating phylogenies" => "man/input.md",
                    "Manipulating and building phylogenies" => "man/manipulating.md",
                    "Traversal and iterators" => "man/traversal.md",
                    "Getting phylogeny attributes" => "man/attributes.md",
                    "Modelling traits" => "man/modelling.md",
                    "Plotting" => "man/plotting.md"
                  ],
                  "List of functions" => "functionlist.md",
                  "API" => "api.md"];
         format = Documenter.HTML(size_threshold_ignore = ["man/plotting.md",
                                                           "man/input.md"]))

deploydocs(repo = "github.com/EcoJulia/Phylo.jl.git",
           devbranch = "dev",
           deps = Deps.pip("pygments",
                           "mkdocs",
                           "mkdocs-material",
                           "python-markdown-math"))
