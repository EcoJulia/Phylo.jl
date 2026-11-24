# SPDX-License-Identifier: BSD-2-Clause

using Pkg
Pkg.update()
"Phylo" ∈ [p.name for p in values(Pkg.dependencies())] &&
    Pkg.rm("Phylo")
# This checkout, not whatever `Pkg.develop("Phylo")` finds by name (on a CI
# runner, a fresh clone of the default branch, so a PR's docs were built
# without its changes)
Pkg.develop(path = dirname(@__DIR__))

using Documenter
using Phylo

makedocs(modules = [Phylo, Phylo.API],
         sitename = "Phylo.jl",
         pages = ["Home" => "index.md",
             "Tutorial / Quick start" => "tutorial.md",
             "Manual" => Any["Phylogeny data types" => "man/treetypes.md",
                             "Creating and writing phylogenies" => "man/io.md",
                             "Manipulating and building phylogenies" => "man/manipulating.md",
                             "Traversal and iterators" => "man/traversal.md",
                             "Getting phylogeny attributes" => "man/attributes.md",
                             "Modelling traits" => "man/modelling.md",
                             "Plotting" => "man/plotting.md"],
             "List of functions" => "functionlist.md",
             "API" => "api.md"];
         format = Documenter.HTML(size_threshold_ignore = ["man/plotting.md",
                                      "man/io.md"]))

deploydocs(repo = "github.com/EcoJulia/Phylo.jl.git",
           devbranch = "dev",
           push_preview = true)
