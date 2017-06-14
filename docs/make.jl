using Documenter
using Phylo

makedocs(modules = [Phylo, Phylo.API],
         clean   = false)

deploydocs(deps = Deps.pip("pygments",
                           "mkdocs",
                           "mkdocs-material",
                           "python-markdown-math"),
           repo = "github.com/richardreeve/Phylo.jl.git",
           julia="0.6",
           osname="linux")
