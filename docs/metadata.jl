# SPDX-License-Identifier: BSD-2-Clause

using Pkg

# Update Phylo folder packages 
Pkg.activate(".")
Pkg.update()

# Update examples folder packages
if isdir("examples")
    if isfile("examples/Project.toml")
        Pkg.activate("examples")
        Pkg.update()
        "Phylo" ∈ [p.name for p in values(Pkg.dependencies())] &&
            Pkg.rm("Phylo")
        Pkg.develop("Phylo")
    end
end

# Update docs folder packages
Pkg.activate("docs")
Pkg.update()
"Phylo" ∈ [p.name for p in values(Pkg.dependencies())] &&
    Pkg.rm("Phylo")
Pkg.develop("Phylo")

# Reformat files in package
using JuliaFormatter
using Phylo
format(Phylo)

# Carry out crosswalk for metadata
using ResearchSoftwareMetadata
ResearchSoftwareMetadata.crosswalk()
