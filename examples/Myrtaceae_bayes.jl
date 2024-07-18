# SPDX-License-Identifier: BSD-2-Clause

using Phylo
using Turing
using DataFrames
using CSV
using Plots

# load the data
df = DataFrame(CSV.File("Data/Myrtaceae.csv"))

# load the tree
const traittree::TraitTree{1} = open(f -> parsenewick(f, TraitTree{1}),
                                     "Data/Qian2016.tree")

# remove missing species from dataframe
dropmissing!(df, :species)

# add underscores to dataframe
df.species = replace.(df.species, " " => "_")

# filter the tree and dataframe for the species that are in  both
keep = intersect(getleafnames(traittree), df.species)
keeptips!(traittree, keep)
filter!(:species => x -> x ∈ keep, df)

# use mean data for trait value
gdf = groupby(df, :species)
dat = combine(gdf,
              [
                  :tmin,
                  :tmax,
                  :trng,
                  :stl1,
                  :stl2,
                  :stl3,
                  :stl4,
                  :swvl1,
                  :swvl2,
                  :swvl3,
                  :swvl4,
                  :ssr,
                  :tp
              ] .=> mean; renamecols = false)

# add the data for tmin to the tmintree
for i in eachrow(dat)
    setnodedata!(traittree, i.species,
                 Phylo.traitdata(eltype(nodedatatype(typeof(traittree))),
                                 ["tmin"], [i.tmin]))
end

trait = ["tmin"]
# add lengths to tree
nodes = getnodes(traittree, postorder)
for node in nodes
    val = getnodedata(traittree, node).value
    if hasinbound(traittree, node)
        len = Phylo.getlength(traittree, Phylo.getinbound(traittree, node))
        td = traitdata(eltype(nodedatatype(typeof(traittree))), trait, val, len)
        setnodedata!(traittree, node, td)
    else
        td = traitdata(eltype(nodedatatype(typeof(traittree))), trait, val)
        setnodedata!(traittree, node, td)
    end
end

@model function βσ()
    β ~ Uniform(-500, 500)
    σ ~ Uniform(0, 1000)
    return β, σ
end

@model function βσ_threepoint(tree, z) # z needs to be for leaves in postorder
    @submodel β, σ = βσ()
    z ~ Phylo.MyDist2(σ, β, tree) # tree.z ~ (implement later)
    return nothing
end

@model function βσλ_threepoint(tree, z, upper = 1.0) # z needs to be for leaves in postorder
    @submodel β, σ = βσ()
    λ ~ Uniform(0, upper)
    z ~ Phylo.MyDist3(σ, β, λ, tree)
    return nothing
end

z = dat.tmin

spl = sample(βσ_threepoint(traittree, z), HMC(0.001, 10), 5000,
             initial_params = [290.0, 8.0]) # add initial_params

spl = sample(βσλ_threepoint(traittree, z), HMC(0.001, 10), 5000,
             initial_params = [290.0, 8.0, 0.8])

for i in eachrow(dat)
    setnodedata!(traittree, i.species,
                 Phylo.traitdata(eltype(nodedatatype(typeof(traittree))),
                                 ["tmax"], [i.tmax]))
end

z = dat.tmax

spl = sample(βσ_threepoint(traittree, z), HMC(0.01, 10), 5000,
             initial_params = [290.0, 8.0])

spl = sample(βσλ_threepoint(traittree, z), HMC(0.01, 10), 5000,
             initial_params = [290.0, 8.0, 0.8])

for i in eachrow(dat)
    setnodedata!(traittree, i.species,
                 Phylo.traitdata(eltype(nodedatatype(typeof(traittree))),
                                 ["trng"], [i.trng]))
end

z = dat.trng

spl = sample(βσ_threepoint(traittree, z), HMC(0.01, 10), 5000,
             initial_params = [4.0, 1.0])

spl = sample(βσλ_threepoint(traittree, z), HMC(0.01, 10), 5000,
             initial_params = [4.0, 1.0, 0.8])

for i in eachrow(dat)
    setnodedata!(traittree, i.species,
                 Phylo.traitdata(eltype(nodedatatype(typeof(traittree))),
                                 ["swvl1"], [i.swvl1]))
end

z = dat.swvl1

estimaterates(traittree, ["swvl1"], lambda = 0.5)

spl = sample(βσ_threepoint(traittree, z), HMC(0.0001, 10), 5000,
             initial_params = [0.1, 1.0])

spl = sample(βσλ_threepoint(traittree, z), HMC(0.0001, 10), 5000,
             initial_params = [0.1, 1.0, 0.5])
plot(spl[:λ])
