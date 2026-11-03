# SPDX-License-Identifier: BSD-2-Clause

using Phylo
using DataFrames
using LinearAlgebra
using CSV
using Statistics
using ForwardDiff
using DynamicPPL

# load tip data
df = DataFrame(CSV.File("Data/Myrtaceae.csv"))

# read new tree
const tree::TraitTree{13} = open(f -> parsenewick(f, TraitTree{13}),
                                   "Data/Qian2016.tree");

# remove missing species from dataframe
dropmissing!(df, :species)

# add underscores to dataframe
df.species = replace.(df.species, " " => "_")

# filter tree to just species we have data for
keep = intersect(getleafnames(tree), df.species)
keeptips!(tree, keep)                                  
filter!(:species => x -> x ∈ keep, df)

# get mean data for trait value
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

# standardise the data
tmin = (dat.tmin .- mean(dat.tmin)) ./ std(dat.tmin)
tmax = (dat.tmax .- mean(dat.tmax)) ./ std(dat.tmax)
trng = (dat.trng .- mean(dat.trng)) ./ std(dat.trng)
stl1 = (dat.stl1 .- mean(dat.stl1)) ./ std(dat.stl1)
stl2 = (dat.stl2 .- mean(dat.stl2)) ./ std(dat.stl2)
stl3 = (dat.stl3 .- mean(dat.stl3)) ./ std(dat.stl3)
stl4 = (dat.stl4 .- mean(dat.stl4)) ./ std(dat.stl4)
swvl1 = (dat.swvl1 .- mean(dat.swvl1)) ./ std(dat.swvl1)
swvl2 = (dat.swvl2 .- mean(dat.swvl2)) ./ std(dat.swvl2)
swvl3 = (dat.swvl3 .- mean(dat.swvl3)) ./ std(dat.swvl3)
swvl4 = (dat.swvl4 .- mean(dat.swvl4)) ./ std(dat.swvl4)
ssr = (dat.ssr .- mean(dat.ssr)) ./ std(dat.ssr)
tp = (dat.tp .- mean(dat.tp)) ./ std(dat.tp)

data = DataFrame(species = keep, tmin = tmin, tmax = tmax, trng = trng, stl1 = stl1, stl2 = stl2, stl3 = stl3, stl4 = stl4, swvl1 = swvl1, swvl2 = swvl2, swvl3 = swvl3, swvl4 = swvl4, ssr = ssr, tp = tp)

# add the trait data to the tree
for i in eachrow(data)
    setnodedata!(tree, i.species, Phylo.traitdata(Union{Float64, ForwardDiff.Dual{ForwardDiff.Tag{DynamicPPL.DynamicPPLTag, Float64}, Float64, 3}}, ["tmin", "tmax", "trng", "stl1", "stl2", "stl3", "stl4", "swvl1", "swvl2", "swvl3", "swvl4", "ssr", "tp"], [i.tmin, i.tmax, i.trng, i.stl1, i.stl2, i.stl3, i.stl4, i.swvl1, i.swvl2, i.swvl3, i.swvl4, i.ssr, i.tp]))
end

trait = ["tmin", "tmax", "trng", "stl1", "stl2", "stl3", "stl4", "swvl1", "swvl2", "swvl3", "swvl4", "ssr", "tp"]
nodes = getnodes(tree, postorder)

# add lengths to tree
for node in nodes
    val = getnodedata(tree, node).value
    if hasinbound(tree, node)
        len = Phylo.getlength(tree, Phylo.getinbound(tree, node))
        td = traitdata(eltype(nodedatatype(typeof(tree))), trait, val, len)
        setnodedata!(tree, node, td)
    else
        td = traitdata(eltype(nodedatatype(typeof(tree))), trait, val)
        setnodedata!(tree, node, td)
    end
end

zvec = [getnodedata(tree, leaf).value for leaf in keep]
z = reduce(vcat, zvec)