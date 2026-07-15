# SPDX-License-Identifier: BSD-2-Clause

using Phylo
using DataFrames
using LinearAlgebra
using CSV
using Statistics
using ForwardDiff
using DynamicPPL

# load tip data
df6 = DataFrame(CSV.File("Data/Myrtaceae.csv"))

# read new tree
const tree6::TraitTree{6} = open(f -> parsenewick(f, TraitTree{6}),
                                   "Data/Qian2016.tree");

# remove missing species from dataframe
dropmissing!(df6, :species)

# add underscores to dataframe
df6.species = replace.(df6.species, " " => "_")

# filter tree to just species we have data for
keep = intersect(getleafnames(tree6), df6.species)
keeptips!(tree6, keep)                                  
filter!(:species => x -> x ∈ keep, df6)

# get mean data for trait value
gdf = groupby(df6, :species)
dat6 = combine(gdf,
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
tmin = (dat6.tmin .- mean(dat6.tmin)) ./ std(dat6.tmin)
tmax = (dat6.tmax .- mean(dat6.tmax)) ./ std(dat6.tmax)
stl1 = (dat6.stl1 .- mean(dat6.stl1)) ./ std(dat6.stl1)
swvl1 = (dat6.swvl1 .- mean(dat6.swvl1)) ./ std(dat6.swvl1)
ssr = (dat6.ssr .- mean(dat6.ssr)) ./ std(dat6.ssr)
tp = (dat6.tp .- mean(dat6.tp)) ./ std(dat6.tp)


data6 = DataFrame(species = keep, tmin = tmin, tmax = tmax, stl1 = stl1, swvl1 = swvl1, ssr = ssr, tp = tp)

# add the trait data to the tree
for i in eachrow(data6)
    setnodedata!(tree6, i.species, Phylo.traitdata(Union{Float64, ForwardDiff.Dual{ForwardDiff.Tag{DynamicPPL.DynamicPPLTag, Float64}, Float64, 3}}, ["tmin", "tmax", "stl1", "swvl1", "ssr", "tp"], [i.tmin, i.tmax, i.stl1, i.swvl1, i.ssr, i.tp]))
end

trait = ["tmin", "tmax", "stl1", "swvl1", "ssr", "tp"]
nodes = getnodes(tree6, postorder)

# add lengths to tree
for node in nodes
    val = getnodedata(tree6, node).value
    if hasinbound(tree6, node)
        len = Phylo.getlength(tree6, Phylo.getinbound(tree6, node))
        td = traitdata(eltype(nodedatatype(typeof(tree6))), trait, val, len)
        setnodedata!(tree6, node, td)
    else
        td = traitdata(eltype(nodedatatype(typeof(tree6))), trait, val)
        setnodedata!(tree6, node, td)
    end
end

zvec = [getnodedata(tree6, leaf).value for leaf in keep]
z6 = reduce(vcat, zvec)