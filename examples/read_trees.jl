# SPDX-License-Identifier: BSD-2-Clause

using Phylo
using DataFrames
using LinearAlgebra
using CSV
using Statistics
using ForwardDiff
using DynamicPPL

# Read the tree, add data to tips and add lengths to branches

# Set up Myrtaceae trees

df = DataFrame(CSV.File("Data/Myrtaceae.csv"))

# load the tree for 1 trait
const bigtree::TraitTree{1} = open(f -> parsenewick(f, TraitTree{1}),
                                   "Data/Qian2016.tree")

# remove missing species from dataframe
dropmissing!(df, :species)

# add underscores to dataframe
df.species = replace.(df.species, " " => "_")

# filter the tree and dataframe for the species that are in  both
keep = intersect(getleafnames(bigtree), df.species)
keeptips!(bigtree, keep)
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

# add the data for tmin to the tree
for i in eachrow(dat)
    setnodedata!(bigtree, i.species, Phylo.traitdata(Union{Float64, ForwardDiff.Dual{ForwardDiff.Tag{DynamicPPL.DynamicPPLTag, Float64}, Float64, 3}}, ["tmin"], [i.tmin]))
end

# trait needs to be a vector of trait names, used for functions later
trait = ["tmin"]
nodes = getnodes(bigtree, postorder)

# add lengths to tree
for node in nodes
    val = getnodedata(bigtree, node).value
    if hasinbound(bigtree, node)
        len = Phylo.getlength(bigtree, Phylo.getinbound(bigtree, node))
        td = traitdata(eltype(nodedatatype(typeof(bigtree))), trait, val, len)
        setnodedata!(bigtree, node, td)
    else
        td = traitdata(eltype(nodedatatype(typeof(bigtree))), trait, val)
        setnodedata!(bigtree, node, td)
    end
end

# vector of trait values for the tips
z1 = dat.tmin

# 2 trait tree
# load the tree
const bigtree2::TraitTree{2} = open(f -> parsenewick(f, TraitTree{2}),
                                   "Data/Qian2016.tree")


# filter tree to just species we have data for
keep = intersect(getleafnames(bigtree2), df.species)
keeptips!(bigtree2, keep)

# add the data for tmin and tmax to the tree
for i in eachrow(dat)
    setnodedata!(bigtree2, i.species, Phylo.traitdata(Union{Float64, ForwardDiff.Dual{ForwardDiff.Tag{DynamicPPL.DynamicPPLTag, Float64}, Float64, 3}}, ["tmin", "tmax"], [i.tmin, i.tmax]))
end

# trait needs to be a vector of trait names, used for functions later
trait = ["tmin", "tmax"]
nodes = getnodes(bigtree2, postorder)

# add lengths to tree
for node in nodes
    val = getnodedata(bigtree2, node).value
    if hasinbound(bigtree2, node)
        len = Phylo.getlength(bigtree2, Phylo.getinbound(bigtree2, node))
        td = traitdata(eltype(nodedatatype(typeof(bigtree2))), trait, val, len)
        setnodedata!(bigtree2, node, td)
    else
        td = traitdata(eltype(nodedatatype(typeof(bigtree2))), trait, val)
        setnodedata!(bigtree2, node, td)
    end
end


leaves = getleafnames(bigtree2)
zvec = [getnodedata(bigtree2, leaf).value for leaf in leaves]
z2 = reduce(vcat, zvec)

# 3 trait tree
# load the tree
const bigtree3::TraitTree{3} = open(f -> parsenewick(f, TraitTree{3}),
                                   "Data/Qian2016.tree")


# filter tree to just species we have data for
keep = intersect(getleafnames(bigtree3), df.species)
keeptips!(bigtree3, keep)

# add the data for tmin and tmax to the tree
for i in eachrow(dat)
    setnodedata!(bigtree3, i.species, Phylo.traitdata(Union{Float64, ForwardDiff.Dual{ForwardDiff.Tag{DynamicPPL.DynamicPPLTag, Float64}, Float64, 3}}, ["tmin", "tmax", "trng"], [i.tmin, i.tmax, i.trng]))
end

# trait needs to be a vector of trait names, used for functions later
trait = ["tmin", "tmax", "trng"]
nodes = getnodes(bigtree3, postorder)

# add lengths to tree
for node in nodes
    val = getnodedata(bigtree3, node).value
    if hasinbound(bigtree3, node)
        len = Phylo.getlength(bigtree3, Phylo.getinbound(bigtree3, node))
        td = traitdata(eltype(nodedatatype(typeof(bigtree3))), trait, val, len)
        setnodedata!(bigtree3, node, td)
    else
        td = traitdata(eltype(nodedatatype(typeof(bigtree3))), trait, val)
        setnodedata!(bigtree3, node, td)
    end
end


leaves = getleafnames(bigtree3)
zvec = [getnodedata(bigtree3, leaf).value for leaf in leaves]
z3 = reduce(vcat, zvec)

# 4 trait tree
# load the tree
const bigtree4::TraitTree{4} = open(f -> parsenewick(f, TraitTree{4}),
                                   "Data/Qian2016.tree")


# filter tree to just species we have data for
keep = intersect(getleafnames(bigtree4), df.species)
keeptips!(bigtree4, keep)

# add the data for tmin and tmax to the tree
for i in eachrow(dat)
    setnodedata!(bigtree4, i.species, Phylo.traitdata(Union{Float64, ForwardDiff.Dual{ForwardDiff.Tag{DynamicPPL.DynamicPPLTag, Float64}, Float64, 3}}, ["tmin", "tmax", "trng", "stl1"], [i.tmin, i.tmax, i.trng, i.stl1]))
end

# trait needs to be a vector of trait names, used for functions later
trait = ["tmin", "tmax", "trng", "stl1"]
nodes = getnodes(bigtree4, postorder)

# add lengths to tree
for node in nodes
    val = getnodedata(bigtree4, node).value
    if hasinbound(bigtree4, node)
        len = Phylo.getlength(bigtree4, Phylo.getinbound(bigtree4, node))
        td = traitdata(eltype(nodedatatype(typeof(bigtree4))), trait, val, len)
        setnodedata!(bigtree4, node, td)
    else
        td = traitdata(eltype(nodedatatype(typeof(bigtree4))), trait, val)
        setnodedata!(bigtree4, node, td)
    end
end


leaves = getleafnames(bigtree4)
zvec = [getnodedata(bigtree4, leaf).value for leaf in leaves]
z4 = reduce(vcat, zvec)

# 5 trait tree
# load the tree
const bigtree5::TraitTree{5} = open(f -> parsenewick(f, TraitTree{5}),
                                   "Data/Qian2016.tree")


# filter tree to just species we have data for
keep = intersect(getleafnames(bigtree5), df.species)
keeptips!(bigtree5, keep)

# add the data for tmin and tmax to the tree
for i in eachrow(dat)
    setnodedata!(bigtree5, i.species, Phylo.traitdata(Union{Float64, ForwardDiff.Dual{ForwardDiff.Tag{DynamicPPL.DynamicPPLTag, Float64}, Float64, 3}}, ["tmin", "tmax", "trng", "stl1", "stl2"], [i.tmin, i.tmax, i.trng, i.stl1, i.stl2]))
end

# trait needs to be a vector of trait names, used for functions later
trait = ["tmin", "tmax", "trng", "stl1", "stl2"]
nodes = getnodes(bigtree5, postorder)

# add lengths to tree
for node in nodes
    val = getnodedata(bigtree5, node).value
    if hasinbound(bigtree5, node)
        len = Phylo.getlength(bigtree5, Phylo.getinbound(bigtree5, node))
        td = traitdata(eltype(nodedatatype(typeof(bigtree5))), trait, val, len)
        setnodedata!(bigtree5, node, td)
    else
        td = traitdata(eltype(nodedatatype(typeof(bigtree5))), trait, val)
        setnodedata!(bigtree5, node, td)
    end
end


leaves = getleafnames(bigtree5)
zvec = [getnodedata(bigtree5, leaf).value for leaf in leaves]
z5 = reduce(vcat, zvec)

# All trait tree
# load the tree
const bigtree13::TraitTree{13} = open(f -> parsenewick(f, TraitTree{13}),
                                   "Data/Qian2016.tree")


# filter tree to just species we have data for
keep = intersect(getleafnames(bigtree13), df.species)
keeptips!(bigtree13, keep)

# add the data for tmin and tmax to the tree
for i in eachrow(dat)
    setnodedata!(bigtree13, i.species, Phylo.traitdata(Union{Float64, ForwardDiff.Dual{ForwardDiff.Tag{DynamicPPL.DynamicPPLTag, Float64}, Float64, 3}}, ["tmin", "tmax", "trng", "stl1", "stl2", "stl3", "stl4", "swvl1", "swvl2", "swvl3", "swvl4", "ssr", "tp"], [i.tmin, i.tmax, i.trng, i.stl1, i.stl2, i.stl3, i.stl4, i.swvl1, i.swvl2, i.swvl3, i.swvl4, i.ssr, i.tp]))
end

# trait needs to be a vector of trait names, used for functions later
trait = ["tmin", "tmax", "trng", "stl1", "stl2", "stl3", "stl4", "swvl1", "swvl2", "swvl3", "swvl4", "ssr", "tp"]
nodes = getnodes(bigtree13, postorder)

# add lengths to tree
for node in nodes
    val = getnodedata(bigtree13, node).value
    if hasinbound(bigtree13, node)
        len = Phylo.getlength(bigtree13, Phylo.getinbound(bigtree13, node))
        td = traitdata(eltype(nodedatatype(typeof(bigtree13))), trait, val, len)
        setnodedata!(bigtree13, node, td)
    else
        td = traitdata(eltype(nodedatatype(typeof(bigtree13))), trait, val)
        setnodedata!(bigtree13, node, td)
    end
end


leaves = getleafnames(bigtree13)
zvec = [getnodedata(bigtree13, leaf).value for leaf in leaves]
z13 = reduce(vcat, zvec)

# get C
n_tips = nleaves(bigtree)
height = nodeheights(bigtree; onlyleaves = true)[1]
C = fill(height, (n_tips, n_tips)) - distances(bigtree) ./ 2
C = abs.(Symmetric(C))
