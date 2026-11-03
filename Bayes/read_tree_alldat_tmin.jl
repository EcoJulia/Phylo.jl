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
const tree::TraitTree{1} = open(f -> parsenewick(f, TraitTree{1}),
                                   "Data/Qian2016.tree");

# remove missing species from dataframe
dropmissing!(df, :species)

# add underscores to dataframe
df.species = replace.(df.species, " " => "_")

# filter tree to just species we have data for
keep = intersect(getleafnames(tree), df.species)
keeptips!(tree, keep)                                  
filter!(:species => x -> x ∈ keep, df)

trait = ["tmin"]
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

# add leaf numbers to dataframe
leaves = getleafnames(tree, postorder)

df.leafnumber = 
indexin(df.species, leaves)

data = DataFrame(species = df.species, tmin = df.tmin, leafnumber = df.leafnumber)

data[data.leafnumber .== nothing, :]

