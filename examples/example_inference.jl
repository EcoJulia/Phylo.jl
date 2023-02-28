using BenchmarkTools
using Phylo
using DataFrames
using Random

Random.seed!(1234)

jtree = open(parsenewick, "test/hummingbirds.tree")

#Create dataframe with leafnames and random trait values
species = getleafnames(jtree)
data = 1000 .* rand(length(species))
dat = DataFrame(species = species, data = data)

@btime jfit = estimaterates(jtree, dat)