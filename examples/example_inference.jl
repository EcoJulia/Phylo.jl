using BenchmarkTools
using Phylo
using DataFrames
using Random
using ProfileView

Random.seed!(1234)

jtree = open(parsenewick, "test/hummingbirds.tree")

#Create dataframe with leafnames and random trait values
species = getleaves(jtree)
data = 1000 .* rand(length(species))
dat = DataFrame(species = species, data = data)

@btime jfit = estimaterates(jtree, dat)

#Profile.clear_malloc_data()

ProfileView.@profview for i in 1:20 
    estimaterates(jtree, dat)
end
