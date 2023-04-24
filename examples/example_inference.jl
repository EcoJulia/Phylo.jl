using BenchmarkTools
using Phylo
using DataFrames
using Random
using RCall

#needed for phylonetworks comapre
using PhyloNetworks
using StatsModels

Random.seed!(1234)

const jtree::TraitTree = open(f -> parsenewick(f, TraitTree), "test/hummingbirds.tree")

#Create dataframe with leafnames and random trait values
species = getleaves(jtree)
data = 1000 .* rand(length(species))
#dat = DataFrame(species = species, data = data)

R"
library('ape')
library('phylolm')"

R"rtree = ape::read.tree(file = 'test/hummingbirds.tree')" 

speciesnames = getleafnames(jtree);
rdat = DataFrame(species = speciesnames, data = data);
setnodedata!.(jtree, rdat.species, Phylo.traitdata.("trait", rdat.data));
@rput rdat

#species need to be the row names in R
R"row.names(rdat) <- rdat$species"

@btime rfit = rcopy(R"phylolm(data~1,data=rdat,phy=rtree)") #2.290 ms
@benchmark rfit = rcopy(R"phylolm(data~1,data=rdat,phy=rtree)")

@btime rfit2 = rcopy(R"phylolm(data~1,data=rdat,phy=rtree, model = c('lambda'))") #13ms
@benchmark rfit2 = rcopy(R"phylolm(data~1,data=rdat,phy=rtree, model = c('lambda'))") 

@btime jfit = estimaterates(jtree, "trait") #2 ms
@benchmark jfit = estimaterates(jtree, "trait")

@btime jfit2 = estimaterates(jtree, "trait", 0.1) #150 ms
@benchmark jfit2 = estimaterates(jtree, "trait", 0.1)

#=
using Profile
using ProfileView
Profile.clear_malloc_data()

ProfileView.@profview for i in 1:100
    estimaterates(jtree, "trait", 0.1)
end
=#

#############################################################
#Compare to phylo networks results 

#load the tree using phylonetworks
tree2 = readTopology("test/hummingbirds.tree");

#get tip names
tipnames = tipLabels(tree2)

#dataframe needs column called tipNames
dat = DataFrame(tipNames = tipnames, data = data)

@btime plm = phylolm(@formula(data ~ 1), dat, tree2) #1 ms also provides different values
@benchmark plm = phylolm(@formula(data ~ 1), dat, tree2)
