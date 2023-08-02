using BenchmarkTools
using Phylo
using DataFrames
using Random
using CSV

#needed for phylolm comaparison
using RCall

#needed for Phyloetworks comaparison
using PhyloNetworks
using StatsModels

Random.seed!(1234)

#Compare Phylo against phylolm and PhyloNetworks on a tree with uniform random tips
#Load the hummingbirds tree from test
const jtree::TraitTree{1} = open(f -> parsenewick(f, TraitTree{1}), "test/hummingbirds.tree")

#Create dataframe with leafnames and random trait values
species = getleaves(jtree)
data = 1000 .* rand(length(species))

#Create dataframe
speciesnames = getleafnames(jtree);
rdat = DataFrame(species = speciesnames, data = data);

#add the data to the tree
for i in eachrow(rdat)
    setnodedata!(jtree, i.species, Phylo.traitdata(["trait"], [i.data]));
end

#run inference
estimaterates(jtree, ["trait"])

#check inference run time
@benchmark estimaterates(jtree, ["trait"])

#run inference to calculate signal
estimaterates(jtree, ["trait"], lambda = 0.1)

#check inference run time
@benchmark estimaterates(jtree, ["trait"], lambda = 0.1)

#Compare against phylolm in R
#load R packages and read tree into R
R"
library('ape')
library('phylolm')"

R"rtree = ape::read.tree(file = 'test/hummingbirds.tree')" 

#send dataframe into R
@rput rdat

#species need to be the row names in R
R"row.names(rdat) <- rdat$species"

#run inference
rcopy(R"phylolm(data~1,data=rdat,phy=rtree)")

#check inference run time
@benchmark rcopy(R"phylolm(data~1,data=rdat,phy=rtree)")

#run inference to calculate signal
rcopy(R"phylolm(data~1,data=rdat,phy=rtree, model = c('lambda'))")

#check inference run time
@benchmark rcopy(R"phylolm(data~1,data=rdat,phy=rtree, model = c('lambda'))")

#Compare against PhyloNetworks
#load the tree using PhyloNetworks
tree2 = readTopology("test/hummingbirds.tree");

#get tip names
tipnames = tipLabels(tree2)

#dataframe needs column called tipNames
dat = DataFrame(tipNames = tipnames, data = data)

#run inference 
phylolm(@formula(data ~ 1), dat, tree2)

#check inference run time
@benchmark plm = phylolm(@formula(data ~ 1), dat, tree2)

#run inference to calculate signal
phylolm(@formula(data ~ 1), dat, tree2, model="lambda")

#check inference run time
@benchmark plm = phylolm(@formula(data ~ 1), dat, tree2, model="lambda")





#Comparison done using the same tree, but with tip trait data generated using Brownian Motion
#Generate data using Brownian Motion
a = BrownianTrait(jtree, "BMtrait")
BMtraits = rand(a)

#create a vector to store trait data
nodes = getnodenames(jtree)
leafnames = getleafnames(jtree);
traitvector = Vector{Float64}();

for leaf in leafnames
    push!(traitvector, BMtraits[leaf])
end

#Crate dataframe
rdat = DataFrame(species = leafnames, data = traitvector)

#add data to the tree
for i in eachrow(rdat)
    setnodedata!(jtree, i.species, Phylo.traitdata(["trait"], [i.data]));
end

#run inference
estimaterates(jtree, ["trait"])

#check inference run time
@benchmark estimaterates(jtree, ["trait"])

#run inference to calculate signal
estimaterates(jtree, ["trait"], lambda = 0.1)

#check inference run time
@benchmark estimaterates(jtree, ["trait"], lambda = 0.1)

#Compare against phylolm in R
#load R packages and read tree into R
R"
library('ape')
library('phylolm')"

R"rtree = ape::read.tree(file = 'test/hummingbirds.tree')" 

#send dataframe into R
@rput rdat

#species need to be the row names in R
R"row.names(rdat) <- rdat$species"

#run inference
rcopy(R"phylolm(data~1,data=rdat,phy=rtree)")

#check inference run time
@benchmark rcopy(R"phylolm(data~1,data=rdat,phy=rtree)")

#run inference to calculate signal
rcopy(R"phylolm(data~1,data=rdat,phy=rtree, model = c('lambda'))")

#check inference run time
@benchmark rcopy(R"phylolm(data~1,data=rdat,phy=rtree, model = c('lambda'))")

#Compare against PhyloNetworks
#dataframe needs column called tipNames
dat = DataFrame(tipNames = tipnames, data = traitvector)

#dataframe needs column called tipNames
dat = DataFrame(tipNames = tipnames, data = data)

#run inference 
phylolm(@formula(data ~ 1), dat, tree2)

#check inference run time
@benchmark plm = phylolm(@formula(data ~ 1), dat, tree2)

#run inference to calculate signal
phylolm(@formula(data ~ 1), dat, tree2, model="lambda")

#check inference run time
@benchmark plm = phylolm(@formula(data ~ 1), dat, tree2, model="lambda")



#Using the Myrtaceae tree, compare packages using real world data
#load the data
df = DataFrame(CSV.File("Data/Myrtaceae.csv"))

#load the tree
const bigtree::TraitTree{1} = open(f -> parsenewick(f, TraitTree{1}), "Data/Qian2016.tree")

#remove missing species from dataframe
dropmissing!(df, :species)

#add underscores to dataframe
df.species = replace.(df.species, " " => "_")

#filter the tree and dataframe for the species that are in  both
keep = intersect(getleafnames(bigtree), df.species)
keeptips!(bigtree, keep)
filter!(:species => x -> x ∈ keep, df)

#use mean data for trait value
gdf = groupby(df, :species)
dat = combine(gdf, [:tmin, :tmax, :trng, :stl1, :stl2, :stl3, :stl4, :swvl1, :swvl2, :swvl3, :swvl4, :ssr, :tp] .=> mean; renamecols=false)

#add the data for tmin to the tree
for i in eachrow(dat)
    setnodedata!(bigtree, i.species, Phylo.traitdata(["tmin"], [i.tmin]));
end

#run inference
estimaterates(bigtree, ["tmin"])

#check inference run time
@benchmark estimaterates(bigtree, ["tmin"])

#run inference to calculate signal
estimaterates(bigtree, ["tmin"], lambda = 0.8)

#check inference run time
@benchmark estimaterates(bigtree, ["tmin"], lambda = 0.8)

#Compare against phylolm in R
#load packages and tree
R"
library('ape')
library('phylolm')"

R"rbigtree = ape::read.tree(file = 'Data/Qian2016.tree')"

#transfer data and list of species in both the tree and data into R
@rput dat
@rput keep

#species need to be the row names in R
R"row.names(dat) <- dat$species"

#filter the tree
R"filteredtree <- keep.tip(rbigtree, keep)"

#run inference
rcopy(R"phylolm(tmin~1,data=dat,phy=filteredtree)") 

#check inference run time
@benchmark rcopy(R"phylolm(tmin~1,data=dat,phy=filteredtree)")

#run inference to calculate signal
rcopy(R"phylolm(tmin~1,data=dat,phy=filteredtree, model = c('lambda'))")

#check inference run time
@benchmark rcopy(R"phylolm(tmin~1,data=dat,phy=filteredtree, model = c('lambda'))") 

#Compare against PhyloNetworks
#load the tree and get the names of the leaves
tree = readTopology("Data/Qian2016.tree");
tipnames = tipLabels(tree)

# Filter for species in the tree and the data
missing_species = setdiff(tipnames, keep)
for i in eachindex(missing_species)
    deleteleaf!(tree, join(split(missing_species[i], " "), "_"))
end

#species column in the data needs to be tipNames
rename!(dat, :species => :tipNames)

#run inference
phylolm(@formula(tmin ~ 1), dat, tree)

#check inference run time
@benchmark phylolm(@formula(tmin ~ 1), dat, tree)

#run inference to calculate signal
phylolm(@formula(tmin ~ 1), dat, tree, model="lambda") 

#check run time
@benchmark phylolm(@formula(tmin ~ 1), dat, tree, model="lambda")


#Preform inference on multiple traits, both generated uniformly

#load the hummingbirds tree that can store data for two trait values
const multtree::TraitTree{2} = open(f -> parsenewick(f, TraitTree{2}), "test/hummingbirds.tree")

#load leafnames and generate uniform random trait values
species = getleaves(multtree)
data1 = 1000 .* rand(length(species))
data2 = 1000 .* rand(length(species))

#create dataframe
speciesnames = getleafnames(multtree);
rdat = DataFrame(species = speciesnames, data1 = data1, data2 = data2);

#add the data to the tree
for i in eachrow(rdat)
    setnodedata!(multtree, i.species, Phylo.traitdata(["trait1", "trait2"], [i.data1, i.data2]));
end

#run inference
estimaterates(multtree, ["trait1", "trait2"]) 

#check inference run time
@benchmark estimaterates(multtree, ["trait1", "trait2"])

#run inference to calculate Signal
estimaterates(multtree, ["trait1", "trait2"], lambda = 0.1)

#check inference run time
@benchmark estimaterates(multtree, ["trait1", "trait2"], lambda = 0.1)

#Preform inference on multiple traits, both generated through Brownian Motion
#load a new tree
const multtree2::TraitTree{2} = open(f -> parsenewick(f, TraitTree{2}), "test/hummingbirds.tree")

#gennerate data
a = BrownianTrait(multtree2, "BMtrait")
b = BrownianTrait(multtree2, "BMtrait")
BMtrait1 = rand(a)
BMtrait2 = rand(b)

#get leaf names
leafnames = getleafnames(multtree2);

#add data to a vector
traitvector1 = Vector{Float64}();
traitvector2 = Vector{Float64}();
for leaf in leafnames
    push!(traitvector1, BMtrait1[leaf])
    push!(traitvector2, BMtrait2[leaf])
end

#create a dataframe
rdat = DataFrame(species = leafnames, data1 = traitvector1, data2 = traitvector2)

#add data to the tree
for i in eachrow(rdat)
    setnodedata!(multtree2, i.species, Phylo.traitdata(["trait1", "trait2"], [i.data1, i.data2]));
end

#run inference
estimaterates(multtree2, ["trait1", "trait2"]) 

#check inference run time
@benchmark estimaterates(multtree2, ["trait1", "trait2"])

#run inference to calculate Signal
estimaterates(multtree2, ["trait1", "trait2"], lambda = 0.1)

#check inference run time
@benchmark estimaterates(multtree2, ["trait1", "trait2"], lambda = 0.1)

#Preform inference on multiple traits using real world data
#load the tree
const multtree3::TraitTree{2} = open(f -> parsenewick(f, TraitTree{2}), "Data/Qian2016.tree")

#load the data
df = DataFrame(CSV.File("Data/Myrtaceae.csv"))

#remove missing species from dataframe
dropmissing!(df, :species)

#add underscores to dataframe
df.species = replace.(df.species, " " => "_")

#filter the tree and dataframe for the species that are in  both
keep = intersect(getleafnames(multtree3), df.species)
keeptips!(multtree3, keep)
filter!(:species => x -> x ∈ keep, df)

#use mean data for trait value
gdf = groupby(df, :species)
dat = combine(gdf, [:tmin, :tmax, :trng, :stl1, :stl2, :stl3, :stl4, :swvl1, :swvl2, :swvl3, :swvl4, :ssr, :tp] .=> mean; renamecols=false)

#add the data for tmin to the tree
for i in eachrow(dat)
    setnodedata!(multtree3, i.species, Phylo.traitdata(["tmin", "tmax"], [i.tmin, i.tmax]));
end

#run inference
estimaterates(multtree3, ["tmin", "tmax"]) 

#check inference run time
@benchmark estimaterates(multtree3, ["tmin", "tmax"])

#run inference to calculate Signal
estimaterates(multtree3, ["tmin", "tmax"], lambda = 0.5)

#check inference run time
@benchmark estimaterates(multtree3, ["tmin", "tmax"], lambda = 0.5)


