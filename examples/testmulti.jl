using Phylo
using DataFrames
using Random
using CSV
using Statistics
using Profile
using ProfileView
using Optim #need to remember to remove

estimaterates(multtree2, ["trait1", "trait2"], lambda = [0.8, 0.8]) 

Random.seed!(110)

x=10
#generate tree
nu = Ultrametric{TraitTree{2}}(x);
tree = rand(nu)

a = BrownianTrait(tree, "BMtrait")
b = BrownianTrait(tree, "BMtrait")
BMtraits1 = rand(a)
BMtraits2 = rand(b)

leafnames = getleafnames(tree);
z1 = Vector{Float64}();
z2 = Vector{Float64}();

for leaf in leafnames
    push!(z1, BMtraits1[leaf])
    push!(z2, BMtraits2[leaf])
end

nodes = getnodes(tree, postorder)

#Crate dataframe
rdat = DataFrame(species = leafnames, data1 = z1, data2 = z2)

#add data to the tree
for i in eachrow(rdat)
    setnodedata!(tree, i.species, Phylo.traitdata(["trait1", "trait2"], [i.data1, i.data2]));
end

estimaterates(tree, ["trait1", "trait2"], lambda = [0.8, 0.8]) 



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
    setnodedata!(multtree3, i.species, Phylo.traitdata(["tmin", "tmax"],
                                                         [i.tmin, i.tmax]));
end

lambda = fill(0.8, 2);

@time estimaterates(multtree3, ["tmin", "tmax"], lambda = lambda)

Profile.clear_malloc_data()

ProfileView.@profview for i in 1:3 
    estimaterates(multtree3, ["tmin", "tmax"], lambda = lambda)
end

#=
nodes = getnodes(multtree3, postorder)

nN = last(nodes)
nd = getnodedata(multtree3, nN)
    
k = length(2)

n = nleaves(multtree3)
