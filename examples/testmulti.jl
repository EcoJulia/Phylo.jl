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

n_tips = 10
#generate tree
nu = Ultrametric{TraitTree{2}}(n_tips);
tree = rand(nu)

a = BrownianTrait(tree, "BMtrait")
b = BrownianTrait(tree, "BMtrait")
bm_traits1 = rand(a)
bm_traits2 = rand(b)

leafnames = getleafnames(tree);
z1 = [bm_traits1[leaf] for leaf in leafnames];
z2 = [bm_traits2[leaf] for leaf in leafnames];

nodes = getnodes(tree, postorder)

#Crate dataframe
rdat = DataFrame(species = leafnames, data1 = z1, data2 = z2)

#add data to the tree
for i in eachrow(rdat)
    setnodedata!(tree, i.species,
                 Phylo.traitdata(["trait1", "trait2"], [i.data1, i.data2]))
end

estimaterates(tree, ["trait1", "trait2"], lambda = [0.8, 0.8])

#Preform inference on multiple traits using real world data
#load the tree
const multtree3::TraitTree{2} = open(f -> parsenewick(f, TraitTree{2}),
                                     "Data/Qian2016.tree");

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

#add the data for tmin to the tree
for i in eachrow(dat)
    setnodedata!(multtree3, i.species,
                 Phylo.traitdata(["tmin", "tmax"],
                                 [i.tmin, i.tmax]))
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
=#

#test w/o signal to see if getting the thing I expect

const multtree2::TraitTree{2} = open(f -> parsenewick(f, TraitTree{2}),
                                     "test/hummingbirds.tree");

#gennerate data
a = BrownianTrait(multtree2, "BMtrait", σ² = 100)
b = BrownianTrait(multtree2, "BMtrait", σ² = 3)
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
rdat = DataFrame(species = leafnames, data1 = traitvector1,
                 data2 = traitvector2)

#add data to the tree
for i in eachrow(rdat)
    setnodedata!(multtree2, i.species,
                 Phylo.traitdata(["trait1", "trait2"], [i.data1, i.data2]))
end

#run inference
estimaterates(multtree2, ["trait1", "trait2"])

estimaterates(multtree2, ["trait1", "trait2"], lambda = [1.0, 1.0])

nodes = getnodes(multtree2, postorder)

nN = last(nodes)
nd = getnodedata(multtree2, nN)

k = 2

n = nleaves(multtree2)

betahat = inv(nd.xx) * nd.Q
sigmahat = (nd.yy - 2 * betahat * nd.Q' + betahat * nd.xx * betahat') / n
