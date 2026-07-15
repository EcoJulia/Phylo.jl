using BenchmarkTools
using Phylo
using DataFrames
using Random
using CSV
using Statistics
using ForwardDiff
using DynamicPPL
using LinearAlgebra

# set seed
Random.seed!(10) #problem when 1234

# set number of tips of the tree
n_tips = 10

# generate a tree
nu = Ultrametric{TraitTree{1}}(n_tips)
tree1 = rand(nu)

# generate trait values using brownian motion, should have signal
a = BrownianTrait(tree1, "BMtrait")
bm_traits = rand(a)

# create a vector to store trait data
nodes = getnodenames(tree1)
leafnames = getleafnames(tree1);
traitvector = [bm_traits[leaf] for leaf in leafnames]

# Crate dataframe
rdat = DataFrame(species = leafnames, data = traitvector)

# add data to the tree
for i in eachrow(rdat)
    setnodedata!(tree1, i.species, Phylo.traitdata(Union{Float64, ForwardDiff.Dual{ForwardDiff.Tag{DynamicPPL.DynamicPPLTag, Float64}, Float64, 3}}, 
                    ["trait"], [i.data]))
end


# set number of tips of the tree
n_tips = 100

# generate a tree
nu = Ultrametric{TraitTree{1}}(n_tips)
tree2 = rand(nu)

# generate trait values using brownian motion, should have signal
a = BrownianTrait(tree2, "BMtrait")
bm_traits = rand(a)

# create a vector to store trait data
nodes = getnodenames(tree2)
leafnames = getleafnames(tree2);
traitvector = [bm_traits[leaf] for leaf in leafnames]

# Crate dataframe
rdat = DataFrame(species = leafnames, data = traitvector)

# add data to the tree
for i in eachrow(rdat)
    setnodedata!(tree2, i.species, Phylo.traitdata(Union{Float64, ForwardDiff.Dual{ForwardDiff.Tag{DynamicPPL.DynamicPPLTag, Float64}, Float64, 3}}, 
                    ["trait"], [i.data]))
end


# set number of tips of the tree
n_tips = 1_000

# generate a tree
nu = Ultrametric{TraitTree{1}}(n_tips)
tree3 = rand(nu)

# generate trait values using brownian motion, should have signal
a = BrownianTrait(tree3, "BMtrait")
bm_traits = rand(a)

# create a vector to store trait data
nodes = getnodenames(tree3)
leafnames = getleafnames(tree3);
traitvector = [bm_traits[leaf] for leaf in leafnames]

# Crate dataframe
rdat = DataFrame(species = leafnames, data = traitvector)

# add data to the tree
for i in eachrow(rdat)
    setnodedata!(tree3, i.species, Phylo.traitdata(Union{Float64, ForwardDiff.Dual{ForwardDiff.Tag{DynamicPPL.DynamicPPLTag, Float64}, Float64, 3}}, 
                    ["trait"], [i.data]))
end


# set number of tips of the tree
n_tips = 10_000

# generate a tree
nu = Ultrametric{TraitTree{1}}(n_tips)
tree4 = rand(nu)

# generate trait values using brownian motion, should have signal
a = BrownianTrait(tree4, "BMtrait")
bm_traits = rand(a)

# create a vector to store trait data
nodes = getnodenames(tree4)
leafnames = getleafnames(tree4);
traitvector = [bm_traits[leaf] for leaf in leafnames]

# Crate dataframe
rdat = DataFrame(species = leafnames, data = traitvector)

# add data to the tree
for i in eachrow(rdat)
    setnodedata!(tree4, i.species, Phylo.traitdata(Union{Float64, ForwardDiff.Dual{ForwardDiff.Tag{DynamicPPL.DynamicPPLTag, Float64}, Float64, 3}}, 
                    ["trait"], [i.data]))
end


# set number of tips of the tree
n_tips = 50_000

# generate a tree
nu = Ultrametric{TraitTree{1}}(n_tips)
tree5 = rand(nu)

# generate trait values using brownian motion, should have signal
a = BrownianTrait(tree5, "BMtrait")
bm_traits = rand(a)

# create a vector to store trait data
nodes = getnodenames(tree5)
leafnames = getleafnames(tree5);
traitvector = [bm_traits[leaf] for leaf in leafnames]

# Crate dataframe
rdat = DataFrame(species = leafnames, data = traitvector)

# add data to the tree
for i in eachrow(rdat)
    setnodedata!(tree5, i.species, Phylo.traitdata(Union{Float64, ForwardDiff.Dual{ForwardDiff.Tag{DynamicPPL.DynamicPPLTag, Float64}, Float64, 3}}, 
                    ["trait"], [i.data]))
end

# set number of tips of the tree
n_tips = 100_000

# generate a tree
nu = Ultrametric{TraitTree{1}}(n_tips)
tree6 = rand(nu)

# generate trait values using brownian motion, should have signal
a = BrownianTrait(tree6, "BMtrait")
bm_traits = rand(a)

# create a vector to store trait data
nodes = getnodenames(tree6)
leafnames = getleafnames(tree6);
traitvector = [bm_traits[leaf] for leaf in leafnames]

# Crate dataframe
rdat = DataFrame(species = leafnames, data = traitvector)

# add data to the tree
for i in eachrow(rdat)
    setnodedata!(tree6, i.species, Phylo.traitdata(Union{Float64, ForwardDiff.Dual{ForwardDiff.Tag{DynamicPPL.DynamicPPLTag, Float64}, Float64, 3}}, 
                    ["trait"], [i.data]))
end

# run inference

# for 10 tips
@benchmark estimaterates(tree1, ["trait"])
# 85.559 μs
@benchmark estimaterates(tree1, ["trait"], lambda = 0.5)
# 12.475 ms

# for 100 tips
@benchmark estimaterates(tree2, ["trait"])
# 878.305 μs
@benchmark estimaterates(tree2, ["trait"], lambda = 0.5)
# 150.723 s

# for 1,000 tips
@benchmark estimaterates(tree3, ["trait"])
# 9.647 ms
@benchmark estimaterates(tree3, ["trait"], lambda = 0.5)
# 68.983 s

# for 10,000 tips
@benchmark estimaterates(tree4, ["trait"])
# 137.764 ms
@benchmark estimaterates(tree4, ["trait"], lambda = 0.5)
# 606.131 s

# for 50,000 tips
@benchmark estimaterates(tree5, ["trait"])
# 845.504 ms

# for 100,000 tips
@benchmark estimaterates(tree6, ["trait"])
# 1.746 s


# Scalinng for multi trait

# set number of tips of the tree
n_tips = 10_000


# 1 trait
# generate a tree
nu = Ultrametric{TraitTree{1}}(n_tips)
multtree1 = rand(nu)

# generate trait values using brownian motion, should have signal
a = BrownianTrait(multtree1, "BMtrait")
bm_traits = rand(a)

# create a vector to store trait data
nodes = getnodenames(multtree1)
leafnames = getleafnames(multtree1);
traitvector = [bm_traits[leaf] for leaf in leafnames]

# Crate dataframe
rdat = DataFrame(species = leafnames, data = traitvector)

# add data to the tree
for i in eachrow(rdat)
    setnodedata!(multtree1, i.species, Phylo.traitdata(Union{Float64, ForwardDiff.Dual{ForwardDiff.Tag{DynamicPPL.DynamicPPLTag, Float64}, Float64, 3}}, 
                    ["trait"], [i.data]))
end

# 2 traits
# generate a tree
nu = Ultrametric{TraitTree{2}}(n_tips)
multtree2 = rand(nu)

# generate trait values using brownian motion, should have signal
a = BrownianTrait(multtree2, "BMtrait")
b = BrownianTrait(multtree2, "BMtrait")
bm_trait1 = rand(a)
bm_trait2 = rand(b)

# create a vector to store trait data
nodes = getnodenames(multtree2)
leafnames = getleafnames(multtree2);
traitvector1 = [bm_trait1[leaf] for leaf in leafnames]
traitvector2 = [bm_trait2[leaf] for leaf in leafnames]

# Crate dataframe
rdat = DataFrame(species = leafnames, data1 = traitvector1, data2 = traitvector2)

# add data to the tree
for i in eachrow(rdat)
    setnodedata!(multtree2, i.species, Phylo.traitdata(Union{Float64, ForwardDiff.Dual{ForwardDiff.Tag{DynamicPPL.DynamicPPLTag, Float64}, Float64, 3}}, 
                    ["trait1", "trait2"], [i.data1, i.data2],))
end

# 3 traits
# generate a tree
nu = Ultrametric{TraitTree{3}}(n_tips)
multtree3 = rand(nu)

# generate trait values using brownian motion, should have signal
a = BrownianTrait(multtree3, "BMtrait")
b = BrownianTrait(multtree3, "BMtrait")
c = BrownianTrait(multtree3, "BMtrait")
bm_trait1 = rand(a)
bm_trait2 = rand(b)
bm_trait3 = rand(c)

# create a vector to store trait data
nodes = getnodenames(multtree3)
leafnames = getleafnames(multtree3);
traitvector1 = [bm_trait1[leaf] for leaf in leafnames]
traitvector2 = [bm_trait2[leaf] for leaf in leafnames]
traitvector3 = [bm_trait3[leaf] for leaf in leafnames]

# Crate dataframe
rdat = DataFrame(species = leafnames, data1 = traitvector1, data2 = traitvector2, data3 = traitvector3)

# add data to the tree
for i in eachrow(rdat)
    setnodedata!(multtree3, i.species, Phylo.traitdata(Union{Float64, ForwardDiff.Dual{ForwardDiff.Tag{DynamicPPL.DynamicPPLTag, Float64}, Float64, 3}}, 
                    ["trait1", "trait2", "trait3"], [i.data1, i.data2, i.data3],))
end

# for 1 trait
@benchmark estimaterates(multtree1, ["trait"])
# 128.392 ms

# for 2 traits
@benchmark estimaterates(multtree2, ["trait1", "trait2"])
# 147.194 ms

# for 3 traits
@benchmark estimaterates(multtree3, ["trait1", "trait2", "trait3"])
# 154.390 ms