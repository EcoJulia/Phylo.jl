# SPDX-License-Identifier: BSD-2-Clause

using Turing
using Phylo
using Distributions
using Random
using BenchmarkTools
using DataFrames
using LinearAlgebra
using Plots
using CSV
using ForwardDiff
using DynamicPPL

# set seed
Random.seed!(678)

@model function βσ()
    β ~ Uniform(-100, 1000)
    σ ~ Uniform(0, 100)
    return β, σ
end

@model function βσ_covariance(z, C)
    @submodel β, σ = βσ()
    z ~ MvNormal(β * ones(length(z)), σ * C)
    return nothing
end

@model function βσ_threepoint(tree, z) # z needs to be for leaves in postorder
    @submodel β, σ = βσ()
    z ~ Phylo.MyDist2(σ, β, tree) # tree.z ~ (implement later)
    return nothing
end

@model function βσλ_threepoint(tree, z, upper = 1.0) # z needs to be for leaves in postorder
    @submodel β, σ = βσ()
    λ ~ Uniform(0, 1.0)
    z ~ Phylo.MyDist3(σ, β, λ, tree)
    return nothing
end

function gen_βσ_covariance(::Type{T}, n_tips) where {T}
    # generate tree
    nu = Ultrametric{T}(n_tips)
    tree = rand(nu)

    # get phylogenetic variance matrix - needed for method using built in Julia functions
    C = fill(1.0, (n_tips, n_tips)) - distances(tree) ./ 2
    C = abs.(Symmetric(C))

    # generate traits, store in vector z
    a = BrownianTrait(tree, "BMtrait", σ² = 1.0)
    bm_traits = rand(a)
    leafnames = getleafnames(tree)
    z = [bm_traits[leaf] for leaf in leafnames]

    # using built in Julia functions 
    return βσ_covariance(z, C)
end

function gen_βσ_threepoint(::Type{T}, n_tips) where {T}
    # generate tree
    nu = Ultrametric{T}(n_tips)
    tree = rand(nu)
    # generate traits, store in vector z
    a = BrownianTrait(tree, "BMtrait", σ² = 1.0)
    bm_traits = rand(a)
    leafnames = getleafnames(tree)
    z = [bm_traits[leaf] for leaf in leafnames]

    # Use threepoint:
    # pop trait data on tree
    nodes = getnodes(tree, postorder)
    rdat = DataFrame(species = leafnames, data = z)
    for i in eachrow(rdat)
        setnodedata!(tree, i.species,
                     Phylo.traitdata(eltype(nodedatatype(typeof(tree))),
                                     ["trait"],
                                     [i.data]))
    end

    # trait needs to be a vector of trait names, used for functions later
    trait = ["trait"]

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

    return βσ_threepoint(tree, z)
end

function gen_βσλ_threepoint(::Type{T}, n_tips) where {T}
    # generate tree
    nu = Ultrametric{T}(n_tips)
    tree = rand(nu)
    # distances(tree)

    # generate traits, store in vector z
    a = BrownianTrait(tree, "BMtrait", σ² = 1.0)
    bm_traits = rand(a)
    leafnames = getleafnames(tree)
    z = [bm_traits[leaf] for leaf in leafnames]

    # Use threepoint:
    # pop trait data on tree
    nodes = getnodes(tree, postorder)
    rdat = DataFrame(species = leafnames, data = z)
    for i in eachrow(rdat)
        setnodedata!(tree, i.species,
                     Phylo.traitdata(eltype(nodedatatype(typeof(tree))),
                                     ["trait"], [i.data]))
    end

    # trait needs to be a vector of trait names, used for functions later
    trait = ["trait"]

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

    return βσλ_threepoint(tree, z)
end

# number of tips on the tree
n_tips = 200
n_samples = 10_000

basemodel = gen_βσ_covariance(TraitTree{1}, n_tips);
spl = sample(basemodel, HMC(0.01, 5), n_samples) # add initial_params
plot(spl[:β])
plot(spl[:σ])

tpmodel = gen_βσ_threepoint(Phylo.TraitTreeNum{1}, n_tips);
tpmodel = gen_βσ_threepoint(TraitTree{1}, n_tips);
tpmodel = gen_βσ_threepoint(Phylo.TraitTreeFloat64{1}, n_tips);
spl = sample(tpmodel, HMC(0.01, 5), n_samples) # add initial_params
plot(spl[:β])
plot(spl[:σ])

tplmodel = gen_βσλ_threepoint(Phylo.TraitTreeNum{1}, n_tips);
tplmodel = gen_βσλ_threepoint(Phylo.TraitTreeDual{1}, n_tips);
tplmodel = gen_βσλ_threepoint(TraitTree{1}, n_tips);
spl = sample(tplmodel, HMC(0.01, 5), n_samples) # add initial_params
plot(spl[:β])
plot(spl[:σ])
plot(spl[:λ])

# real data
# Using the Myrtaceae tree, compare packages using real world data
# load the data
df = DataFrame(CSV.File("Data/Myrtaceae.csv"))

# load the tree
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

#check data in right order
n_samples = 1_000
model1 = βσ_threepoint(bigtree, dat.tmin);
spl1 = sample(model1, HMC(0.01, 5), n_samples)#; initial_params = [290, 8])
# 29.01 seconds
estimaterates(bigtree, ["tmin"])
loglikelihood(model1, (β=290, σ=8))

# get C
n_tips = nleaves(bigtree)
height = nodeheights(bigtree; onlyleaves = true)[1]
C = fill(height, (n_tips, n_tips)) - distances(bigtree) ./ 2
C = abs.(Symmetric(C))
model2 = βσ_covariance(dat.tmin, C);
spl2 = sample(model2, HMC(0.01, 5), n_samples)
# 279.69 seconds
loglikelihood(model2, (β=290, σ=8))

#data w/ lambda
model3 = βσλ_threepoint(bigtree, dat.tmin);
spl3 = sample(model3, HMC(0.01, 5), n_samples)
# 41.93 seconds
estimaterates(bigtree, ["tmin"], lambda = 0.5)



#Testing for scaling
n_samples = 10_000

#10 tips
tpmodel1 = gen_βσ_threepoint(TraitTree{1}, 10);
spl1 = sample(tpmodel1, HMC(0.01, 5), n_samples)
#9.53s

#100 tips
tpmodel2 = gen_βσ_threepoint(TraitTree{1}, 100);
spl2 = sample(tpmodel2, HMC(0.01, 5), n_samples)
#41.94s

#250 tips
tpmodel3 = gen_βσ_threepoint(TraitTree{1}, 250);
spl3 = sample(tpmodel3, HMC(0.01, 5), n_samples)
#104.46s

#500 tips
tpmodel4 = gen_βσ_threepoint(TraitTree{1}, 500);
spl4 = sample(tpmodel4, HMC(0.01, 5), n_samples)
#215.17s

#750 tips
tpmodel5 = gen_βσ_threepoint(TraitTree{1}, 750);
spl5 = sample(tpmodel5, HMC(0.01, 5), n_samples)
#327.25s

#1,000 tips
tpmodel6 = gen_βσ_threepoint(TraitTree{1}, 1_000);
spl6 = sample(tpmodel6, HMC(0.01, 5), n_samples)
#445.84s

#2,000 tips
tpmodel7 = gen_βσ_threepoint(TraitTree{1}, 2_000);
spl7 = sample(tpmodel7, HMC(0.01, 5), n_samples)
#984.99s

#10,000 tips
tpmodel8 = gen_βσ_threepoint(TraitTree{1}, 10_000);
spl8 = sample(tpmodel8, HMC(0.01, 5), n_samples)
#6191.94s



# with signal 

n_samples = 10_000

#10 tips
tpmodel1 = gen_βσλ_threepoint(TraitTree{1}, 10);
spl1 = sample(tpmodel1, HMC(0.01, 5), n_samples)
# 12.83 seconds

#100 tips
tpmodel2 = gen_βσλ_threepoint(TraitTree{1}, 100);
spl2 = sample(tpmodel2, HMC(0.01, 5), n_samples)
# 50.39 seconds

#250 tips
tpmodel3 = gen_βσλ_threepoint(TraitTree{1}, 250);
spl3 = sample(tpmodel3, HMC(0.01, 5), n_samples)
# 134.47 seconds

#500 tips
tpmodel4 = gen_βσλ_threepoint(TraitTree{1}, 500);
spl4 = sample(tpmodel4, HMC(0.01, 5), n_samples)
# 274.9 seconds

#750 tips
tpmodel5 = gen_βσλ_threepoint(TraitTree{1}, 750);
spl5 = sample(tpmodel5, HMC(0.01, 5), n_samples)
# 438.18 seconds

#1,000 tips
tpmodel6 = gen_βσλ_threepoint(TraitTree{1}, 1_000);
spl6 = sample(tpmodel6, HMC(0.01, 5), n_samples)
# 601.1 seconds
