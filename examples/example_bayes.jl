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
using PDMats


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

@model function βσ_threepoint(tree, z) 
    @submodel β, σ = βσ()
    z ~ Phylo.MyDist2(σ, β, tree) 
    return nothing
end

@model function βσλ_threepoint(tree, z, upper = 1.0) 
    @submodel β, σ = βσ()
    λ ~ Uniform(0.0, 1.0)
    z ~ Phylo.MyDist3(σ, β, λ, tree)
    return nothing
end

@model function βσ_multthreepoint(tree, z) 
    @submodel β, σ = βσ_mult()
    z ~ Phylo.MyDist4(σ, β, tree) 
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
