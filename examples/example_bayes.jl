# SPDX-License-Identifier: BSD-2-Clause

using Turing
using Phylo
using Distributions
using Random
using BenchmarkTools
using DataFrames
using LinearAlgebra
using Plots

# set seed
Random.seed!(678)

@model function βσ()
    β ~ Uniform(-100, 100)
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

model = gen_βσ_covariance(TraitTree{1}, n_tips);
sample(model, HMC(0.01, 5), 1);
spl = sample(model, HMC(0.01, 5), n_samples) # add initial_params
plot(spl[:β])
plot(spl[:σ])

model = gen_βσ_threepoint(Phylo.TraitTreeNum{1}, n_tips);
model = gen_βσ_threepoint(TraitTree{1}, n_tips);
model = gen_βσ_threepoint(Phylo.TraitTreeFloat64{1}, n_tips);
sample(model, HMC(0.01, 5), 1);
spl = sample(model, HMC(0.01, 5), n_samples) # add initial_params
plot(spl[:β])
plot(spl[:σ])

model = gen_βσλ_threepoint(Phylo.TraitTreeNum{1}, n_tips);
model = gen_βσλ_threepoint(Phylo.TraitTreeDual{1}, n_tips);
model = gen_βσλ_threepoint(TraitTree{1}, n_tips);
sample(model, HMC(0.01, 5), 1);
spl = sample(model, HMC(0.01, 5), n_samples) # add initial_params
plot(spl[:β])
plot(spl[:σ])
plot(spl[:λ])
