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

@model function phyloinftree(tree, z) #z needs to be for leaves in postorder
    β ~ Uniform(-100, 100)
    σ ~ Uniform(0, 100)

    return z ~ Phylo.MyDist2(σ, β, tree) #tree.z ~ (impliment later)
end

@model function phyloinftreelambda(tree, z, upper = 1.0) #z needs to be for leaves in postorder
    lambda ~ Uniform(0, 1.0)
    beta ~ Uniform(-100, 100)
    sigma ~ Uniform(0, 100)

    return z ~ Phylo.MyDist3(sigma, beta, lambda, tree)
end

function run_βσ_inference(::Type{T}, n_tips, n_samples) where {T}
    # generate tree
    nu = Ultrametric{T}(n_tips)
    tree = rand(nu)
    # distances(tree)

    # get phylogenetic variance matrix - needed for method using built in Julia functions
    # C = fill(1.0, (n_tips,n_tips)) - distances(tree)./2
    # C = abs.(Symmetric(C)) 

    # generate traits, store in vector z
    a = BrownianTrait(tree, "BMtrait", σ² = 1.0)
    bm_traits = rand(a)
    leafnames = getleafnames(tree)
    z = [bm_traits[leaf] for leaf in leafnames]

    # using built in Julia functions 
    #=
    @model function phyloinf(z,C)
        beta ~ Uniform(-100, 100) 
        sigma ~ Uniform(0, 100)

        z ~ MvNormal(beta * ones(length(z)), sigma * C) 
    end

    c = sample(phyloinf(z, C), HMC(0.1, 5), 100000)

    plot(c[:beta])
    plot(c[:sigma])

    =#

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

    return sample(phyloinftree(tree, z), HMC(0.01, 5), n_samples)
end

function run_βσλ_inference(::Type{T}, n_tips, n_samples) where {T}
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
                                     ["trait"],
                                     [i.data]))
    end

    #trait needs to be a vector of trait names, used for functions later
    trait = ["trait"]

    #add lengths to tree
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

    intnodeheights = nodeheights(tree, noleaves = true)
    longnodeheight = maximum(intnodeheights)

    leafnodeheights = nodeheights(tree, onlyleaves = true)
    shortleafheight = minimum(leafnodeheights)

    # calculate upper bound for signal
    upper = shortleafheight / longnodeheight

    return sample(phyloinftreelambda(tree, z), HMC(0.01, 5), n_samples) #add initial_params 
end

# number of tips on the tree
n_tips = 400
n_samples = 10_000
# spl = run_βσ_inference(TraitTree{1}, n_tips, n_samples)
spl = run_βσ_inference(Phylo.TraitTreeFloat64{1}, n_tips, n_samples)
plot(spl[:β])
plot(spl[:σ])

# spl = run_βσλ_inference(Phylo.TraitTreeNum{1}, n_tips, n_samples)
spl = run_βσλ_inference(TraitTree{1}, n_tips, n_samples)
plot(spl[:β])
plot(spl[:σ])
plot(spl[:λ])
