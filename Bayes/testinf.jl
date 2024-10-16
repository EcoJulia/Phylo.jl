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

# generate a random tree
# number of tips on the tree
n_tips = 100

# generate tree
nu = Ultrametric{TraitTree{1}}(n_tips);
tree = rand(nu);
distances(tree)

# get phylogentic variance matrix - needed for method using built in Julia functions
C = fill(1.0, (n_tips, n_tips)) - distances(tree) ./ 2
C = abs.(Symmetric(C))

# generate traits, store in vector z
a = BrownianTrait(tree, "BMtrait", σ² = 1.0);
bm_traits = rand(a)

leafnames = getleafnames(tree);
z = Vector{Float64}();

for leaf in leafnames
    push!(z, bm_traits[leaf])
end

# using built in Julia functions
@model function phyloinf(z, C)
    beta ~ Uniform(-100, 100)
    sigma ~ Uniform(0, 100)

    return z ~ MvNormal(beta * ones(length(z)), sigma * C)
end

c = sample(phyloinf(z, C), HMC(0.1, 5), 100000)

plot(c[:beta])
plot(c[:sigma])

# Use threepoint:
# pop trait data on tree
nodes = getnodes(tree, postorder)
rdat = DataFrame(species = leafnames, data = z)
for i in eachrow(rdat)
    setnodedata!(tree, i.species, Phylo.traitdata(["trait"], [i.data]))
end

# trait needs to be a vector of trait names, used for functions later
trait = ["trait"]

# add lengths to tree
for node in nodes
    val = getnodedata(tree, node).value
    if hasinbound(tree, node)
        len = Phylo.getlength(tree, Phylo.getinbound(tree, node))
        td = traitdata(trait, val, len)
        setnodedata!(tree, node, td)
    else
        td = traitdata(trait, val)
        setnodedata!(tree, node, td)
    end
end

# Need to create own distribution to use threepoint to calculate likelihood
# needs to be a mutable struct as threepoint changes tree
struct MyDist2{T <: AbstractTree, N <: Number} <:
       ContinuousMultivariateDistribution
    sigma::N
    beta::N
    tree::T
end

# rand creates a vector of tip trait values dependent on the tree, sigma (rate of evolution) and beta (root trait value)
function Distributions.rand(rng::AbstractRNG, d::MyDist2)
    a = BrownianTrait(d.tree, "BMtrait", σ² = d.sigma)
    bm_traits = rand(a)

    leafnames = getleafnames(d.tree, postorder)
    z = Vector{Float64}()

    for leaf in leafnames
        push!(z, bm_traits[leaf])
    end
    return z
end

# define loglikelihood function, n is number of leaves, nd is node data from the root node, sigma is the rate of evolution, beta is the trait data for the root node
function loglik(n, nd, sigma, beta)
    return -(1.0 / 2.0) * (n * log(2π) + nd.logV + n * log(abs(sigma)) +
            abs(sigma)^(-1) * (nd.yy[] - 2 * nd.Q[] * beta + nd.xx * beta^2))
end

# define logpdf for my dist
function Distributions.logpdf(d::MyDist2, z::Vector{Float64})
    # add errors for if tree doesnt have right data

    n = nleaves(d.tree)
    nodes = getnodes(d.tree, postorder)
    trait = getnodedata(d.tree, nodes[1]).name

    threepoint!(d.tree, trait, nodes)

    nN = last(nodes)
    nd = getnodedata(tree, nN)

    return loglik(n, nd, d.sigma, d.beta)
end

@model function phyloinftree(tree, z) # z needs to be for leaves in postorder
    beta ~ Uniform(-100, 100)
    sigma ~ Uniform(0, 100)

    return z ~ MyDist2(sigma, beta, tree) # tree.z ~ (impliment later)
end

c2 = sample(phyloinftree(tree, z), HMC(0.01, 5), 100000)
plot(c2[:beta])
plot(c2[:sigma])

# define loglikelihood function, used in Bayes methods
function loglik(n, nd, sigma, beta)
    return -(1.0 / 2.0) * (n * log(2π) + nd.logV + n * log(abs(sigma)) +
            abs(sigma)^(-1) * (nd.yy[] - 2 * nd.Q[] * beta + nd.xx * beta^2))
end

# Need to create own distribution to use threepoint to calculate likelihood
# needs to be a mutable struct as threepoint changes tree
# also needs renamed
struct MyDist3{T <: AbstractTree, N <: Number} <:
       ContinuousMultivariateDistribution
    sigma::N
    beta::N
    lambda::N
    tree::T
end

# rand creates a vector of tip trait values dependent on the tree, sigma (rate of evolution) and beta (root trait value)
function Distributions.rand(rng::AbstractRNG, d::MyDist3) # incorrect but can fix later
    a = BrownianTrait(d.tree, "BMtrait", σ² = d.sigma)
    bm_traits = rand(a)

    leafnames = getleafnames(d.tree, postorder)
    z = Vector{Float64}()

    for leaf in leafnames
        push!(z, bm_traits[leaf])
    end
    return z
end

# define logpdf for my dist
function Distributions.logpdf(d::MyDist3, z::Vector{Float64})
    # add errors for if tree doesnt have right data

    n = nleaves(d.tree)
    nodes = getnodes(d.tree, postorder)
    trait = getnodedata(d.tree, nodes[1]).name

    # add lengths to tree - must be a better way
    for node in nodes
        val = getnodedata(tree, node).value
        if hasinbound(tree, node)
            len = Phylo.getlength(tree, Phylo.getinbound(tree, node))
            td = traitdata(trait, val, len)
            setnodedata!(tree, node, td)
        else
            td = traitdata(trait, val)
            setnodedata!(tree, node, td)
        end
    end

    # multiply internal branches by lambda
    t = [getnodedata(d.tree, node).t for node in nodes]

    for (i, node) in enumerate(nodes)
        if !isleaf(d.tree, node)
            tupdate = d.lambda * t[i]
            getnodedata(d.tree, node).t = tupdate
        end
    end

    threepoint!(d.tree, trait, nodes)

    nN = last(nodes)
    nd = getnodedata(d.tree, nN)

    return loglik(n, nd, d.sigma, d.beta)
end

@model function phyloinftreelambda(tree, z) # z needs to be for leaves in postorder
    intnodeheights = nodeheights(tree, noleaves = true)
    longnodeheight = maximum(intnodeheights)

    leafnodeheights = nodeheights(tree, onlyleaves = true)
    shortleafheight = minimum(leafnodeheights)

    upper = shortleafheight / longnodeheight

    lambda ~ Uniform(0, upper)
    beta ~ Uniform(-100, 100)
    sigma ~ Uniform(0, 100)

    return z ~ MyDist3(sigma, beta, lambda, tree) # tree.z ~ (impliment later)
end

c3 = sample(phyloinftreelambda(tree, z), HMC(0.01, 5), 10000) # add initial_params 

loglikelihood(phyloinftreelambda(tree, z), c3)
loglikelihood(phyloinftreelambda(tree, z),
              (lambda = 1.0, beta = 0.0, sigma = 1.0))
# is the problem when I set lambda to 0?

n = nleaves(tree)
nodes = getnodes(tree, postorder)
trait = getnodedata(tree, nodes[1]).name

threepoint!(tree, trait, nodes)

nN = last(nodes)
nd = getnodedata(tree, nN)

ndat = getnodedata(tree, nodes[5])
ndat.t

loglik(n, nd, 1.0, 0.0)

plot(c3[:beta])
plot(c3[:sigma])
plot(c3[:lambda])

# Other method - will likely delete but keep for now incase any of it needed 
#=
struct MyDist <: ContinuousMultivariateDistribution 
    sigma
    beta
end

# Distributions.rand(rng::AbstractRNG, d::MyDist) = # create random tree w/ beta sigma

function Distributions.logpdf(d::MyDist, tree::TraitTree) 
    # add errors for if tree doesnt have right data

    nodes = getnodes(tree, postorder)
    trait = getnodedata(tree, nodes[1]).name

    threepoint!(tree, trait, nodes)

    nN = last(nodes)
    nd = getnodedata(tree, nN)

    return loglik(n, nd, d.sigma, d.beta) # my defined - may need to have threepoint in function
end

function Distributions.loglikelihood(d::MyDist, tree::TraitTree) 
    # add errors for if tree doesnt have right data

    nodes = getnodes(tree, postorder)
    trait = getnodedata(tree, nodes[1]).name

    threepoint!(tree, trait, nodes)

    nN = last(nodes)
    nd = getnodedata(tree, nN)

    return loglik(n, nd, d.sigma, d.beta) # my defined - may need to have threepoint in function
end

@model function phyloinftree(tree)
    beta ~ Uniform(-100, 100) 
    sigma ~ Uniform(0, 100)

    tree ~ MyDist(sigma, beta)
    # return beta, sigma
end

c2 = sample(phyloinftree(tree), HMC(0.01, 5), 100000)
plot(c2[:beta])
plot(c2[:sigma])

=#
