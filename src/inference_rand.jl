# SPDX-License-Identifier: BSD-2-Clause

using LinearAlgebra
using Distributions
using ForwardDiff
using DynamicPPL
using ReverseDiff
using Bijectors

import Distributions: _logpdf, loglikelihood

# used for Bayes calculations 

# define loglikelihood function, used in Bayes methods
#=
function loglik(n, nd, sigma, beta)
    return -(1.0 / 2.0) * (n * log(2π) + nd.logV + n * log(abs(sigma)) +
            abs(sigma)^(-1) * (nd.yy[] - 2 * nd.Q[] * beta + nd.xx * beta^2))
end
=#
function loglik(n, nd, sigma2, beta)
    return -(1/2) * (
        n * log(2π) +
        nd.logV +
        n * log(sigma2) +
        (1/sigma2) * (nd.yy[] - 2*nd.Q[]*beta + nd.xx * beta^2)
    )
end

# Need to create own distribution to use threepoint to calculate likelihood
mutable struct MyDist2{T <: AbstractTree, N <: Number} <:
               ContinuousMultivariateDistribution
    sigma::N
    beta::Union{N, ReverseDiff.TrackedReal}
    tree::T
end

eltype(::MD) where {T, N <: Number, MD <: MyDist2{T, N}} = N
eltype(::Type{MD}) where {T, N <: Number, MD <: MyDist2{T, N}} = N

# rand creates a vector of tip trait values dependent on the tree, sigma (rate of evolution) and beta (root trait value)
function Distributions.rand(rng::AbstractRNG, d::MyDist2)
    a = BrownianTrait(d.tree, "BMtrait", d.beta; σ² = d.sigma)
    bm_traits = rand(a)
    z = [bm_traits[leaf] for leaf in getleafnames(d.tree, postorder)]
    return z
end

# define logpdf for my dist
function Distributions.logpdf(d::MyDist2, z::AbstractVector)#z::Vector{Float64})
    # add errors for if tree doesnt have right data

    nodes = getnodes(d.tree, postorder)
    trait = getnodedata(d.tree, nodes[1]).name
    leaves = getleaves(d.tree);
     # add tipdata to tree
    for i in 1:length(leaves)
        len = Phylo.getlength(d.tree, Phylo.getinbound(d.tree, leaves[i]))
        td = traitdata(eltype(nodedatatype(typeof(d.tree))), trait, [z[i]], len)
        setnodedata!(d.tree, leaves[i], td)
    end

    n = nleaves(d.tree)

    threepoint!(d.tree, trait, nodes)

    nN = last(nodes)
    nd = getnodedata(d.tree, nN)

    return loglik(n, nd, d.sigma, d.beta)
end

Distributions.length(d::MyDist2) = nleaves(d.tree)
Distributions.dim(d::MyDist2) = nleaves(d.tree)
Base.size(d::MyDist2) = (length(d),)
Base.size(d::MyDist2, i::Integer) = (i == 1 ? length(d) : 1)

Bijectors.bijector(d::MyDist2) = identity

#BrownianTrait that takes lambda into acount
struct BrownianTraitSignal{T <: AbstractTree, N <: Number} <:
       Sampleable{Univariate, EvolvedTrait{T}}
    tree::T
    trait::String
    start::N
    λ::N
    _rand::Function
    f::Function
end

function BrownianTraitSignal(tree::T, trait::String, start::N = 0.0, λ::N = 1.0;
                       σ² = missing, σ = missing,
                       f::Function = identity) where
         {T <: AbstractTree, N <: Number}
    iscontinuous(N) ||
        throw(TypeError(:BrownianTrait,
                        "Type $N must be continuous for a Gaussian Trait",
                        AbstractFloat, N))
    if ismissing(σ)
        σ = sqrt(σ²)
    end

    f ≢ identity && f(start) ≉ start &&
        @warn "Note that the third argument (the starting state) for trait '$trait' is untransformed - transformed version is $(f(start))"

    if ismissing(σ)
        dimension(N) ≡
        dimension(sqrt(getlength(tree, first(getbranches(tree))))) ||
            throw(DimensionMismatch("Dimensions of start, σ[²] and branch lengths must combine correctly if using Unitful"))
        return BrownianTraitSignal{T, N}(tree, trait, start, λ,
                                   ((rng::AbstractRNG, start::N, length) -> start +
                                                                            randn(rng,
                                                                                  typeof(one(start))) *
                                                                            sqrt(length)),
                                   f)
    else
        dimension(N) ≡
        dimension(σ * sqrt(getlength(tree,
                                     first(getbranches(tree))))) ||
            throw(DimensionMismatch("Dimensions of start, σ[²] and branch lengths must combine correctly if using Unitful"))
        return BrownianTraitSignal{T, N}(tree, trait, start, λ,
                                   ((rng::AbstractRNG, start::N, length) -> start +
                                                                            σ *
                                                                            randn(rng,
                                                                                  typeof(one(start))) *
                                                                            sqrt(length)),
                                   f)
    end
end

function rand!(rng::AbstractRNG,
               bm::BrownianTraitSignal{TREE, N},
               tree::TREE) where {TREE <: AbstractTree, N <: Number}
    trait = Dict{nodetype(TREE), N}()
    use_dict = (bm.f ≢ identity)
    for node in traversal(tree, preorder)
        if isroot(tree, node)
            if use_dict
                trait[node] = bm.start
                setnodedata!(tree, node, bm.trait, bm.f(bm.start))
            else
                setnodedata!(tree, node, bm.trait, bm.start)
            end
        elseif isleaf(tree, node)
            inb = getinbound(tree, node)
            prt = src(tree, inb)
            previous = use_dict ? trait[prt] :
                       getnodedata(tree, prt, bm.trait)
            h = getheight(tree, node)
            length = λ * N(getlength(tree, inb)) + (1 - λ) * h
            value = bm._rand(rng, previous, length) 
            if use_dict
                trait[node] = value
                setnodedata!(tree, node, bm.trait, bm.f(value))
            else
                setnodedata!(tree, node, bm.trait, value)
            end
        else
            inb = getinbound(tree, node)
            prt = src(tree, inb)
            previous = use_dict ? trait[prt] :
                       getnodedata(tree, prt, bm.trait)
            value = bm._rand(rng, previous, λ * N(getlength(tree, inb)))
            if use_dict
                trait[node] = value
                setnodedata!(tree, node, bm.trait, bm.f(value))
            else
                setnodedata!(tree, node, bm.trait, value)
            end
        end
    end
    return tree
end

function rand(rng::AbstractRNG,
              bm::BrownianTraitSignal{TREE, N}) where {TREE <: AbstractTree,
                                                 N <: Number}
    untrait = Dict{nodetype(TREE), N}()
    traitbyname = Dict{nodenametype(TREE), typeof(bm.f(bm.start))}()
    for node in traversal(bm.tree, preorder)
        if isroot(bm.tree, node)
            untrait[node] = bm.start
            traitbyname[getnodename(bm.tree, node)] = bm.f(bm.start)
        elseif isleaf(bm.tree, node)
            inb = getinbound(bm.tree, node)
            prt = src(bm.tree, inb)
            previous = untrait[prt]
            h = getheight(bm.tree, node)
            length = bm.λ * N(getlength(bm.tree, inb)) + (1 - bm.λ) * h
            value = bm._rand(rng, previous, length)
            untrait[node] = value
            traitbyname[getnodename(bm.tree, node)] = bm.f(value)
        else
            inb = getinbound(bm.tree, node)
            prt = src(bm.tree, inb)
            previous = untrait[prt]
            value = bm._rand(rng, previous, bm.λ * N(getlength(bm.tree, inb)))
            untrait[node] = value
            traitbyname[getnodename(bm.tree, node)] = bm.f(value)
        end
    end
    return traitbyname
end

# Distribution for Bayes using threepoint algorithm to find phylogenetic signal
struct MyDist3{T <: AbstractTree, N <: Number} <:
       ContinuousMultivariateDistribution
    sigma::N
    beta::N
    lambda::N
    tree::T
end

function eltype(::MD) where {T <: AbstractTree, N <: Number, MD <:
                                                             MyDist3{T, N}}
    return eltype(nodedatatype(T))
end
function eltype(::Type{MD}) where {T <: AbstractTree, N <: Number,
                                   MD <: MyDist3{T, N}}
    return eltype(nodedatatype(T))
end

# rand creates a vector of tip trait values dependent on the tree, sigma (rate of evolution) and beta (root trait value)
function Distributions.rand(rng::AbstractRNG, d::MyDist3) # incorrect but can fix later
    a = BrownianTraitSignal(d.tree, "BMtrait", start = d.beta, λ = d.lambda, σ² = d.sigma)
    bm_traits = rand(rng, a)

    z = [bm_traits[leaf] for leaf in getleafnames(d.tree, postorder)]
    return z
end

# define logpdf for my dist
function Distributions.logpdf(d::MD, z::Vector{Float64}) where {MD <: MyDist3}
    # add errors for if tree doesnt have right data

    n = nleaves(d.tree)
    nodes = getnodes(d.tree, postorder)
    trait = getnodedata(d.tree, nodes[1]).name

    # multiply internal branches by lambda
    for node in nodes
        if isleaf(d.tree, node)
            getnodedata(d.tree, node).t = d.lambda *
                                          getlength(d.tree,
                                                    getinbound(d.tree, node)) +
                                          (1.0 - d.lambda) *
                                          heighttoroot(d.tree, node)
        elseif isroot(d.tree, node)
            getnodedata(d.tree, node).t = zero(d.lambda)
        else
            getnodedata(d.tree, node).t = d.lambda *
                                          getlength(d.tree,
                                                    getinbound(d.tree, node))
        end
    end

    threepoint!(d.tree, trait, nodes)

    nN = last(nodes)
    nd = getnodedata(d.tree, nN)

    return loglik(n, nd, d.sigma, d.beta)
end

# Adapt BrownianTrait for multiple traits
struct BrownianTraitMult{T <: AbstractTree, N <: Number} <:
       Sampleable{Multivariate, EvolvedTrait{T}}
    tree::T
    trait::Vector{String}
    start::Vector{N}
    _rand::Function
    f::Function
end

eltype(::MD) where {T, N <: Number, MD <: BrownianTraitMult{T, N}} = N
eltype(::Type{MD}) where {T, N <: Number, MD <: BrownianTraitMult{T, N}} = N

function BrownianTraitMult(tree::T, trait::Vector{String}, start::Vector{N};
                       σ = missing,
                       f::Function = identity) where
         {T <: AbstractTree, N <: Number}
    iscontinuous(N) ||
        throw(TypeError(:BrownianTraitMult,
                        "Type $N must be continuous for a Gaussian Trait",
                        AbstractFloat, N))

    ntraits = length(trait)
    
    if ismissing(start)
        start = zeros(ntraits)
    end

    f ≢ identity && f(start) ≉ start &&
        @warn "Note that the third argument (the starting state) for trait '$trait' is untransformed - transformed version is $(f(start))"


    if ismissing(σ)
        dimension(N) ≡
        dimension(sqrt(getlength(tree, first(getbranches(tree))))) ||
            throw(DimensionMismatch("Dimensions of start, σ[²] and branch lengths must combine correctly if using Unitful"))
        return BrownianTraitMult{T, N}(tree, trait, start,
                                   ((rng::AbstractRNG, start::Vector{N}, length) -> start + randn(rng, N, size(start)) * length),
                                   f)
    else
        return BrownianTraitMult{T, N}(tree, trait, start,
                                   ((rng::AbstractRNG, start::Vector{N}, length) -> start +
                                                                            σ *
                                                                            randn(rng, N, size(start)) *
                                                                            length),
                                   f)
    end
    
end

function rand!(rng::AbstractRNG,
               bm::BrownianTraitMult{TREE, N},
               tree::TREE) where {TREE <: AbstractTree, N <: Number}
    trait = Dict{nodetype(TREE), N}()
    use_dict = (bm.f ≢ identity)
    for node in traversal(tree, preorder)
        if isroot(tree, node)
            if use_dict
                trait[node] = bm.start
                setnodedata!(tree, node, bm.trait, bm.f(bm.start))
            else
                setnodedata!(tree, node, bm.trait, bm.start)
            end
        else
            inb = getinbound(tree, node)
            prt = src(tree, inb)
            previous = use_dict ? trait[prt] :
                       getnodedata(tree, prt, bm.trait)
            value = bm._rand(rng, previous, N(getlength(tree, inb)))
            if use_dict
                trait[node] = value
                setnodedata!(tree, node, bm.trait, bm.f(value))
            else
                setnodedata!(tree, node, bm.trait, value)
            end
        end
    end
    return tree
end

function rand(rng::AbstractRNG,
              bm::BrownianTraitMult{TREE, N}) where {TREE <: AbstractTree,
                                                 N <: Number}
    untrait = Dict{nodetype(TREE), Vector{N}}()
    traitbyname = Dict{nodenametype(TREE), typeof(bm.f(bm.start))}()
    for node in traversal(bm.tree, preorder)
        if isroot(bm.tree, node)
            untrait[node] = bm.start
            traitbyname[getnodename(bm.tree, node)] = bm.f(bm.start)
        else
            inb = getinbound(bm.tree, node)
            prt = src(bm.tree, inb)
            previous = untrait[prt]
            value = bm._rand(rng, previous, N(getlength(bm.tree, inb)))
            untrait[node] = value
            traitbyname[getnodename(bm.tree, node)] = bm.f(value)
        end
    end
    return traitbyname
end




# Distribution for multiple traits

mutable struct MyDist4{T <: AbstractTree, N <: Number} <:
               ContinuousMultivariateDistribution
    sigma::Matrix{N} 
    beta::AbstractVector #Union{Vector{N}, ReverseDiff.TrackedArray}
    tree::T
end


eltype(::MD) where {T, N <: Number, MD <: MyDist4{T, N}} = N
eltype(::Type{MD}) where {T, N <: Number, MD <: MyDist4{T, N}} = N



function Distributions.rand(rng::AbstractRNG, d::MyDist4)
    traitnames = getnodedata(d.tree, getroot(d.tree)).name
    a = BrownianTraitMult(d.tree, traitnames, d.beta, σ = d.sigma)
    bm_traits = rand(a)
    z_vec = [bm_traits[leaf] for leaf in getleafnames(d.tree, postorder)]
    z = reduce(vcat, z_vec)
    return z
end

idx(tip::Int, trait::Int, n_traits::Int) = ((tip - 1) * n_traits) + trait

# define logpdf for my dist
function Distributions._logpdf(d::MD, z::AbstractArray) where {MD <: MyDist4}

    n = nleaves(d.tree)
    nodes = getnodes(d.tree, postorder)
    trait = getnodedata(d.tree, nodes[1]).name
    m = size(trait)[1]

    leaves = getleaves(d.tree);
     # add tipdata to tree
    for i in 1:length(leaves)
        # get the m means for tip i from vector z
        idx = (i-1)*m + 1 : i*m
        len = Phylo.getlength(d.tree, Phylo.getinbound(d.tree, leaves[i]))
        td = traitdata(eltype(nodedatatype(typeof(d.tree))), trait, z[idx], len)
        setnodedata!(d.tree, leaves[i], td)
    end

    threepoint!(d.tree, trait, nodes)

    nN = last(nodes)
    nd = getnodedata(d.tree, nN)

    return -(1.0 / 2.0) * (n * m * log(2π) + m * nd.logV + n * log(abs(det(d.sigma))) + tr((nd.yy .- 2 * nd.Q * d.beta' .+ d.beta * nd.xx * d.beta') * inv(d.sigma)))
end

loglikelihood(d::MyDist4, z::AbstractVector{<:Number}) = _logpdf(d, z)

Distributions.length(d::MyDist4) = nleaves(d.tree) * length(d.beta)

Distributions.dim(d::MyDist4) = nleaves(d.tree) * length(d.beta)


Base.size(d::MyDist4) = (length(d),)
Base.size(d::MyDist4, i::Integer) = (i == 1 ? length(d) : 1)

Bijectors.bijector(d::MyDist4) = identity