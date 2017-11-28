import Distributions: ValueSupport, Sampleable
import Base: eltype, rand
using Phylo
using Distributions

mutable struct Phylogenetics{T <: AbstractTree} <: ValueSupport end
Base.eltype(::Type{Phylogenetics{T}}) where T <: AbstractTree = T

"""
    Nonultrametric{T <: AbstractTree,
                   RNG <: Sampleable}(n::Int,
                                      rng::RNG = Exponential())
    Nonultrametric{T <: AbstractTree,
                   RNG <: Sampleable}(tiplabels::Vector{String},
                                      rng::RNG = Exponential())

The sampler for non-ultrametric phylogenetic trees of size `n` or with
tip labels `tiplabels`. Generate random trees by calling rand().
Currently only works for `NamedTree`s.
"""
struct Nonultrametric{T <: AbstractTree,
                              RNG <: Sampleable{Univariate, Continuous}} <:
    Sampleable{Univariate, Phylogenetics{T}}
    n::Int
    tiplabels::Vector{String}
    rng::RNG
    
    function Nonultrametric{T, RNG}(n::Int, rng::RNG) where {T, RNG}
        return new{T, RNG}(n, map(i -> "tip $i", 1:n), rng)
    end
    
    function Nonultrametric{T, RNG}(tiplabels::Vector{String}, rng::RNG) where {T, RNG}
        return new{T, RNG}(length(tiplabels), tiplabels, rng)
    end
end

function Nonultrametric{T}(n::Int) where T <: AbstractTree
    return Nonultrametric{T, Exponential}(n, Exponential())
end

Nonultrametric(n::Int) = Nonultrametric{NamedTree}(n)

function Nonultrametric{T}(tiplabels::Vector{String}) where T <: AbstractTree
    return Nonultrametric{T, Exponential}(tiplabels, Exponential())
end

Nonultrametric(tiplabels::Vector{String}) = Nonultrametric{NamedTree}(tiplabels)

function rand(t::Nonultrametric{T, RNG}) where {T, RNG}
    t.n >= 2 || error("A tree must have at least 2 tips")
    tree = T(t.tiplabels)
    roots = nodenamefilter(isroot, tree)
    while length(roots) > 1
        children = sample(collect(roots), 2, replace=false)
        parent = addnode!(tree)
        addbranch!(tree, parent, children[1], rand(t.rng))
        addbranch!(tree, parent, children[2], rand(t.rng))
        roots = nodenamefilter(isroot, tree)
    end
    return tree
end

"""
    Ultrametric{T <: AbstractTree,
                RNG <: Sampleable}(n::Int,
                                   rng::RNG = Exponential())
    Ultrametric{T <: AbstractTree,
                RNG <: Sampleable}(tiplabels::Vector{String},
                                   rng::RNG = Exponential())

The sampler for ultrametric phylogenetic trees of size `n` or with
tip labels `tiplabels`. Generate random trees by calling rand().
Currently only works for `NamedTree`s.
"""
struct Ultrametric{T <: AbstractTree,
                 RNG <: Sampleable{Univariate, Continuous}} <:
    Sampleable{Univariate, Phylogenetics{T}}
    n::Int
    tiplabels::Vector{String}
    rng::RNG
    
    function Ultrametric{T, RNG}(n::Int, rng::RNG) where {T, RNG}
        return new{T, RNG}(n, map(i -> "tip $i", 1:n), rng)
    end
    
    function Ultrametric{T, RNG}(tiplabels::Vector{String}, rng::RNG) where {T, RNG}
        return new{T, RNG}(length(tiplabels), tiplabels, rng)
    end
end

function Ultrametric{T}(n::Int) where T <: AbstractTree
    return Ultrametric{T, Exponential}(n, Exponential())
end

Ultrametric(n::Int) = Ultrametric{NamedTree}(n)

function Ultrametric{T}(tiplabels::Vector{String}) where T <: AbstractTree
    return Ultrametric{T, Exponential}(tiplabels, Exponential())
end

Ultrametric(tiplabels::Vector{String}) = Ultrametric{NamedTree}(tiplabels)

function rand(t::Ultrametric{T, RNG}) where {T, RNG}
    t.n >= 2 || error("A tree must have at least 2 tips")
    tree = T(t.tiplabels)
    roots = nodenamefilter(isroot, tree)
    depth = zero(rand(t.rng))
    leaves = getleafnames(tree)
    while length(roots) > 1
        tocoalesce = collect(roots)
        coalescers = sample(tocoalesce, 2, replace=false)
        parent = addnode!(tree)
        depth += rand(t.rng) * 2.0 / length(tocoalesce)
        c1 = filter(x -> coalescers[1] ∈ nodehistory(tree, x), leaves)
        d1 = map(x -> getheight(tree, x), c1)
        c2 = filter(x -> coalescers[2] ∈ nodehistory(tree, x), leaves)
        d2 = map(x -> getheight(tree, x), c2)
        addbranch!(tree, parent, coalescers[1], depth - d1[1])
        addbranch!(tree, parent, coalescers[2], depth - d2[1])
        roots = nodenamefilter(isroot, tree)
    end
    return tree
end
