using Compat.Random
import Distributions: ValueSupport, Sampleable
import Base: eltype, rand
using Phylo
using Distributions
using Missings
using IterableTables: getiterator

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
    leafinfo::Any

    function Nonultrametric{T, RNG}(n::Int, tiplabels::Vector{String}, rng::RNG, leafinfo) where {T, RNG}
        return new{T, RNG}(n, tiplabels, rng, leafinfo)
    end

    function Nonultrametric{T, RNG}(n::Int, rng::RNG) where {T, RNG}
        return new{T, RNG}(n, ["tip $i" for i in Base.OneTo(n)], rng, missing)
    end

    function Nonultrametric{T, RNG}(tiplabels::Vector{String}, rng::RNG) where {T, RNG}
        return new{T, RNG}(length(tiplabels), tiplabels, rng, missing)
    end
end

function Nonultrametric{T}(n::Int) where T <: AbstractTree
    return Nonultrametric{T, Exponential}(n, Exponential())
end

function Nonultrametric{T}(tiplabels::Vector{String}) where T <: AbstractTree
    return Nonultrametric{T, Exponential}(tiplabels, Exponential())
end

function Nonultrametric{T}(leafinfo) where T <: AbstractTree
    tipnames = unique(collect(keys(leafinfo)))
    return Nonultrametric{T, Exponential}(length(tipnames), tipnames,
                                          Exponential(), leafinfo)
end

Nonultrametric(info::LI) where LI =
    Nonultrametric{NamedTree}(info)

Nonultrametric(n::Int) = Nonultrametric{NamedTree}(n)
Nonultrametric(tiplabels::Vector{String}) =
    Nonultrametric{NamedTree}(tiplabels)

function rand(t::Nonultrametric{T, RNG}) where {T, RNG}
    t.n >= 2 || error("A tree must have at least 2 tips")
    if ismissing(t.leafinfo)
        tree = T(t.tiplabels; rootheight = 0.0)
    else
        tree = T(t.leafinfo; rootheight = 0.0)
    end
    while nroots(tree) > 1
        roots = getroots(tree)
        children = sample(collect(roots), 2, replace=false)
        parent = createnode!(tree)
        createbranch!(tree, parent, children[1], rand(t.rng))
        createbranch!(tree, parent, children[2], rand(t.rng))
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
    leafinfo::Any

    function Ultrametric{T, RNG}(n::Int, tiplabels::Vector{String}, rng::RNG, leafinfo) where {T, RNG}
        return new{T, RNG}(n, tiplabels, rng, leafinfo)
    end

    function Ultrametric{T, RNG}(n::Int, rng::RNG) where {T, RNG}
        return new{T, RNG}(n, ["tip $i" for i in Base.OneTo(n)], rng, missing)
    end

    function Ultrametric{T, RNG}(tiplabels::Vector{String}, rng::RNG) where {T, RNG}
        return new{T, RNG}(length(tiplabels), tiplabels, rng, missing)
    end
end

function Ultrametric{T}(n::Int) where T <: AbstractTree
    return Ultrametric{T, Exponential}(n, Exponential())
end

function Ultrametric{T}(tiplabels::Vector{String}) where T <: AbstractTree
    return Ultrametric{T, Exponential}(tiplabels, Exponential())
end

function Ultrametric{T}(leafinfo) where T <: AbstractTree
    tipnames = unique(collect(keys(leafinfo)))
    return Ultrametric{T, Exponential}(length(tipnames), tipnames,
                                       Exponential(), leafinfo)
end

Ultrametric(info::LI) where LI = Ultrametric{NamedTree}(info)

Ultrametric(n::Int) = Ultrametric{NamedTree}(n)
Ultrametric(tiplabels::Vector{String}) = Ultrametric{NamedTree}(tiplabels)

function rand(t::Ultrametric{T, RNG}) where {T, RNG}
    t.n >= 2 || error("A tree must have at least 2 tips")
    if ismissing(t.leafinfo)
        tree = T(t.tiplabels; rootheight = 0.0)
    else
        tree = T(t.leafinfo; rootheight = 0.0)
    end
    depth = zero(rand(t.rng))
    leaves = getleaves(tree)
    while nroots(tree) > 1
        roots = getroots(tree)
        tocoalesce = collect(roots)
        coalescers = sample(tocoalesce, 2, replace=false)
        parent = createnode!(tree)
        depth += rand(t.rng) * 2.0 / length(tocoalesce)
        d1 = getheight(tree, first(nodefuture(tree, coalescers[1]) ∩
                                   leaves))
        d2 = getheight(tree, first(nodefuture(tree, coalescers[2]) ∩
                                   leaves))
        createbranch!(tree, parent, coalescers[1], depth - d1)
        createbranch!(tree, parent, coalescers[2], depth - d2)
    end
    return tree
end

function rand(s::S, treenames) where
    {TREE <: AbstractTree,
     S <: Sampleable{Univariate, Phylogenetics{TREE}}}
    trees = Dict{eltype(treenames), TREE}()
    for name in treenames
        trees[name] = rand(s)
    end
    return TreeSet(trees, Dict{String, Dict{String, Any}}())
end

function rand(s::S, n::Int) where
    {TREE <: AbstractTree,
     S <: Sampleable{Univariate, Phylogenetics{TREE}}}
    return rand(s, 1:n)
end
