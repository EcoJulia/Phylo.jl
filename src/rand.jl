using Compat.Random
import Distributions: ValueSupport, Sampleable
import Distributions: _rand!
import Base: eltype, rand
using Phylo
using Distributions
using Missings
using IterableTables: getiterator

struct Phylogenetics{T <: AbstractTree} <: ValueSupport end
Base.eltype(::Type{Phylogenetics{T}}) where T <: AbstractTree = T
Base.eltype(::Sampleable{S, P}) where {S, P <: Phylogenetics} = eltype(P)

"""
    Nonultrametric{T <: AbstractTree,
                   SAMP <: Sampleable}(n::Int,
                                       sampleable::SAMP = Exponential())
    Nonultrametric{T <: AbstractTree,
                   SAMP <: Sampleable}(tiplabels::Vector{String},
                                       sampleable::SAMP = Exponential())

The sampler for non-ultrametric phylogenetic trees of size `n` or with
tip labels `tiplabels`. Generate random trees by calling rand().
"""
struct Nonultrametric{T <: AbstractTree,
                      SAMP <: Sampleable{Univariate, Continuous}} <:
                          Sampleable{Univariate, Phylogenetics{T}}
    n::Int
    tiplabels::Vector{String}
    sampleable::SAMP
    leafinfo::Any

    function Nonultrametric{T, SAMP}(n::Int, tiplabels::Vector{String},
                                     sampleable::SAMP, leafinfo) where {T, SAMP}
        return new{T, SAMP}(n, tiplabels, sampler(sampleable), leafinfo)
    end

    function Nonultrametric{T, SAMP}(n::Int, sampleable::SAMP) where {T, SAMP}
        return new{T, SAMP}(n, ["tip $i" for i in Base.OneTo(n)],
                            sampler(sampleable), missing)
    end

    function Nonultrametric{T, SAMP}(tiplabels::Vector{String},
                                     sampleable::SAMP) where {T, SAMP}
        return new{T, SAMP}(length(tiplabels), tiplabels,
                            sampler(sampleable), missing)
    end
end

function Nonultrametric{T}(n::Int) where T <: AbstractTree
    return Nonultrametric{T, Exponential}(n, Exponential())
end

function Nonultrametric{T}(tiplabels::Vector{String}) where T <: AbstractTree
    return Nonultrametric{T, Exponential}(tiplabels, Exponential())
end

function Nonultrametric{T}(leafinfo) where T <: AbstractTree
    tipnames = unique([info[1] for info in getiterator(leafinfo)])
    return Nonultrametric{T, Exponential}(length(tipnames), tipnames,
                                          Exponential(), leafinfo)
end

Nonultrametric(info::LI) where LI =
    Nonultrametric{Phylo.LT{OneRoot, LI}}(info)

Nonultrametric(n::Int) = Nonultrametric{RootedTree}(n)
Nonultrametric(tiplabels::Vector{String}) =
    Nonultrametric{RootedTree}(tiplabels)

function _rand!(rng::AbstractRNG, t::Nonultrametric{T, SAMP}) where {T, SAMP}
    t.n >= 2 || error("A tree must have at least 2 tips")
    if ismissing(t.leafinfo)
        tree = T(t.tiplabels)
    else
        tree = T(t.leafinfo)
    end
    while nroots(tree) > 1
        roots = getroots(tree)
        children = sample(collect(roots), 2, replace=false)
        parent = createnode!(tree)
        createbranch!(tree, parent, children[1], rand(rng, t.sampleable))
        createbranch!(tree, parent, children[2], rand(rng, t.sampleable))
    end
    return tree
end

"""
    Ultrametric{T <: AbstractTree,
                SAMP <: Sampleable}(n::Int,
                                    sampleable::SAMP = Exponential())
    Ultrametric{T <: AbstractTree,
                SAMP <: Sampleable}(tiplabels::Vector{String},
                                    sampleable::SAMP = Exponential())

The sampler for ultrametric phylogenetic trees of size `n` or with
tip labels `tiplabels`. Generate random trees by calling rand().
"""
struct Ultrametric{T <: AbstractTree,
                   SAMP <: Sampleable{Univariate, Continuous}} <:
                       Sampleable{Univariate, Phylogenetics{T}}
    n::Int
    tiplabels::Vector{String}
    sampleable::SAMP
    leafinfo::Any

    function Ultrametric{T, SAMP}(n::Int, tiplabels::Vector{String},
                                  sampleable::SAMP, leafinfo) where {T, SAMP}
        return new{T, SAMP}(n, tiplabels, sampler(sampleable), leafinfo)
    end

    function Ultrametric{T, SAMP}(n::Int, sampleable::SAMP) where {T, SAMP}
        return new{T, SAMP}(n, ["tip $i" for i in Base.OneTo(n)],
                            sampler(sampleable), missing)
    end

    function Ultrametric{T, SAMP}(tiplabels::Vector{String},
                                  sampleable::SAMP) where {T, SAMP}
        return new{T, SAMP}(length(tiplabels), tiplabels,
                            sampler(sampleable), missing)
    end
end

function Ultrametric{T}(n::Int) where T <: AbstractTree
    return Ultrametric{T, Exponential}(n, Exponential())
end

function Ultrametric{T}(tiplabels::Vector{String}) where T <: AbstractTree
    return Ultrametric{T, Exponential}(tiplabels, Exponential())
end

function Ultrametric{T}(leafinfo) where T <: AbstractTree
    tipnames = unique([info[1] for info in getiterator(leafinfo)])
    return Ultrametric{T, Exponential}(length(tipnames), tipnames,
                                       Exponential(), leafinfo)
end

Ultrametric(info::LI) where LI =
    Ultrametric{Phylo.LT{OneRoot, LI}}(info)

Ultrametric(n::Int) = Ultrametric{RootedTree}(n)
Ultrametric(tiplabels::Vector{String}) = Ultrametric{RootedTree}(tiplabels)

function _rand!(rng::AbstractRNG, t::Ultrametric{T, SAMP}) where {T, SAMP}
    t.n >= 2 || error("A tree must have at least 2 tips")
    if ismissing(t.leafinfo)
        tree = T(t.tiplabels)
    else
        tree = T(t.leafinfo)
    end
    depth = zero(rand(rng, t.sampleable))
    leaves = getleaves(tree)
    while nroots(tree) > 1
        roots = getroots(tree)
        tocoalesce = collect(roots)
        coalescers = sample(tocoalesce, 2, replace=false)
        parent = createnode!(tree)
        depth += rand(rng, t.sampleable) * 2.0 / length(tocoalesce)
        d1 = getheight(tree, first(nodefuture(tree, coalescers[1]) ∩
                                   leaves))
        d2 = getheight(tree, first(nodefuture(tree, coalescers[2]) ∩
                                   leaves))
        createbranch!(tree, parent, coalescers[1], depth - d1)
        createbranch!(tree, parent, coalescers[2], depth - d2)
    end
    return tree
end

rand(s::S, treenames::AbstractVector{LABEL}) where
{LABEL, TREE <: AbstractTree,
 S <: Sampleable{Univariate, Phylogenetics{TREE}}} =
     rand(Random.GLOBAL_RNG, s, treenames)

function rand(rng::AbstractRNG, s::S, treenames::AbstractVector{LABEL}) where
    {LABEL, TREE <: AbstractTree{OneTree},
     S <: Sampleable{Univariate, Phylogenetics{TREE}}}
    trees = Dict{eltype(treenames), TREE}()
    for name in treenames
        trees[name] = rand(rng, s)
    end
    return trees
end
