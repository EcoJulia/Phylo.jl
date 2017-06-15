import Distributions: ValueSupport, Sampleable
import Base: eltype, rand
using Phylo
using Distributions

type Phylogenetics{T <: AbstractTree} <: ValueSupport end
Base.eltype{T <: AbstractTree}(::Type{Phylogenetics{T}}) = T

"""
    Nonultrametric{T <: AbstractTree,
                   RNG <: Sampleable}(n::Int,
                                      treetype::Type{T} = NamedTree,
                                      rng::RNG = Exponential())
    Nonultrametric{T <: AbstractTree,
                   RNG <: Sampleable}(tiplabels::Vector{String},
                                      treetype::Type{T} = NamedTree,
                                      rng::RNG = Exponential())

The sampler for non-ultrametric phylogenetic trees of size `n` or with
tip labels `tiplabels`. Generate random trees by calling rand().
Currently only works for `NamedTree`s.
"""
type Nonultrametric{T <: AbstractTree,
                    RNG <: Sampleable{Univariate, Continuous}} <:
    Sampleable{Univariate, Phylogenetics{T}}
    n::Int
    tiplabels::Vector{String}
    rng::RNG
    
    function (::Type{Nonultrametric{T, RNG}}){T, RNG}(n::Int,
                                                      ::Type{T},
                                                      rng::RNG)
        return new{T, RNG}(n, map(i -> "tip $i", 1:n), rng)
    end
    
    function (::Type{Nonultrametric{T, RNG}}){T, RNG}(tiplabels::Vector{String},
                                                      ::Type{T},
                                                      rng::RNG)
        return new{T, RNG}(length(tiplabels), tiplabels, rng)
    end
end

Nonultrametric{T}(n::Int, treetype::Type{T} = NamedTree) =
    Nonultrametric{T, Exponential}(n, treetype, Exponential())
Nonultrametric{T}(tiplabels::Vector{String}, treetype::Type{T} = NamedTree) =
    Nonultrametric{T, Exponential}(tiplabels, treetype, Exponential())

function rand{T, RNG}(t::Nonultrametric{T, RNG})
    t.n >= 2 || error("A tree must have at least 2 tips")
    tree = T(t.tiplabels)
    roots = NodeNameIterator(tree, isroot)
    while length(roots) > 1
        children = sample(collect(roots), 2, replace=false)
        parent = addnode!(tree)
        addbranch!(tree, parent, children[1], rand(t.rng))
        addbranch!(tree, parent, children[2], rand(t.rng))
        roots = NodeNameIterator(tree, isroot)
    end
    return tree
end

"""
    Ultrametric{T <: AbstractTree,
                RNG <: Sampleable}(n::Int,
                                   treetype::Type{T} = NamedTree,
                                   rng::RNG = Exponential())
    Ultrametric{T <: AbstractTree,
                RNG <: Sampleable}(tiplabels::Vector{String},
                                   treetype::Type{T} = NamedTree,
                                   rng::RNG = Exponential())

The sampler for ultrametric phylogenetic trees of size `n` or with
tip labels `tiplabels`. Generate random trees by calling rand().
Currently only works for `NamedTree`s.
"""
type Ultrametric{T <: AbstractTree,
                 RNG <: Sampleable{Univariate, Continuous}} <:
    Sampleable{Univariate, Phylogenetics{T}}
    n::Int
    tiplabels::Vector{String}
    rng::RNG
    
    function (::Type{Ultrametric{T, RNG}}){T, RNG}(n::Int, ::Type{T}, rng::RNG)
        return new{T, RNG}(n, map(i -> "tip $i", 1:n), rng)
    end
    
    function (::Type{Ultrametric{T, RNG}}){T, RNG}(tiplabels::Vector{String},
                                                   ::Type{T}, rng::RNG)
        return new{T, RNG}(length(tiplabels), tiplabels, rng)
    end
end

Ultrametric{T}(n::Int, treetype::Type{T} = NamedTree) =
    Ultrametric{T, Exponential}(n, treetype, Exponential())

Ultrametric{T}(tiplabels::Vector{String}, treetype::Type{T} = NamedTree) =
    Ultrametric{T, Exponential}(tiplabels, treetype, Exponential())

function rand{T, RNG}(t::Ultrametric{T, RNG})
    t.n >= 2 || error("A tree must have at least 2 tips")
    tree = T(t.tiplabels)
    roots = NodeNameIterator(tree, isroot)
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
        roots = NodeNameIterator(tree, isroot)
    end
    return tree
end
