using Compat.Random
import Distributions: ValueSupport, Sampleable
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

function rand(rng::AbstractRNG, t::Nonultrametric{T, SAMP}) where {T, SAMP}
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

function rand(rng::AbstractRNG, t::Ultrametric{T, SAMP}) where {T, SAMP}
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
    return TreeSet(trees)
end

rand(rng::AbstractRNG, s::S, n::Tuple{Int}) where
    {TREE <: AbstractTree{OneTree},
     S <: Sampleable{Univariate, Phylogenetics{TREE}}} =
    rand(rng, s, 1:n[1])

abstract type EvolvedTrait{T <: AbstractTree} <: ValueSupport end
Base.eltype(::Type{EvolvedTrait{T}}) where T <: AbstractTree = T
Base.eltype(::Sampleable{S, P}) where {S, P <: EvolvedTrait} = eltype(P)

struct BMTrait{T <: AbstractTree, N <: Number} <:
    Sampleable{Univariate, EvolvedTrait{T}}
    tree::T
    trait::String
    start::N
    _rand::Function
    f::Function
end

function BMTrait(tree::T, trait::String, start::N = 0.0;
                 σ² = missing, σ = missing, f::Function = identity) where
    {T <: AbstractTree, N <: Number}
    if ismissing(σ)
        σ = sqrt(σ²)
    end
    if f ≢ identity && f(start) ≉ start
            @warn "Note that the third argument (the starting state) for trait '$trait' is untransformed, so transformed version is $(f(start))"
    end
    
    return ismissing(σ) ?
          BMTrait{T, N}(tree, trait, start,
                               ((rng::AbstractRNG, start, length) ->
                                start + randn(rng) * sqrt(length)), f) :
          BMTrait{T, N}(tree, trait, start,
                               ((rng::AbstractRNG, start, length) ->
                                start + σ * randn(rng) * sqrt(length)), f)
end

function rand(rng::AbstractRNG, bm::BMTrait{TREE}) where TREE <: AbstractTree
    trait = Dict{nodetype(TREE), typeof(bm.start)}()
    use_dict = (bm.f ≢ identity)
    for node in traversal(bm.tree, preorder)
        if isroot(bm.tree, node)
            if use_dict
                trait[node] = bm.start
                setnodedata!(bm.tree, node, bm.trait, bm.f(bm.start))
            else
                setnodedata!(bm.tree, node, bm.trait, bm.start)
            end                
        else
            inb = getinbound(bm.tree, node)
            prt = src(bm.tree, inb)
            previous = use_dict ? trait[prt] :
                getnodedata(bm.tree, prt, bm.trait)
            value = bm._rand(rng, previous, getlength(bm.tree, inb))
            if use_dict
                trait[node] = value
                setnodedata!(bm.tree, node, bm.trait, bm.f(value))
            else
                setnodedata!(bm.tree, node, bm.trait, value)
            end
        end
    end
    return bm.tree
end

abstract type AbstractDiscreteTrait{T <: AbstractTree, E <: Enum} <:
    Sampleable{Univariate, EvolvedTrait{T}} end

struct DiscreteTrait{T <: AbstractTree, E <: Enum} <:
    Sampleable{Univariate, EvolvedTrait{T}}
    tree::T
    transition_matrix::Matrix{Float64}
    trait::String
end

function DiscreteTrait(tree::TREE, ttype::Type{TRAIT},
                       transition_matrix::AbstractMatrix{Float64},
                       trait::String = "$ttype") where
    {TREE <: AbstractTree, TRAIT <: Enum}
    return DiscreteTrait{TREE, TRAIT}(tree, transition_matrix, trait)
end

function enum_rand(rng::AbstractRNG, current::E,
                   mat::AbstractMatrix{Float64}) where E <: Enum
    p = @view mat[Int(current) + 1, :]
    km1 = length(p) - 1

    rp = 1.0  # remaining total probability
    i = 0
    
    while i < km1
        i += 1
        @inbounds pi = p[i]
        if pi < rp
            if rand(rng) < pi / rp
                return E(i-1)
            end
            rp -= pi
        else
            return E(i-1)
        end
    end
    
    return E(km1)
end

function rand(rng::AbstractRNG, dt::DiscreteTrait{TREE, TRAIT}) where
    {TREE <: AbstractTree, TRAIT <: Enum}
    num = Int(typemax(E))
    for node in traversal(dt.tree, preorder)
        if isroot(dt.tree, node)
            setnodedata!(dt.tree, node, dt.trait, TRAIT(sample(0:num)))
        else
            inb = getinbound(dt.tree, node)
            prt = src(dt.tree, inb)
            previous = getnodedata(dt.tree, prt, dt.trait)
            value = enum_rand(rng, previous,
                              exp(getlength(dt.tree, inb) .*
                                  dt.transition_matrix))
            setnodedata!(dt.tree, node, dt.trait, value)
        end
    end
    return dt.tree
end

struct SymmetricDiscreteTrait{T <: AbstractTree, TRAIT <: Enum} <:
    Sampleable{Univariate, EvolvedTrait{T}}
    tree::T
    transition_rate::Float64
    trait::String
end

function SymmetricDiscreteTrait(tree::TREE, ttype::Type{TRAIT},
                                transition_rate::Float64,
                                trait::String = "$ttype") where
    {TREE <: AbstractTree, TRAIT <: Enum}
    return SymmetricDiscreteTrait{TREE, TRAIT}(tree, transition_rate, trait)
end

function enum_rand(rng::AbstractRNG, current::TRAIT,
                   p_stay::Float64) where TRAIT <: Enum
    if rand(rng) < p_stay
        return current
    end

    km1 = Int(typemax(TRAIT))
    rp = km1  # remaining total options
    e = Int(current)
    passed = false
    i = Int(typemin(TRAIT))
    
    while i < km1
        passed |= (e == i)
        if rp > 1
            if rand(rng) < 1.0 / rp
                return TRAIT(i+passed)
            end
            i += 1
            rp -= 1
        else
            return TRAIT(i+passed)
        end
    end
end

function rand(rng::AbstractRNG, dt::SymmetricDiscreteTrait{TREE, TRAIT}) where
    {TREE <: AbstractTree, TRAIT <: Enum}
    num = Int(typemax(TRAIT))
    frac = 1.0 / num
    for node in traversal(dt.tree, preorder)
        if isroot(dt.tree, node)
            setnodedata!(dt.tree, node, dt.trait, TRAIT(sample(0:num)))
        else
            inb = getinbound(dt.tree, node)
            prt = src(dt.tree, inb)
            previous = getnodedata(dt.tree, prt, dt.trait)
            value = enum_rand(rng, previous,
                              frac + (1.0 - frac) *
                              exp(-num * getlength(dt.tree, inb) *
                                  dt.transition_rate))
            setnodedata!(dt.tree, node, dt.trait, value)
        end
    end
    return dt.tree
end
