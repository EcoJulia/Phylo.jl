using Compat.Random
import Distributions: ValueSupport, Sampleable
import Base: eltype, rand
using Phylo
using Distributions
using Missings
using IterableTables: getiterator
using Unitful
import Unitful: numtype
using Statistics

@inline numtype(::Type{R}) where R <: Real = R
iscontinuous(::Type{N}) where N <: Number = numtype(N) <: AbstractFloat

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
                      SAMP <: Sampleable{Univariate, Continuous},
                      LenUnits <: Number} <:
                          Sampleable{Univariate, Phylogenetics{T}}
    n::Int
    tiplabels::Vector{String}
    sampleable::SAMP
    leafinfo::Any
    height::LenUnits

    function Nonultrametric{T, SAMP,
                            LenUnits}(n::Int, tiplabels::Vector{String},
                                      sampleable::SAMP, leafinfo,
                                      height::LenUnits) where {T, SAMP,
                                                               LenUnits}
        return new{T, SAMP, LenUnits}(n, tiplabels, sampler(sampleable),
                                      leafinfo, height)
    end
    
    function Nonultrametric{T, SAMP,
                            LenUnits}(n::Int, sampleable::SAMP,
                                      height::LenUnits) where {T, SAMP,
                                                               LenUnits}
        return new{T, SAMP, LenUnits}(n, ["tip $i" for i in Base.OneTo(n)],
                                      sampler(sampleable), missing, height)
    end

    function Nonultrametric{T, SAMP,
                            LenUnits}(tiplabels::Vector{String},
                                      sampleable::SAMP,
                                      height::LenUnits) where {T, SAMP,
                                                               LenUnits}
        return new{T, SAMP, LenUnits}(length(tiplabels), tiplabels,
                                      sampler(sampleable), missing, height)
    end
end

function Nonultrametric{T}(n::Int, height = 1.0) where T <: AbstractTree
    return Nonultrametric{T, Exponential, typeof(height)}(n, Exponential(),
                                                          height)
end

function Nonultrametric{T}(tiplabels::Vector{String},
                           height = 1.0) where T <: AbstractTree
    return Nonultrametric{T, Exponential, typeof(height)}(tiplabels,
                                                          Exponential(),
                                                          height)
end

function Nonultrametric{T}(leafinfo,
                           height = 1.0) where T <: AbstractTree
    tipnames = unique([info[1] for info in getiterator(leafinfo)])
    return Nonultrametric{T, Exponential,
                          typeof(height)}(length(tipnames), tipnames,
                                          Exponential(), leafinfo, height)
end

Nonultrametric(info::LI, height::U = 1.0) where {LI, U} =
    Nonultrametric{Phylo.LT{OneRoot, LI, U}}(info)

Nonultrametric(n::Int, height::U = 1.0) where U =
    Nonultrametric{Phylo.LTD{OneRoot, U}}(n, height)
Nonultrametric(tiplabels::Vector{String}, height::U = 1.0) where U =
    Nonultrametric{Phylo.LTD{OneRoot, U}}(tiplabels, height)

function rand(rng::AbstractRNG, t::Nonultrametric{T, SAMP, U}) where
    {T, SAMP, U}
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
        createbranch!(tree, parent, children[1],
                      rand(rng, t.sampleable) * t.height)
        createbranch!(tree, parent, children[2],
                      rand(rng, t.sampleable) * t.height)
    end
    scale = t.height / mean(getheight.(tree, getleaves(tree)))
    for branch in collect(getbranches(tree))
        source = src(tree, branch)
        destination = dst(tree, branch)
        len = uconvert(unit(t.height), getlength(tree, branch) * scale)
        deletebranch!(tree, branch)
        createbranch!(tree, source, destination, len)
    end
    return tree
end

"""
    Ultrametric{T <: AbstractTree,
                SAMP <: Sampleable,
                LenUnits <: Number}(n::Int,
                                    sampleable::SAMP = Exponential())
    Ultrametric{T <: AbstractTree,
                SAMP <: Sampleable,
                LenUnits <: Number}(tiplabels::Vector{String},
                                    sampleable::SAMP = Exponential())

The sampler for ultrametric phylogenetic trees of size `n` or with
tip labels `tiplabels`. Generate random trees by calling rand().
"""
struct Ultrametric{T <: AbstractTree,
                   SAMP <: Sampleable{Univariate, Continuous},
                   LenUnits <: Number} <:
                       Sampleable{Univariate, Phylogenetics{T}}
    n::Int
    tiplabels::Vector{String}
    sampleable::SAMP
    leafinfo::Any
    height::LenUnits

    function Ultrametric{T, SAMP, U}(n::Int, tiplabels::Vector{String},
                                     sampleable::SAMP, leafinfo,
                                     height::U = 1.0) where {T, SAMP, U}
        return new{T, SAMP, U}(n, tiplabels, sampler(sampleable),
                               leafinfo, height)
    end

    function Ultrametric{T, SAMP, U}(n::Int, sampleable::SAMP,
                                     height::U = 1.0) where {T, SAMP, U}
        return new{T, SAMP, U}(n, ["tip $i" for i in Base.OneTo(n)],
                               sampler(sampleable), missing, height)
    end

    function Ultrametric{T, SAMP, U}(tiplabels::Vector{String},
                                     sampleable::SAMP,
                                     height::U = 1.0) where
        {T, SAMP, U}
        return new{T, SAMP, U}(length(tiplabels), tiplabels,
                               sampler(sampleable), missing, height)
    end
end

function Ultrametric{T}(n::Int, height = 1.0) where T <: AbstractTree
    return Ultrametric{T, Exponential, typeof(height)}(n, Exponential(), height)
end

function Ultrametric{T}(tiplabels::Vector{String},
                        height = 1.0) where T <: AbstractTree
    return Ultrametric{T, Exponential, typeof(height)}(tiplabels, Exponential(),
                                                       height)
end

function Ultrametric{T}(leafinfo, height = 1.0) where T <: AbstractTree
    tipnames = unique([info[1] for info in getiterator(leafinfo)])
    return Ultrametric{T, Exponential, typeof(height)}(length(tipnames),
                                                       tipnames,
                                                       Exponential(), leafinfo,
                                                       height)
end

Ultrametric(info::LI, height::U = 1.0) where {LI, U} =
    Ultrametric{Phylo.LT{OneRoot, LI, U}}(info)
Ultrametric(n::Int, height::U = 1.0) where U =
    Ultrametric{Phylo.LTD{OneRoot, U}}(n, height)
Ultrametric(tiplabels::Vector{String}, height::U = 1.0) where U =
    Ultrametric{Phylo.LTD{OneRoot, U}}(tiplabels, height)

function rand(rng::AbstractRNG, t::Ultrametric{T, SAMP}) where {T, SAMP}
    t.n >= 2 || error("A tree must have at least 2 tips")
    if ismissing(t.leafinfo)
        tree = T(t.tiplabels)
    else
        tree = T(t.leafinfo)
    end
    depth = zero(rand(rng, t.sampleable)) * t.height
    leaves = getleaves(tree)
    while nroots(tree) > 1
        roots = getroots(tree)
        tocoalesce = collect(roots)
        coalescers = sample(tocoalesce, 2, replace=false)
        parent = createnode!(tree)
        depth += rand(rng, t.sampleable) * 2.0 / length(tocoalesce) * t.height
        d1 = getheight(tree, first(nodefuture(tree, coalescers[1]) ∩
                                   leaves))
        d2 = getheight(tree, first(nodefuture(tree, coalescers[2]) ∩
                                   leaves))
        createbranch!(tree, parent, coalescers[1], depth - d1)
        createbranch!(tree, parent, coalescers[2], depth - d2)
    end
    scale = t.height / getheight(tree, first(getleaves(tree)))
    for branch in collect(getbranches(tree))
        source = src(tree, branch)
        destination = dst(tree, branch)
        len = uconvert(unit(t.height), getlength(tree, branch) * scale)
        deletebranch!(tree, branch)
        createbranch!(tree, source, destination, len)
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

"""
    BrownianTrait{T <: AbstractTree, N <: Number}

A continuous trait evolved on a phylogenetic tree. This is a
`Sampleable` type, so a random trait can be created using
`rand()`. The trait to be evolved can be any continuous numeric type,
including `Unitful` types for instance, and in the simplest case is
determined by the third argument to the constructor `start`:

function BrownianTrait(tree::AbstractTree, trait::String, start::Number = 0.0;
                       σ² = missing, σ = missing, f::Function = identity)

Note that when `Unitful` is being used, either here or in branch
lengths, `σ`/`σ²` keyword argument units must be appropriate. The final keyword argument, `f`, is a function to transform the evolved gaussian trait into its true value. By default this is the identity function, but can, for instance, be `abs` to force a positive value on the trait, or more complex functions as required, such as a transformation to turn a continuous variable into a discrete trait
"""
struct BrownianTrait{T <: AbstractTree, N <: Number} <:
    Sampleable{Univariate, EvolvedTrait{T}}
    tree::T
    trait::String
    start::N
    _rand::Function
    f::Function
end

function BrownianTrait(tree::T, trait::String, start::N = 0.0;
                       σ² = missing, σ = missing, f::Function = identity) where
    {T <: AbstractTree, N <: Number}
    iscontinuous(N) ||
        throw(TypeError("Type $N must be continuous for a Gaussian Trait"))
    if ismissing(σ)
        σ = sqrt(σ²)
    end
    
    f ≢ identity && f(start) ≉ start &&
        @warn "Note that the third argument (the starting state) for trait '$trait' is untransformed - transformed version is $(f(start))"
    
    if ismissing(σ)
        dimension(start) ==
            dimension(sqrt(getlength(tree, first(getbranches(tree))))) ||
            throw(DimensionMismatch("Dimensions of start, σ[²] and branch lengths must combine correctly if using Unitful"))
        return BrownianTrait{T, N}(tree, trait, start,
                                   ((rng::AbstractRNG, start, length) ->
                                    start + randn(rng) * sqrt(length)), f)
    else
        dimension(start) == dimension(σ*sqrt(getlength(tree,
                                                       first(getbranches(tree))))) ||
                                                      throw(DimensionMismatch("Dimensions of start, σ[²] and branch lengths must combine correctly if using Unitful"))
        return BrownianTrait{T, N}(tree, trait, start,
                                   ((rng::AbstractRNG, start, length) ->
                                    start + σ * randn(rng) * sqrt(length)), f)
    end
end

function rand(rng::AbstractRNG,
              bm::BrownianTrait{TREE}) where TREE <: AbstractTree
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

"""
    DiscreteTrait{T <: AbstractTree, E <: Enum}

A discrete trait evolved on a phylogenetic tree. This is a
`Sampleable` type, so a random trait can be created using
`rand(dt)`. The trait to be evolved must be an Enum (generally created
using `@enum`), and is the second argument to the constructor:

function DiscreteTrait(tree::AbstractTree, ttype::Type{<:Enum},
                       transition_matrix::AbstractMatrix{Float64},
                       trait::String = "\$ttype")

The transition matrix holds transition rates from row to column (so
row sums must be zero), and the transition probabilities in a branch
are calculated as `exp(transition_matrix .* branch_length)`.
"""
struct DiscreteTrait{T <: AbstractTree, E <: Enum, N <: Number} <:
    Sampleable{Univariate, EvolvedTrait{T}}
    tree::T
    transition_matrix::Matrix{N}
    trait::String
end

function DiscreteTrait(tree::TREE, ttype::Type{TRAIT},
                       transition_matrix::AbstractMatrix{N},
                       trait::String = "$ttype") where
    {TREE <: AbstractTree, TRAIT <: Enum, N <: Number}
    all(size(transition_matrix) .== length(instances(TRAIT))) ||
        error("Transition matrix size must match number of traits in $TRAIT")
    dimension(first(transition_matrix)) ==
        NoDims / dimension(getlength(tree, first(getbranches(tree)))) ||
        throw(DimensionMismatch("Dimensions of transition_matrix and branch lengths must cancel if using Unitful units"))
    dimension(getlength(tree, first(getbranches(tree))) *
              first(transition_matrix)) == NoDims ||
              throw(DimensionMismatch("Dimensions of transition_matrix and branch lengths must cancel correctly if using Unitful"))
    all(sum(transition_matrix .* getlength(tree, first(getbranches(tree))),
            dims=2) .+ 1 .≈ 1) ||
                error("Row sums of transition matrix must be zero")
    
    return DiscreteTrait{TREE, TRAIT, N}(tree, transition_matrix, trait)
end

function enum_rand(rng::AbstractRNG, current::E,
                   mat::AbstractMatrix{Float64}) where E <: Enum
    traits = instances(E)
    row = findfirst(x -> x == current, traits)
    p = @view mat[row, :]
    tot = length(traits)

    rp = 1.0  # remaining total probability
    i = 1
    while i < tot
        @inbounds pi = p[i]
        if pi < rp
            if rand(rng) < pi / rp
                return traits[i]
            end
            rp -= pi
        else
            return traits[i]
        end
        i += 1
    end
    
    return traits[tot]
end

function rand(rng::AbstractRNG, dt::DiscreteTrait{TREE, TRAIT}) where
    {TREE <: AbstractTree, TRAIT <: Enum}
    num = Int(typemax(TRAIT))
    for node in traversal(dt.tree, preorder)
        if isroot(dt.tree, node)
            setnodedata!(dt.tree, node, dt.trait, TRAIT(sample(0:num)))
        else
            inb = getinbound(dt.tree, node)
            prt = src(dt.tree, inb)
            previous = getnodedata(dt.tree, prt, dt.trait)
            value = enum_rand(rng, previous,
                              exp(uconvert.(NoUnits,
                                            getlength(dt.tree, inb) .*
                                            dt.transition_matrix)))
            setnodedata!(dt.tree, node, dt.trait, value)
        end
    end
    return dt.tree
end

"""
    SymmetricDiscreteTrait{T <: AbstractTree, E <: Enum}

The simplest possible discrete trait evolved on a phylogenetic
tree. This is a `Sampleable` type, so a random trait can be created
using `rand(sdt)`. The trait to be evolved must be an Enum (generally
created using `@enum`), and is the second argument to the constructor:

function DiscreteTrait(tree::AbstractTree, ttype::Type{<:Enum},
                       transition_rate::Number,
                       trait::String = "\$ttype")

The transition matrix holds transition rates from row to column (so
row sums must be zero), and the transition probabilities in a branch
are calculated as `exp(transition_matrix .* branch_length)`.
"""
struct SymmetricDiscreteTrait{T <: AbstractTree, E <: Enum, N <: Number} <:
    Sampleable{Univariate, EvolvedTrait{T}}
    tree::T
    transition_rate::N
    trait::String
end

function SymmetricDiscreteTrait(tree::TREE, ttype::Type{TRAIT},
                                transition_rate::N,
                                trait::String = "$ttype") where
    {TREE <: AbstractTree, TRAIT <: Enum, N <: Number}
    dimension(getlength(tree, first(getbranches(tree))) *
              transition_rate) == NoDims ||
              throw(DimensionMismatch("Dimensions of transition_matrix and branch lengths must cancel correctly if using Unitful"))

    return SymmetricDiscreteTrait{TREE, TRAIT, N}(tree, transition_rate, trait)
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
