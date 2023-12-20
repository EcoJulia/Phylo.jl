using Random
import Distributions: ValueSupport, Sampleable
import Base: eltype, rand
import Random: rand!
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

Random.rand!(s::Sampleable, t::AbstractTree) = rand!(Random.GLOBAL_RNG, s, t)

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
    Nonultrametric{Phylo.ReT{OneRoot, LI, PolytomousBranching, U}}(info)

Nonultrametric(n::Int, height::U = 1.0) where U =
    Nonultrametric{Phylo.ReTD{OneRoot, PolytomousBranching, U}}(n, height)
Nonultrametric(tiplabels::Vector{String}, height::U = 1.0) where U =
    Nonultrametric{Phylo.ReTD{OneRoot, PolytomousBranching, U}}(tiplabels, height)

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
        children = sample(rng, collect(roots), 2, replace=false)
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
    Ultrametric{Phylo.ReT{OneRoot, LI, PolytomousBranching, U}}(info)
Ultrametric(n::Int, height::U = 1.0) where U =
    Ultrametric{Phylo.ReTD{OneRoot, PolytomousBranching, U}}(n, height)
Ultrametric(tiplabels::Vector{String}, height::U = 1.0) where U =
    Ultrametric{Phylo.ReTD{OneRoot, PolytomousBranching, U}}(tiplabels, height)

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
        coalescers = sample(rng, tocoalesce, 2, replace=false)
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
        throw(TypeError(:BrownianTrait, "Type $N must be continuous for a Gaussian Trait", AbstractFloat, N))
    if ismissing(σ)
        σ = sqrt(σ²)
    end
    
    f ≢ identity && f(start) ≉ start &&
        @warn "Note that the third argument (the starting state) for trait '$trait' is untransformed - transformed version is $(f(start))"
    
    if ismissing(σ)
        dimension(N) ≡
            dimension(sqrt(getlength(tree, first(getbranches(tree))))) ||
            throw(DimensionMismatch("Dimensions of start, σ[²] and branch lengths must combine correctly if using Unitful"))
        return BrownianTrait{T, N}(tree, trait, start,
                                   ((rng::AbstractRNG, start::N, length) ->
                                    start + randn(rng, typeof(one(start))) *
                                    sqrt(length)), f)
    else
        dimension(N) ≡
            dimension(σ*sqrt(getlength(tree,
                                       first(getbranches(tree))))) ||
                                           throw(DimensionMismatch("Dimensions of start, σ[²] and branch lengths must combine correctly if using Unitful"))
        return BrownianTrait{T, N}(tree, trait, start,
                                   ((rng::AbstractRNG, start::N, length) ->
                                    start +
                                    σ * randn(rng, typeof(one(start))) *
                                    sqrt(length)), f)
    end
end

function rand!(rng::AbstractRNG,
               bm::BrownianTrait{TREE, N},
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
              bm::BrownianTrait{TREE, N}) where {TREE <: AbstractTree,
                                                 N <: Number}
    untrait = Dict{nodetype(TREE), N}()
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
                       trait::String)

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

function enum_rand(rng::AbstractRNG, current::TRAIT,
                   mat::AbstractMatrix{Float64}) where TRAIT <: Enum
    traits = instances(TRAIT)
    row = findfirst(x -> x == current, traits)
    p = @view mat[row, :]
    tot = length(traits)

    rp = 1.0  # remaining total probability
    i = 1
    while i < tot
        @inbounds pi = p[i]
        if pi < rp
            if rand(rng, Float64) < pi / rp
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

function rand!(rng::AbstractRNG, dt::DiscreteTrait{TREE, TRAIT},
               tree::TREE) where {TREE <: AbstractTree, TRAIT <: Enum}
    traits = instances(TRAIT)
    num = length(traits)
    for node in traversal(tree, preorder)
        if isroot(tree, node)
            setnodedata!(tree, node, dt.trait,
                         sample(rng, collect(traits)))
        else
            inb = getinbound(tree, node)
            prt = src(tree, inb)
            previous = getnodedata(tree, prt, dt.trait)
            value = enum_rand(rng, previous,
                              exp(uconvert.(NoUnits,
                                            getlength(tree, inb) .*
                                            dt.transition_matrix)))
            setnodedata!(tree, node, dt.trait, value)
        end
    end
    return tree
end

"""
    SymmetricDiscreteTrait{T <: AbstractTree, E <: Enum}

The simplest possible discrete trait evolved on a phylogenetic
tree. This is a `Sampleable` type, so a random trait can be created
using `rand(sdt)`. The trait to be evolved must be an Enum (generally
created using `@enum`), and is the second argument to the constructor:

function DiscreteTrait(tree::AbstractTree, ttype::Type{<:Enum},
                       transition_rate::Number,
                       trait::String)

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
    if rand(rng, Float64) < p_stay
        return current
    end
    
    traits = instances(TRAIT)
    num = length(traits)
    remain = num - 1  # remaining total options
    from = findfirst(x -> x == current, traits)
    passed = false
    i = 1
    while i < num
        passed |= (from == i)
        if remain > 1
            if rand(rng, Float64) < inv(remain)
                return traits[i+passed]
            end
            i += 1
            remain -= 1
        else
            return traits[i+passed]
        end
    end
end

function rand!(rng::AbstractRNG, dt::SymmetricDiscreteTrait{TREE, TRAIT},
               tree::TREE) where {TREE <: AbstractTree, TRAIT <: Enum}
    traits = instances(TRAIT)
    num = length(traits)
    frac = inv(num)
    for node in traversal(tree, preorder)
        if isroot(tree, node)
            setnodedata!(tree, node, dt.trait,
                         sample(rng, collect(instances(TRAIT))))
        else
            inb = getinbound(tree, node)
            prt = src(tree, inb)
            previous = getnodedata(tree, prt, dt.trait)
            value = enum_rand(rng, previous,
                              frac + (1.0 - frac) *
                              exp(-num * getlength(tree, inb) *
                                  dt.transition_rate))
            setnodedata!(tree, node, dt.trait, value)
        end
    end
    return tree
end
