using Compat
import Base.start, Base.next, Base.done, Base.getindex
import Compat.IteratorSize, Base.length, Compat.IteratorEltype, Base.eltype

import Phylo.API: _nodetype, _branchtype

mutable struct TreeSet{NL, BL, TREE <: AbstractTree{NL, BL}} <: AbstractTree{NL, BL}
    trees::Dict{String, TREE}
    treeinfo::Dict{String, Dict{String, Any}}
end

TreeSet(trees::Dict{String, TREE}, treeinfo) where
        {NL, BL, TREE <: AbstractTree{NL, BL}} =
    TreeSet{NL, BL, TREE}(trees, treeinfo)

function IteratorSize(::Type{TreeSet})
    return HasLength()
end

function IteratorEltype(::Type{TreeSet})
    return HasEltype()
end

_nodetype(ts::TreeSet) = _nodetype(first(treeiter(ts)))
_branchtype(ts::TreeSet) = _branchtype(first(treeiter(ts)))
_branchtype(::TreeSet{NL, BL, TREE}) where {NL, BL, TREE <: AbstractBranchTree{NL, BL}} = Branch{NL}
ntrees(ts::TreeSet) = length(ts.trees)

struct TreeIterator{NL, BL, TREE <: AbstractTree{NL, BL},
                    TREESET <: TreeSet{NL, BL, TREE}} <: AbstractTreeIterator{TREE}
    ts::TREESET
end
treeiter(ts::TREESET) where {NL, BL, TREE <: AbstractTree{NL, BL},
                    TREESET <: TreeSet{NL, BL, TREE}} =
    TreeIterator{NL, BL, TREE, TREESET}(ts)

struct TreeNameIterator{NL, BL, TREE <: AbstractTree{NL, BL},
                    TREESET <: TreeSet{NL, BL, TREE}} <: AbstractTreeIterator{TREE}
    ts::TREESET
end
treenameiter(ts::TREESET) where {NL, BL, TREE <: AbstractTree{NL, BL},
                    TREESET <: TreeSet{NL, BL, TREE}} =
    TreeNameIterator{NL, BL, TREE, TREESET}(ts)

struct TreeInfoIterator{NL, BL, TREE <: AbstractTree{NL, BL},
                    TREESET <: TreeSet{NL, BL, TREE}} <: AbstractTreeIterator{TREE}
    ts::TREESET
end
treeinfoiter(ts::TREESET) where {NL, BL, TREE <: AbstractTree{NL, BL},
                    TREESET <: TreeSet{NL, BL, TREE}} =
    TreeInfoIterator{NL, BL, TREE, TREESET}(ts)

eltype(::TreeSet{NL, BL, TREE}) where {NL, BL, TREE <: AbstractTree{NL, BL}} = TREE
eltype(::TreeIterator{NL, BL, TREE, TREESET}) where {NL, BL, TREE <: AbstractTree{NL, BL}, TREESET <: TreeSet{NL, BL, TREE}} = TREE
eltype(::TreeNameIterator) = String
eltype(::TreeInfoIterator) = Dict{String, Any}

length(ts::TreeSet) = length(ts.trees)
length(ti::TreeIterator) = length(ti.ts.trees)
length(tni::TreeNameIterator) = length(tni.ts.trees)
length(tii::TreeInfoIterator) = length(tii.ts.treeinfo)

start(ts::TreeSet) = start(treeiter(ts))
start(ti::TreeIterator) = start(ti.ts.trees)
start(tni::TreeNameIterator) = start(keys(tni.ts.trees))
start(tii::TreeInfoIterator) = start(tii.ts.treeinfo)

next(ts::TreeSet, state) = next(treeiter(ts), state)
function next(ti::TreeIterator, state)
    v, s = next(ti.ts.trees, state)
    return v[2], s
end
next(tni::TreeNameIterator, state) = next(keys(tni.ts.trees), state)

function next(tii::TreeInfoIterator, state)
    v, s = next(tii.ts.treeinfo, state)
    return v[2], s
end

done(ts::TreeSet, state) = done(treeiter(ts), state)
done(ti::TreeIterator, state) = done(ti.ts.trees, state)
done(ti::TreeNameIterator, state) = done(keys(ti.ts.trees), state)
done(ti::TreeInfoIterator, state) = done(ti.ts.treeinfo, state)

getindex(ts::TreeSet, idx) = ts.trees[idx]
