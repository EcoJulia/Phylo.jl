using Compat
import Base.getindex
import Compat.IteratorSize, Base.length, Compat.IteratorEltype, Base.eltype

import Phylo.API: _nodetype, _branchtype
import Phylo.API: _getnodenames, _getbranchnames, _getleafnames

mutable struct TreeSet{LABEL, NL, BL, TREE <: AbstractTree{NL, BL}} <: AbstractTree{NL, BL}
    trees::Dict{LABEL, TREE}
    treeinfo::Dict{String, Dict{String, Any}}
end

TreeSet(trees::Dict{LABEL, TREE}, treeinfo) where
        {LABEL, NL, BL, TREE <: AbstractTree{NL, BL}} =
    TreeSet{LABEL, NL, BL, TREE}(trees, treeinfo)

function IteratorSize(::Type{TreeSet})
    return HasLength()
end

function IteratorEltype(::Type{TreeSet})
    return HasEltype()
end

_nodetype(ts::TreeSet) = _nodetype(first(treeiter(ts)))
_branchtype(ts::TreeSet) = _branchtype(first(treeiter(ts)))
_branchtype(::TreeSet{LABEL, NL, BL, TREE}) where {LABEL, NL, BL, TREE <: AbstractBranchTree{NL, BL}} = Branch{NL}
ntrees(ts::TreeSet) = length(ts.trees)

struct TreeIterator{LABEL, NL, BL, TREE <: AbstractTree{NL, BL},
                    TREESET <: TreeSet{LABEL, NL, BL, TREE}} <: AbstractTreeIterator{TREE}
    ts::TREESET
end
treeiter(ts::TREESET) where {LABEL, NL, BL, TREE <: AbstractTree{NL, BL},
                    TREESET <: TreeSet{LABEL, NL, BL, TREE}} =
    TreeIterator{LABEL, NL, BL, TREE, TREESET}(ts)

struct TreeNameIterator{LABEL, NL, BL, TREE <: AbstractTree{NL, BL},
                    TREESET <: TreeSet{LABEL, NL, BL, TREE}} <: AbstractTreeIterator{TREE}
    ts::TREESET
end
treenameiter(ts::TREESET) where {LABEL, NL, BL, TREE <: AbstractTree{NL, BL},
                    TREESET <: TreeSet{LABEL, NL, BL, TREE}} =
    TreeNameIterator{LABEL, NL, BL, TREE, TREESET}(ts)

struct TreeInfoIterator{LABEL, NL, BL, TREE <: AbstractTree{NL, BL},
                    TREESET <: TreeSet{LABEL, NL, BL, TREE}} <: AbstractTreeIterator{TREE}
    ts::TREESET
end
treeinfoiter(ts::TREESET) where {LABEL, NL, BL, TREE <: AbstractTree{NL, BL},
                    TREESET <: TreeSet{LABEL, NL, BL, TREE}} =
    TreeInfoIterator{LABEL, NL, BL, TREE, TREESET}(ts)

eltype(::TreeSet{LABEL, NL, BL, TREE}) where {LABEL, NL, BL, TREE <: AbstractTree{NL, BL}} = TREE
eltype(::TreeIterator{LABEL, NL, BL, TREE, TREESET}) where {LABEL, NL, BL, TREE <: AbstractTree{NL, BL}, TREESET <: TreeSet{LABEL, NL, BL, TREE}} = TREE
eltype(tni::TreeNameIterator) = eltype(keys(tni.ts.trees))
eltype(tii::TreeInfoIterator) = eltype(values(tii.ts.treeinfo))

length(ts::TreeSet) = length(ts.trees)
length(ti::TreeIterator) = length(ti.ts.trees)
length(tni::TreeNameIterator) = length(tni.ts.trees)
length(tii::TreeInfoIterator) = length(tii.ts.treeinfo)

if VERSION >= v"0.7.0-"
import Base.iterate
iterate(ts::TreeSet) = iterate(treeiter(ts))
iterate(ts::TreeSet, state) = iterate(treeiter(ts), state)
iterate(ti::TreeIterator) = iterate(values(ti.ts.trees))
iterate(ti::TreeIterator, state) = iterate(values(ti.ts.trees), state)
iterate(tni::TreeNameIterator) = iterate(keys(tni.ts.trees))
iterate(tni::TreeNameIterator, state) = iterate(keys(tni.ts.trees), state)
iterate(tii::TreeInfoIterator) = iterate(values(tii.ts.treeinfo))
iterate(tii::TreeInfoIterator, state) = iterate(values(tii.ts.treeinfo), state)
else
import Base.start, Base.next, Base.done
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
end
getindex(ts::TreeSet, idx) = ts.trees[idx]

function _getleafnames(ts::TREESET) where {LABEL, NL, BL, TREE <: AbstractTree{NL, BL}, TREESET <: TreeSet{LABEL, NL, BL, TREE}}
    lns = unique(map(t -> _getleafnames(t), values(ts.trees)))
    if length(lns) > 1
        error("Inconsistent leaf names in TreeSet")
    elseif isempty(lns)
        return NL[]
    end
    return first(lns)
end

function _getnodenames(ts::TREESET) where {LABEL, NL, BL, TREE <: AbstractTree{NL, BL}, TREESET <: TreeSet{LABEL, NL, BL, TREE}}
    nns = unique(map(t -> _getnodenames(t), values(ts.trees)))
    if length(nns) > 1
        error("Inconsistent node names in TreeSet")
    elseif isempty(nns)
        return NL[]
    end
    return first(nns)
end

function _getbranchnames(ts::TREESET) where {LABEL, NL, BL, TREE <: AbstractTree{NL, BL}, TREESET <: TreeSet{LABEL, NL, BL, TREE}}
    bns = unique(map(t -> _getbranchnames(t), values(ts.trees)))
    if length(bns) > 1
        error("Inconsistent branch names in TreeSet")
    elseif isempty(bns)
        return BL[]
    end
    return first(bns)
end

function _getleafinfo(ts::TREESET) where {LABEL, NL, BL, TREE <: AbstractTree{NL, BL}, TREESET <: TreeSet{LABEL, NL, BL, TREE}}
    lis = unique(map(t -> _getleafinfo(t), values(ts.trees)))
    if length(lis) > 1
        error("Inconsistent leafinfo in TreeSet")
    elseif isempty(lis)
        return eltype(lis)()
    end
    return first(lis)
end

function _nleaves(ts::TREESET) where {LABEL, NL, BL, TREE <: AbstractTree{NL, BL}, TREESET <: TreeSet{LABEL, NL, BL, TREE}}
    nls = unique(map(t -> _nleaves(t), values(ts.trees)))
    if length(nls) > 1
        error("Inconsistent leafinfo in TreeSet")
    elseif isempty(nls)
        return 0
    end
    return first(nls)
end
