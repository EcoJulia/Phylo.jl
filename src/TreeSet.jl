using Compat

mutable struct TreeSet{LABEL, RT, NL, N, B, TREE <:
                       AbstractTree{OneTree, RT, NL, N, B}} <:
    AbstractTree{ManyTrees, RT, NL, N, B}
    trees::Dict{LABEL, TREE}
    treeinfo::Dict{String, Dict{String, Any}}
end

import Compat.IteratorSize
function IteratorSize(::Type{TreeSet})
    return HasLength()
end

import Compat.IteratorEltype
function IteratorEltype(::Type{TreeSet})
    return HasEltype()
end

import Phylo.API: _ntrees
_ntrees(ts::TreeSet) = length(ts.trees)

import Phylo.API: _gettrees
_gettrees(ts::TreeSet) = values(ts.trees)

import Phylo.API: _gettreenames
_gettreenames(ts::TreeSet) = keys(ts.trees)

gettreeinfo(ts::TreeSet) = ts.treeinfo
gettreeinfo(ts::TreeSet, name) = ts.treeinfo[name]

import Phylo.API: _gettree
_gettree(ts::TreeSet, name) = ts.trees[name]

import Base.eltype
eltype(::TreeSet{LABEL, RT, NL, N, B, TREE}) where
    {LABEL, RT, NL, N, B,
     TREE <: AbstractTree{OneTree, RT, NL, N, B}} = TREE

import Base.length
length(ts::TreeSet) = length(ts.trees)

import Base.getindex
getindex(ts::TreeSet, idx) = ts.trees[idx]

import Phylo.API: _getleafnames
function _getleafnames(ts::TREESET) where {LABEL, RT, NL, N, B,
                    TREE <: AbstractTree{OneTree, RT, NL, N, B},
                    TREESET <: TreeSet{LABEL, RT, NL, N, B, TREE}}
    lns = unique(_getleafnames(t) for t in values(ts.trees))
    if length(lns) > 1
        error("Inconsistent leaf names in TreeSet")
    elseif isempty(lns)
        return NL[]
    end
    return first(lns)
end

import Phylo.API: _nleaves
function _nleaves(ts::TREESET) where {LABEL, RT, NL, N, B,
                    TREE <: AbstractTree{OneTree, RT, NL, N, B},
                    TREESET <: TreeSet{LABEL, RT, NL, N, B, TREE}}
    nls = unique(_nleaves(t) for t in values(ts.trees))
    if length(nls) > 1
        error("Inconsistent leafinfo in TreeSet")
    elseif isempty(nls)
        return 0
    end
    return first(nls)
end

import Base.keys
keys(ts::TreeSet) = keys(ts.trees)
