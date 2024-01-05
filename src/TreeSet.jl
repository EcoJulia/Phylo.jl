"""
    TreeSet

A collection of trees with the same tips.
"""
mutable struct TreeSet{LABEL, RT, NL, N, B,
                       TREE <:
                       AbstractTree{OneTree, RT, NL, N, B}} <:
               AbstractTree{ManyTrees, RT, NL, N, B}
    trees::Dict{LABEL, TREE}
    treeinfo::Dict{LABEL, Dict{String, Any}}

    function TreeSet(treedict::Dict{LABEL, TREE},
                     info::Dict{LABEL, Dict{String, Any}} =
                     Dict(broadcast(key -> (key, Dict{String, Any}()),
                                    keys(treedict)))) where
             {LABEL, RT, NL, N, B, TREE <: AbstractTree{OneTree, RT, NL, N, B}}
        return new{LABEL, RT, NL, N, B, TREE}(treedict, info)
    end
end

function TreeSet(trees::AbstractVector{T}) where {T <: AbstractTree{OneTree}}
    return TreeSet(Dict(Pair.(Base.OneTo(length(trees)), trees)))
end

"""
    treesettype(::Type{AbstractTree}, ::Type{LABEL} = String)

Returns type of a TreeSet containing a collection of trees, from those trees' type
and the type of label used to identify trees.
"""
function treesettype(::Type{TREE},
                     ::Type{LABEL} = String) where
         {RT, NL, N, B, TREE <: AbstractTree{OneTree, RT, NL, N, B}, LABEL}
    return TreeSet{LABEL, RT, NL, N, B, TREE}
end

import Base.IteratorSize
function IteratorSize(::Type{TreeSet})
    return HasLength()
end

import Base.IteratorEltype
function IteratorEltype(::Type{TreeSet})
    return HasEltype()
end

import Phylo.API: _getleafinfo
_getleafinfo(ts::TreeSet) = _getleafinfo(first(gettrees(ts)))
function _getleafinfo(ts::TreeSet, leafname)
    return _getleafinfo(first(gettrees(ts)), leafname)
end

import Phylo.API: _ntrees
_ntrees(ts::TreeSet) = length(ts.trees)

import Phylo.API: _gettrees
_gettrees(ts::TreeSet) = values(ts.trees)

import Phylo.API: _gettreenames
_gettreenames(ts::TreeSet) = keys(ts.trees)

import Phylo.API: _gettreeinfo
_gettreeinfo(ts::TreeSet) = ts.treeinfo
_gettreeinfo(ts::TreeSet{LABEL}, name::LABEL) where {LABEL} = ts.treeinfo[name]

import Phylo.API: _gettree
_gettree(ts::TreeSet, name) = ts.trees[name]

import Base.eltype
function eltype(::TreeSet{LABEL, RT, NL, N, B, TREE}) where
         {LABEL, RT, NL, N, B,
          TREE <: AbstractTree{OneTree, RT, NL, N, B}}
    return TREE
end
function eltype(::Type{TreeSet{LABEL, RT, NL, N, B, TREE}}) where
         {LABEL, RT, NL, N, B,
          TREE <: AbstractTree{OneTree, RT, NL, N, B}}
    return TREE
end

import Base.length
length(ts::TreeSet) = length(ts.trees)

import Base.getindex
getindex(ts::TreeSet{LABEL}, idx::LABEL) where {LABEL} = ts.trees[idx]

import Phylo.API: _getleafnames
function _getleafnames(ts::TREESET) where {LABEL, RT, NL, N, B,
                                           TREE <:
                                           AbstractTree{OneTree, RT, NL, N, B},
                                           TREESET <:
                                           TreeSet{LABEL, RT, NL, N, B, TREE}}
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
                                      TREE <:
                                      AbstractTree{OneTree, RT, NL, N, B},
                                      TREESET <:
                                      TreeSet{LABEL, RT, NL, N, B, TREE}}
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
