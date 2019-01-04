using LightGraphs, MetaGraphs, Missings
using LightGraphs: SimpleEdge

struct SimpleNode{RT, NL} <: AbstractNode{RT, NL}
    id::Int
end
==(e1::SimpleNode, e2::SimpleNode) = (e1.id == e2.id)

struct SimpleBranch{RT, NL} <: AbstractBranch{RT, NL}
    src::Int
    dst::Int
end
SimpleBranch{RT, NL}(edge::SimpleEdge{Int}) where {RT, NL} =
    SimpleBranch{RT, NL}(src(edge), dst(edge))
==(e1::SimpleBranch, e2::SimpleBranch) = (e1.src == e2.src && e1.dst == e2.dst)

abstract type MetaTree{RT, NL, TD} <:
    AbstractTree{OneTree, RT, NL, SimpleNode{RT, NL}, SimpleBranch{RT, NL}} end

_noderecordtype(::Type{<:MetaTree}) = Dict{Symbol, Any}
_branchrecordtype(::Type{<:MetaTree}) = Dict{Symbol, Any}

mutable struct UnrootedMetaTree{NL, TD} <: MetaTree{Unrooted, NL, TD}
    mg::MetaGraph{Int, Float64}
    isvalid::Union{Bool, Missing}
    name::Union{String, Missing}
    tipdata::TD
end

mutable struct RootedMetaTree{RT <: Rooted, NL, TD} <: MetaTree{RT, NL, TD}
    mg::MetaDiGraph{Int, Float64}
    roots::Vector{SimpleNode{RT, NL}}
    isvalid::Union{Bool, Missing}
    name::Union{String, Missing}
    tipdata::TD
end

const RootedTree = RootedMetaTree{OneRoot, String, Dict{String, Any}}
const ManyRootTree = RootedMetaTree{ManyRoots, String, Dict{String, Any}}
const UnrootedTree = UnrootedMetaTree{String, Dict{String, Any}}
function RootedTree(n::Int = 0)
    tree = RootedTree(MetaDiGraph(), Vector{SimpleNode{OneRoot, String}}(),
                      missing, missing, Dict{String, Any}())
    set_indexing_prop!(tree.mg, :name)
    createnodes!(tree, n)

    return tree
end
function RootedTree(names::Vector{String})
    tree = RootedTree(MetaDiGraph(), Vector{SimpleNode{OneRoot, String}}(),
                      missing, missing, Dict{String, Any}())
    set_indexing_prop!(tree.mg, :name)
    createnodes!(tree, names)

    return tree
end

import Phylo.API._nnodes
_nnodes(tree::MetaTree) = nv(tree.mg)
import Phylo.API._getnodes
_getnodes(tree::T) where {RT, NL, TD, T <: MetaTree{RT, NL, TD}} =
    SimpleNode{RT, NL}.(vertices(tree.mg))
import Phylo.API._getnodenames
_getnodenames(tree::MetaTree) =
    String[_getnodename(tree, v) for v in _getnodes(tree)]
import Phylo.API._getnode
_getnode(tree::MetaTree{RT, NL, TD}, nodename::NL) where {RT, NL, TD} =
    SimpleNode{RT, NL}(tree[nodename, :name])
_getnode(tree::MetaTree{RT, NL, TD}, nodename::NL) where {RT, NL<:Integer, TD} =
    SimpleNode{RT, NL}(tree[:name][nodename])
import Phylo.API._getnodename
_getnodename(tree::MetaTree, node::SimpleNode) = tree.mg[node.id, :name]
import Phylo.API._hasnode
_hasnode(tree::MetaTree, node::SimpleNode) = 1 ≤ node.id ≤ _nnodes(tree)
_hasnode(tree::MetaTree, id::Int) = 1 ≤ id ≤ _nnodes(tree)

import Phylo.API._nbranches
_nbranches(tree::MetaTree) = ne(tree.mg)
import Phylo.API._getbranches
_getbranches(tree::T) where {RT, NL, TD, T <: MetaTree{RT, NL, TD}} =
    SimpleBranch{RT, NL}.(edges(tree.mg))
import Phylo.API._getbranchnames
_getbranchnames(tree::MetaTree) =
    [node + (idx - 1) * _nnodes(tree)
     for node in vertices(tree)
     for idx in Base.OneTo(length(fadj(tree, node)))]
import Phylo.API._getbranch
function _getbranch(tree::MetaTree{RT, NL, TD}, id::Int) where {RT, NL, TD}
    node = id % _nnodes(tree)
    idx = id ÷ _nnodes(tree)
    return SimpleBranch{RT, NL}(node, fadj(tree, node)[idx])
end
import Phylo.API._getbranchname
function _getbranchname(tree::MetaTree, branch::SimpleBranch)
    idx = findfirst(dest -> dest == branch.dst, fadj(tree, branch.src))
    idx > 0 || error("Branch not found")
    return branch.src + (idx - 1) * _nnodes(tree)
end
import Phylo.API._hasbranch
function _hasbranch(tree::MetaTree, branch::SimpleBranch)
    idx = findfirst(dest -> dest == branch.dst, fadj(tree, branch.src))
    return idx > 0
end
function _hasbranch(tree::MetaTree{RT, NL, TD}, id::Int) where {RT, NL, TD}
    node = id % _nnodes(tree)
    idx = id ÷ _nnodes(tree)
    return SimpleBranch{RT, NL}(node, fadj(tree, node)[idx])
end

import Phylo.API._src
_src(tree::RootedMetaTree, branch::SimpleBranch{<:Rooted}) =
     _getnode(tree, branch.src)
_src(tree::RootedMetaTree, id::Int) =
     _getnodename(tree, _getbranch(tree, id).src)

 import Phylo.API._dst
_dst(tree::RootedMetaTree, branch::SimpleBranch{<:Rooted}) =
     _getnode(tree, branch.dst)
_dst(tree::RootedMetaTree, id::Int) =
     _getnodename(tree, _getbranch(tree, id).dst)

 import Phylo.API._getroots
_getroots(tree::RootedMetaTree) = tree.roots

import Phylo.API._indegree
_indegree(tree::RootedMetaTree, node) =
    length(badj(tree.mg, _getnode(tree, node).id))
import Phylo.API._outdegree
_outdegree(tree::RootedMetaTree, node) =
    length(fadj(tree.mg, _getnode(tree, node).id))
import Phylo.API._degree
_degree(tree::UnrootedMetaTree, node) =
    length(fadj(tree.mg, _getnode(tree, node).id))

import Phylo.API._getparent
_getparent(tree::RootedMetaTree{RT, NL, TD},
           node::SimpleNode{RT, NL}) where {RT, NL, TD} =
    SimpleNode{RT, NL}(first(badj(tree.mg, node.id)))
_getparent(tree::RootedMetaTree{RT, NL, TD}, node::NL) where {RT, NL, TD} =
    _getnodename(tree, _getparent(tree, _getnode(tree, node)))
import Phylo.API._getinbound
_getinbound(tree::RootedMetaTree{RT, NL, TD},
            node::SimpleNode{RT, NL}) where {RT, NL, TD} =
    SimpleBranch{RT, NL}(first(badj(tree.mg, node.id)), node.id)
_getinbound(tree::RootedMetaTree{RT, NL, TD}, node::NL) where {RT, NL, TD} =
    _getbranchname(tree, _getinbound(tree, _getnode(tree, node)))

import Phylo.API._getchildren
_getchildren(tree::RootedMetaTree{RT, NL, TD},
              node::SimpleNode{RT, NL}) where {RT, NL, TD} =
    [SimpleNode{RT, NL}(out) for out in fadj(tree.mg, node.id)]
_getchildren(tree::RootedMetaTree{RT, NL, TD}, node::NL) where {RT, NL, TD} =
    [_getnodename(tree, out)
     for out in _getchildren(tree, _getnode(tree, node))]
import Phylo.API._getoutbounds
_getoutbounds(tree::RootedMetaTree{RT, NL, TD},
              node::SimpleNode{RT, NL}) where {RT, NL, TD} =
    [SimpleBranch{RT, NL}(node, out)
     for out in fadj(tree.mg, _getnode(tree, node).id)]
_getoutbounds(tree::RootedMetaTree{RT, NL, TD}, node::NL) where {RT, NL, TD} =
    [_getbranchname(tree, branch)
     for branch in _getoutbounds(tree, _getnode(tree, node))]

import Phylo.API._getsiblings
_getsiblings(tree::RootedMetaTree{RT, NL, TD},
             node::SimpleNode{RT, NL}) where {RT, NL, TD} =
    [SimpleNode{RT, NL}(out) for out in fadj(tree.mg, node.id)]
_getsiblings(tree::RootedMetaTree{RT, NL, TD}, node::NL) where {RT, NL, TD} =
    [_getnodename(tree, out)
     for out in _getchildren(tree, _getnode(tree, node))]
import Phylo.API._getconnections
_getconnections(tree::UnrootedMetaTree{NL, TD}, node) where {NL, TD} =
    SimpleNode{Unrooted, NL}.(fadj(tree.mg, _getnode(tree, node).id))
_getconnections(tree::RootedMetaTree{RT, NL, TD}, node::NL) where {RT, NL, TD} =
    [_getbranchname(tree, branch)
     for branch in _getconnections(tree, _getnode(tree, node))]

import Phylo.API._addinbound!
_addinbound!(::RootedMetaTree, node) = Nothing()
import Phylo.API._removeinbound!
_removeinbound!(::RootedMetaTree, node) = Nothing()
import Phylo.API._addoutbound!
_addoutbound!(::RootedMetaTree, node) = Nothing()
import Phylo.API._removeoutbound!
_removeoutbound!(::RootedMetaTree, node) = Nothing()
import Phylo.API._addconnection!
_addconnection!(::UnrootedMetaTree, node) = Nothing()
import Phylo.API._removeconnection!
_removeconnection!(::UnrootedMetaTree, node) = Nothing()

import Phylo.API._createbranch!
function _createbranch!(tree::MetaTree{RT, NL, TD},
                        source::SimpleNode{RT, NL},
                        dest::SimpleNode{RT, NL},
                        length::Float64, data) where {RT, NL, TD}
    if ismissing(data)
        add_edge!(tree, source.id, dest.id)
        add_prop!(tree, source.id, dest.id, :length, length)
    else
        data[:length] = length
        add_edge!(tree, source.id, dest.id, data)
    end
    branch = SimpleBranch{RT, NL}(source.id, dest.id)
    filter!(n -> n != dest.id, tree.roots)
    return branch
end
_createbranch!(tree::MetaTree{RT, NL, TD}, source::NL, dest::NL,
               length::Float64, data) where {RT, NL, TD} =
    _createbranch!(tree, _getnode(tree, source), _getnode(tree, dest),
                   length, data)

import Phylo.API._removebranch!
function _removebranch!(tree::MetaTree{RT, NL, TD},
                        source::SimpleNode{RT, NL},
                        dest::SimpleNode{RT, NL},
                        length::Float64, data) where {RT, NL, TD}
    rem_edge!(tree, source.id, dest.id, data)
    push!(tree.roots, dest)
end
_removebranch!(tree::MetaTree{RT, NL, TD}, source::NL, dest::NL) where
    {RT, NL, TD} =
    _removebranch!(tree, _getnode(tree, source), _getnode(tree, dest))

import Phylo.API._createnode!
function _createnode!(tree::T, nodename::NL, data) where
    {RT, NL, TD, T <: MetaTree{RT, NL, TD}}
    if ismissing(data)
        add_vertex!(tree.mg, noderecordtype(T)())
    else
        add_vertex!(tree.mg, data)
    end
    node = SimpleNode{RT, NL}(nv(tree.mg))
    set_indexing_prop!(tree.mg, node.id, :name, nodename)
    return node
end

import Phylo.API._deletenode!
function _deletenode!(tree::RootedMetaTree, node::SimpleNode)
    rem_vertex!(tree.mg, node.id)
    filter!(n -> n != node.id, tree.roots)
end
function _deletenode!(tree::UnrootedMetaTree, node::SimpleNode)
    rem_vertex!(tree.mg, node.id)
end
_deletenode!(tree::MetaTree{RT, NL, TD}, nodename::NL) where {RT, NL, TD} =
    _deletenode!(tree, _getnode(tree, nodename))
