using Missings

struct LinkBranch{RT, NL} <: AbstractBranch{RT, NL}
    inout::Tuple{AbstractNode{RT, NL}, AbstractNode{RT, NL}}
    length::Float64
    data::Dict{String, Any}
end

mutable struct LinkNode{RT, NL} <: AbstractNode{RT, NL}
    name::NL
    inbound::Union{AbstractBranch{RT, NL}, Missing}
    other::Vector{AbstractBranch{RT, NL}}
    data::Dict{String, Any}
end

mutable struct LinkTree{RT, NL, N <: LinkNode{RT, NL},
                        B <: LinkBranch{RT, NL}, TD} <:
               AbstractTree{OneTree, RT, NL, N, B}
    name::Union{String, Missing}
    nodes::Dict{NL, N}
    roots::Vector{N}
    branches::Vector{Union{B, Missing}}
    data::Dict{String, Any}
    tipdata::TD
    rootheight::Float64
    isvalid::Union{Bool, Missing}

    function LinkTree{RT, NL, N, B, TD}(tipnames::Vector{NL} = NL[];
        treename::Union{String, Missing} = missing,
        tipdata::TD = TD(), rootheight = NaN) where {RT, NL, N, B, TD}
        tree = new{RT, NL, N, B, TD}(treename, Dict{NL, N}(), Vector{N}(),
                                     Vector{B}(), Dict{String, Any}(),
                                     tipdata, rootheight, missing)
        if !isempty(tipnames)
            createnodes!(tree, tipnames)
        elseif !isempty(tipdata)
            tips = unique(keys(tipdata))
            createnodes!(tree, tips)
        end
        return tree
    end
end
const LB{RT} = LinkBranch{RT, String}
const LN{RT} = LinkNode{RT, String}
const LT{RT, TD} = LinkTree{RT, String, LN{RT}, LB{RT}, TD}
const RootedTree = LT{OneRoot, Dict{String, Any}}
const ManyRootTree = LT{ManyRoots, Dict{String, Any}}
const UnrootedTree = LT{Unrooted, Dict{String, Any}}

# LinkBranch methods
function LinkBranch(from::LinkNode{RT, NL},
                    to::LinkNode{RT, NL},
                    length::Float64 = NaN,
                    data::Dict{String, Any} = Dict{String, Any}()) where
    {RT, NL}
    return LinkBranch{RT, NL}((from, to), length, data)
end

import Phylo.API._src
Phylo.API._src(branch::LinkBranch{RT, NL}) where {RT, NL} =
    branch.inout[1]

import Phylo.API._dst
Phylo.API._dst(branch::LinkBranch{RT, NL}) where {RT, NL} =
    branch.inout[2]

import Phylo.API._getlength
Phylo.API._getlength(branch::LinkBranch{RT, NL}) where {RT, NL} =
    branch.length

import Phylo.API._branchrecordtype
Phylo.API._branchrecordtype(::Type{LinkBranch{RT, NL}}) where {RT, NL} =
    Dict{String, Any}

import Phylo.API._getbranchinfo
Phylo.API._getbranchinfo(branch::LinkBranch{RT, NL}) where {RT, NL} =
    branch.data

import Phylo.API._setbranchinfo!
Phylo.API._setbranchinfo!(branch::LinkBranch{RT, NL},
                          data::Dict{String, Any}) where {RT, NL} =
    branch.data = data

# LinkNode methods
function LinkNode(tree::T, name::NL,
                  data::Dict{String, Any} = Dict{String, Any}()) where
    {RT, NL, N <: LinkNode, B, T <: AbstractTree{OneTree, RT, NL, N, B}}
    return N(name, missing, Vector{B}(), data)
end

import Phylo.API._noderecordtype
Phylo.API._noderecordtype(::Type{LinkNode{RT, NL}}) where {RT, NL} =
    Dict{String, Any}

import Phylo.API._getnoderecord
Phylo.API._getnoderecord(node::LinkNode{RT, NL}) where {RT, NL} =
    node.data

import Phylo.API._setnoderecord!
Phylo.API._setnoderecord!(node::LinkNode{RT, NL},
                        data::Dict{String, Any}) where {RT, NL} =
    node.data = data

import Phylo.API._hasinbound
_hasinbound(node::LinkNode{RT, NL}) where {RT <: Rooted, NL} =
    !ismissing(node.inbound)

import Phylo.API._degree
_degree(node::LinkNode{Unrooted, NL}) where {NL} =
    length(node.other)

import Phylo.API._getinbound
_getinbound(node::LinkNode) = node.inbound

import Phylo.API._addinbound!
function _addinbound!(node::LinkNode{RT, NL}, branch::B) where
         {RT <: Rooted, NL, B <: AbstractBranch{RT, NL}}
    !_hasinbound(node) ||
        error("LinkNode already has an inbound connection")
    node.inbound = branch
end

import Phylo.API._removeinbound!
function _removeinbound!(node::LinkNode{RT, NL}, branch::B) where
    {RT <: Rooted, NL, B <: LinkBranch{RT, NL}}
    _hasinbound(node) ||
        error("Node has no inbound connection")
    node.inbound === branch ||
        error("Node has no inbound connection from branch $inbound")
    node.inbound = missing
end

import Phylo.API._getoutbounds
_getoutbounds(node::LinkNode{RT, NL}) where
    {RT <: Rooted, NL} = node.other

import Phylo.API._addoutbound!
function _addoutbound!(node::LinkNode{RT, NL}, branch::B) where
         {RT <: Rooted, NL, B <: LinkBranch{RT, NL}}
    push!(node.other, branch)
end

import Phylo.API._removeoutbound!
function _removeoutbound!(node::LinkNode{RT, NL}, branch::B) where
         {RT <: Rooted, NL, B <: LinkBranch{RT, NL}}
    branch ∈ _getoutbounds(node) ? filter!(p -> p !== branch, node.other) :
         error("Node does not have outbound connection to branch $branch")
end

import Phylo.API._getconnections
_getconnections(node::LinkNode{Unrooted, NL}) where
    {NL} = node.other

import Phylo.API._addconnection!
function _addconnection!(node::LinkNode{Unrooted, NL}, branch::B) where
    {NL, B <: LinkBranch{Unrooted, NL}}
    push!(node.other, branch)
end

import Phylo.API._removeconnection!
function _removeconnection!(node::LinkNode{Unrooted, NL}, branch::B) where
    {NL, B <: LinkBranch{Unrooted, NL}}
    outbound ∈ _getoutbounds(node) ? filter!(p -> p !== branch, node.other) :
         error("Node does not have a connection to branch $branch")
end

# LinkTree methods
const TREENAME = "Tree"
import Phylo.API._treenametype
_treenametype(::Type{T}) where
    {RT, NL, N, B, T <: LinkTree{RT, NL, N, B}} = typeof(TREENAME)
import Phylo.API._gettreename
_gettreename(tree::LinkTree{RT, NL, N, B, TD}) where
    {RT, NL, N, B, TD} = ismissing(tree.name) ? TREENAME : tree.name

import Phylo.API._addnode!
function Phylo.API._addnode!(tree::LinkTree{RT, NL, N, B, TD},
    node::N, name::NL) where {RT, NL, N, B, TD}
    tree.nodes[name] = node
    push!(tree.roots, node)
    isvalid = missing
    return name
end

import Phylo.API._removenode!
function Phylo.API._removenode!(tree::LinkTree{RT, NL, N, B, TD},
                                node::N) where {RT, NL, N, B, TD}
    delete!(tree.nodes, node.name)
    filter!(n -> n !== node, tree.roots)
    isvalid = missing
end

import Phylo.API._deletenode!
_deletenode!(tree::LinkTree{RT, NL, N, B, TD}, name::NL) where
    {RT, NL, N, B, TD} = _deletenode!(tree, tree.nodes[name])
function _deletenode!(tree::LinkTree{RT, NL, N, B, TD}, node::N) where
    {RT <: Rooted, NL, N, B, TD}
    _hasinbound(node) && _deletebranch!(tree, node.inbound)
    while _outdegree(node) > 0
        _deletebranch!(tree, first(node.other))
    end
    _removenode!(tree, node)
end
function _deletenode!(tree::LinkTree{Unrooted, NL, N, B, TD}, branch::B) where
    {NL, N, B, TD}
    while _degree(node) > 0
        _deletebranch!(tree, first(node.other))
    end
    _removenode!(tree, node)
end

import Phylo.API._getroots
_getroots(tree::LinkTree{RT, NL, N, B, TD}) where
    {RT <: Rooted, NL, N, B, TD} = tree.roots
_getroots(tree::LinkTree{Unrooted, NL, N, B, TD}) where
    {NL, N, B, TD} = error("Unrooted trees do not have roots")

import Phylo.API._getnodes
_getnodes(tree::LinkTree{RT, NL, N, B, TD}) where
    {RT, NL, N, B, TD} = values(tree.nodes)

import Phylo.API._getnodenames
_getnodenames(tree::LinkTree{RT, NL, N, B, TD}) where
    {RT, NL, N, B, TD} = keys(tree.nodes)

import Phylo.API._getnodename
_getnodename(::LinkTree{RT, NL, N, B, TD}, node::N) where
    {RT, NL, N, B, TD} = node.name

import Phylo.API._hasnode
_hasnode(tree::LinkTree{RT, NL, N, B, TD}, name::NL) where
    {RT, NL, N, B, TD} = haskey(tree.nodes, name)
_hasnode(tree::LinkTree{RT, NL, N, B, TD}, node::N) where
    {RT, NL, N, B, TD} = node ∈ values(tree.nodes)

import Phylo.API._getnode
_getnode(tree::LinkTree{RT, NL, N, B, TD}, name::NL) where
    {RT, NL, N, B, TD} = tree.nodes[name]

import Phylo.API._getbranches
_getbranches(tree::LinkTree{RT, NL, N, B, TD}) where
    {RT, NL, N, B, TD} = skipmissing(tree.branches)

import Phylo.API._getbranchnames
_getbranchnames(tree::LinkTree{RT, NL, N, B, TD}) where
    {RT, NL, N, B, TD} =
    filter(i -> !ismissing(tree.branches[i]), 1:length(tree.branches))

import Phylo.API._getbranchname
Phylo.API._getbranchname(tree::LinkTree{RT, NL, N, B, TD}, id::Int) where
    {RT, NL, N, B, TD} = id
Phylo.API._getbranchname(tree::LinkTree{RT, NL, N, B, TD}, branch::B) where
    {RT, NL, N, B, TD} = first(indexin([branch], tree.branches))

import Phylo.API._hasbranch
_hasbranch(tree::LinkTree{RT, NL, N, B, TD}, id) where
    {RT, NL, N, B, TD} =
    1 ≤ id ≤ length(tree.branches) && !ismissing(tree.branches[id])
_hasbranch(tree::LinkTree{RT, NL, N, B, TD}, branch::B) where
    {RT, NL, N, B, TD} = branch ∈ tree.branches

import Phylo.API._getbranch
_getbranch(tree::LinkTree{RT, NL, N, B, TD}, id::Int) where
    {RT, NL, N, B, TD} = tree.branches[id]

import Phylo.API._createbranch!
function _createbranch!(tree::LinkTree{RT, NL, N, B, TD},
                        from::LinkNode{RT, NL},
                        to::LinkNode{RT, NL},
                        length::Float64 = NaN,
                        data::Dict{String, Any} = Dict{String, Any}()) where
    {RT <: Rooted, NL, N, B, TD}
    branch = LinkBranch(from, to, length, data)
    _addoutbound!(from, branch)
    _addinbound!(to, branch)
    return _addbranch!(tree, branch)
end
function _createbranch!(tree::LinkTree{Unrooted, NL, N, B, TD},
                        from::LinkNode{Unrooted, NL},
                        to::LinkNode{Unrooted, NL},
                        length::Float64 = NaN,
                        data::Dict{String, Any} = Dict{String, Any}()) where
    {NL, N, B, TD}
    branch = LinkBranch(tree, from, to, length, data)
    _addconnection!(from, branch)
    _addconnection!(to, branch)
    return _addbranch!(tree, branch)
end

import Phylo.API._addbranch!
function _addbranch!(tree::LinkTree{RT, NL, N, B, TD},
                     branch::B) where {RT <: Rooted, NL, N, B, TD}
    id = _newbranchlabel(tree)
    len = length(tree.branches)
    if id <= len
        tree.branches[id] = branch
    elseif id == len + 1
        push!(tree.branches, branch)
    else
        error("Adding branch beyond end of branch vector")
    end

    filter!(n -> n !== _dst(branch), tree.roots)
    isvalid = missing
    return id
end
function _addbranch!(tree::LinkTree{Unrooted, NL, N, B, TD},
                     branch::B, id::Int) where {NL, N, B, TD}
    id = _newbranchlabel(tree)
    tree.branches[id] = branch
    isvalid = missing
    return id
end

import Phylo.API._removebranch!
function Phylo.API._removebranch!(tree::LinkTree{RT, NL, N, B, TD},
                                  branch::B,
                                  id::Int = first(indexin([branch],
                                                  tree.branches))) where
    {RT <: Rooted, NL, N, B, TD}
    tree.branches[id] = missing
    push!(tree.roots, _dst(branch))
    isvalid = missing
end
function Phylo.API._removebranch!(tree::LinkTree{Unrooted, NL, N, B, TD},
                                  branch::B,
                                  id::Int = first(indexin([branch],
                                                  tree.branches))) where
    {NL, N, B, TD}
    tree.branches[id] = missing
    isvalid = missing
end

import Phylo.API._deletebranch!
_deletebranch!(tree::LinkTree{RT, NL, N, B, TD}, id::Int) where
    {RT, NL, N, B, TD} = _deletebranch!(tree, tree.branches[id], id)
function _deletebranch!(tree::LinkTree{RT, NL, N, B, TD}, branch::B,
                        id::Int = first(indexin([branch], tree.branches))) where
    {RT <: Rooted, NL, N, B, TD}
    _removeoutbound!(branch.inout[1], branch)
    _removeinbound!(branch.inout[2], branch)
    _removebranch!(tree, branch, id)
end
function _deletebranch!(tree::LinkTree{Unrooted, NL, N, B, TD}, branch::B,
                        id::Int = first(indexin([branch], tree.branches))) where
    {NL, N, B, TD}
    _removeconnection!(branch.inout[1], branch)
    _removeconnection!(branch.inout[2], branch)
    _removebranch!(tree, branch, id)
end

import Phylo.API._createnode!
function _createnode!(tree::LinkTree{RT, NL, N, B, TD}, nodename::NL,
                      data::Dict{String, Any} = Dict{String, Any}()) where
    {RT <: Rooted, NL, N, B, TD}
    node = LinkNode(tree, nodename, data)
    return _addnode!(tree, node, nodename)
end

import Phylo.API._getleafinfo
_getleafinfo(tree::LinkTree{RT, NL, N, B, TD}) where
    {RT, NL, N, B, TD} = tree.tipdata

import Phylo.API._noderecordtype
Phylo.API._noderecordtype(::Type{LinkTree{RT, NL, N, B, TD}}) where
    {RT, NL, N, B, TD} = _noderecordtype(N)

import Phylo.API._getnoderecord
Phylo.API._getnoderecord(tree::LinkTree{RT, NL, N, B, TD}, name::NL) where
    {RT, NL, N, B, TD} = _getnoderecord(_getnode(tree, name))
Phylo.API._getnoderecord(tree::LinkTree{RT, NL, N, B, TD}, node::N) where
    {RT, NL, N, B, TD} = _getnoderecord(node)

import Phylo.API._branchrecordtype
Phylo.API._branchrecordtype(::Type{LinkTree{RT, NL, N, B, TD}}) where
    {RT, NL, N, B, TD} = _branchrecordtype(N)

import Phylo.API._getbranchinfo
Phylo.API._getbranchinfo(tree::LinkTree{RT, NL, N, B, TD}, id::Int) where
    {RT, NL, N, B, TD} = _getbranchinfo(_getbranch(tree, id))
Phylo.API._getbranchinfo(tree::LinkTree{RT, NL, N, B, TD}, branch::B) where
    {RT, NL, N, B, TD} = _getbranchinfo(branch)
