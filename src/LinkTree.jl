using Missings
using Compat
using Compat: @warn
using SimpleTraits

newempty(::Type{Data}) where Data = Data()

struct LinkBranch{RT, NL, Data} <: AbstractBranch{RT, NL}
    name::Int
    inout::Tuple{AbstractNode{RT, NL}, AbstractNode{RT, NL}}
    length::Float64
    data::Data
end

mutable struct LinkNode{RT, NL, Data,
                        B <: AbstractBranch{RT, NL}} <: AbstractNode{RT, NL}
    name::NL
    inbound::Union{B, Missing}
    other::Vector{B}
    data::Data
    LinkNode{RT, NL, Data, B}(name::NL, data::Data = newempty(Data)) where
    {RT, NL, Data, B} = new{RT, NL, Data, B}(name, missing, Vector{B}(), data)

end

import Phylo.API: _prefernodeobjects
_prefernodeobjects(::Type{<:LinkNode}) = true

mutable struct LinkTree{RT, NL, N <: LinkNode{RT, NL},
                        B <: LinkBranch{RT, NL}, TD} <:
                            AbstractTree{OneTree, RT, NL, N, B}
    name::Union{String, Missing}
    nodedict::Dict{NL, Int}
    roots::Vector{N}
    nodes::Vector{Union{N, Missing}}
    branches::Vector{Union{B, Missing}}
    data::Dict{String, Any}
    tipdata::TD
    rootheight::Float64
    isvalid::Union{Bool, Missing}
    cache::Dict{TraversalOrder, Vector{N}}
    
    function LinkTree{RT, NL, N, B, TD}(tipnames::Vector{NL} = NL[];
                                        treename::Union{String, Missing} = missing,
                                        tipdata::TD = newempty(TD),
                                        rootheight = NaN) where
        {RT, NL, N, B, TD}
        tree = new{RT, NL, N, B, TD}(treename, Dict{NL, N}(), Vector{N}(),
                                     Vector{Union{N, Missing}}(), Vector{B}(),
                                     Dict{String, Any}(),
                                     tipdata, rootheight, missing,
                                     Dict{TraversalOrder, Vector{N}}())
        if !isempty(tipnames)
            createnodes!(tree, tipnames)
        elseif !isempty(tipdata)
            tips = unique(keys(tipdata))
            createnodes!(tree, tips)
        end
        return tree
    end
end
function LinkTree{RT, NL, N, B, TD}(leafinfos::TD) where {RT, NL, N, B, TD}
    leafnames = unique(info[1] for info in getiterator(leafinfos))
    return LinkTree{RT, NL, N, B, TD}(leafnames; tipdata = leafinfos)
end

const LB{RT} = LinkBranch{RT, String, Dict{String, Any}}
const LN{RT} = LinkNode{RT, String, Dict{String, Any}, LB{RT}}
const LT{RT, TD} = LinkTree{RT, String, LN{RT}, LB{RT}, TD}
const RootedTree = LT{OneRoot, Dict{String, Any}}
const ManyRootTree = LT{ManyRoots, Dict{String, Any}}
const UnrootedTree = LT{Unrooted, Dict{String, Any}}

# LinkBranch methods
function LinkBranch(name::Int,
                    from::LinkNode{RT, NL},
                    to::LinkNode{RT, NL},
                    len::Float64 = NaN,
                    data::Data = nothing) where {RT, NL, Data}
    len >= 0.0 || isnan(len) ||
        error("Branch length must be positive or NaN (no recorded length), not $len")
    return LinkBranch{RT, NL, Data}(name, (from, to), len, data)
end

import Phylo.API: _preferbranchobjects
_preferbranchobjects(::Type{<:LinkBranch}) = true

import Phylo.API: _src
_src(::AbstractTree, branch::LinkBranch{<:Rooted}) = branch.inout[1]
import Phylo.API: _dst
_dst(::AbstractTree, branch::LinkBranch{<:Rooted}) = branch.inout[2]
import Phylo.API: _conn
function _conn(::AbstractTree, branch::LinkBranch{RT, NL, D},
               exclude::AbstractNode{RT, NL}) where {RT, NL, D}
    return exclude ≡ branch.inout[1] ? branch.inout[2] :
        (exclude ≡ branch.inout[2] ? branch.inout[1] :
         error("Branch $(branch.name) not connected to $(exclude.name)"))
end
import Phylo.API: _getlength
_getlength(::AbstractTree, branch::LinkBranch) = branch.length
import Phylo.API: _getbranchname
_getbranchname(tree::LinkTree, id::Int) = id
function _getbranchname(tree::LinkTree{RT, NL, N, B}, branch::B) where
    {RT, NL, N, B}
    id = branch.name
    branch ≡ tree.branches[id] || error("Branch $id not in tree")
    return id
end
import Phylo.API: _branchdatatype
_branchdatatype(::Type{LinkBranch{RT, NL, Data}}) where {RT, NL, Data} = Data
import Phylo.API: _getbranchdata
_getbranchdata(::AbstractTree, branch::LinkBranch) = branch.data
import Phylo.API: _setbranchdata!
_setbranchdata!(::AbstractTree, branch::LinkBranch{RT, NL, Data},
                data::Data) where {RT, NL, Data} = (branch.data = data)

# LinkNode methods
function LinkNode(tree::T, name::NL,
                  data::Dict{String, Any} = Dict{String, Any}()) where
    {RT, NL, N <: LinkNode, B, T <: AbstractTree{OneTree, RT, NL, N, B}}
    return N(name, missing, Vector{B}(), data)
end

import Phylo.API: _nodedatatype
_nodedatatype(::Type{LinkNode{RT, NL, Data, B}}) where {RT, NL, B, Data} = Data

import Phylo.API: _getnodedata
_getnodedata(::LinkTree, node::LinkNode) = node.data

import Phylo.API: _setnodedata!
_setnodedata!(::LinkTree, node::LinkNode{RT, NL, Data, B},
              data::Data) where {RT, NL, B, Data} = (node.data = data)

import Phylo.API: _hasinbound
_hasinbound(::LinkTree, node::LinkNode{<: Rooted}) = !ismissing(node.inbound)

import Phylo.API: _degree
_degree(::LinkTree, node::LinkNode{Unrooted}) = length(node.other)

import Phylo.API: _getinbound
_getinbound(::LinkTree, node::LinkNode) = node.inbound

import Phylo.API: _addinbound!
function _addinbound!(tree::LinkTree,
                      node::LinkNode{<: Rooted, NL, Data, B},
                      branch::B) where {NL, Data, B}
    _hasinbound(tree, node) &&
        error("LinkNode already has an inbound connection")
    node.inbound = branch
    filter!(n -> n ≠ node, tree.roots)
    tree.isvalid = missing
end

import Phylo.API: _removeinbound!
function _removeinbound!(tree::LinkTree, node::LinkNode{<: Rooted, NL, Data, B},
                         branch::B) where {NL, Data, B}
    _hasinbound(tree, node) || error("Node has no inbound connection")
    node.inbound ≡ branch ||
        error("Node has no inbound connection from branch $inbound")
    node.inbound = missing
    push!(tree.roots, node)
    tree.isvalid = missing
end

import Phylo.API: _getoutbounds
_getoutbounds(::LinkTree, node::LinkNode{<: Rooted}) = node.other

import Phylo.API: _addoutbound!
_addoutbound!(::LinkTree,
              node::LinkNode{<: Rooted, NL, Data, B},
              branch::B) where {NL, Data, B} = push!(node.other, branch)

import Phylo.API: _removeoutbound!
function _removeoutbound!(tree::LinkTree,
                          node::LinkNode{<: Rooted, NL, Data, B},
                          branch::B) where {NL, Data, B}
    branch ∈ _getoutbounds(tree, node) ? filter!(p -> p ≢ branch, node.other) :
        error("Node does not have outbound connection to branch $branch")
end

import Phylo.API: _getconnections
_getconnections(::AbstractTree, node::LinkNode{Unrooted}) = node.other

import Phylo.API: _addconnection!
function _addconnection!(::AbstractTree, node::LinkNode{Unrooted, NL, Data, B},
                         branch::B) where {NL, Data, B}
    push!(node.other, branch)
end

import Phylo.API: _removeconnection!
function _removeconnection!(tree::AbstractTree,
                            node::LinkNode{Unrooted, NL}, branch::B) where
    {NL, B <: LinkBranch{Unrooted, NL}}
    outbound ∈ _getoutbounds(tree, node) ?
        filter!(p -> p ≢ branch, node.other) :
        error("Node does not have a connection to branch $branch")
end

# LinkTree methods
const TREENAME = "Tree"
const NODENAME = "Node"

import Phylo.API: _validate!
function _validate!(tree::LinkTree{RT, NL, N, B, TD}) where {RT, NL, N, B, TD}
    tree.isvalid = true
    if TD <: Dict
        if !isempty(tree.tipdata)
            tree.isvalid &= (Set(keys(tree.tipdata)) == Set(getleafnames(tree)))
        end
    end
    nr = nroots(tree)
    if RT == OneRoot
        if nr != 1
            @warn "Wrong number of roots for $RT tree ($nr)"
            tree.isvalid = false
        end
    elseif RT == ManyRoots
        if nr < 1
            @warn "Wrong number of roots for $RT tree ($nr)"
            tree.isvalid = false
        end
    end
    return tree.isvalid
end

import Phylo.API: _treenametype
_treenametype(::Type{<: LinkTree}) = String
import Phylo.API: _gettreename
_gettreename(tree::LinkTree) = ismissing(tree.name) ? TREENAME : tree.name
import Phylo.API: _newnodelabel
function _newnodelabel(tree::LinkTree{RT, String}) where RT
    id = length(tree.nodes) + 1
    name = NODENAME * " $id"
    return haskey(tree.nodedict, name) ?
        _newlabel(collect(getnodenames(tree)), NODENAME) : name
end
import Phylo.API: _newbranchlabel
_newbranchlabel(tree::LinkTree) = length(tree.branches) + 1
import Phylo.API: _deletenode!
function _deletenode!(tree::LinkTree{RT, NL, N}, node::N) where {RT, NL, N}
    (haskey(tree.nodedict, node.name) &&
     tree.nodes[tree.nodedict[node.name]] ≡ node) ||
     error("Node $(node.name) is not in tree $(_gettreename(tree)), cannot be deleted")
    # Does nothing for unrooted tree, as inbound is never set
    !ismissing(node.inbound) && _deletebranch!(tree, node.inbound)
    # Delete outbound connections of rooted tree, all connections of unrooted
    while _degree(tree, node) > 0
        _deletebranch!(tree, first(node.other))
    end
    tree.nodes[tree.nodedict[node.name]] = missing
    delete!(tree.nodedict, node.name)
    filter!(n -> n ≢ node, tree.roots)
    tree.isvalid = missing
    return true
end

import Phylo.API: _getroots
_getroots(tree::LinkTree{<: Rooted}) = tree.roots
_getroots(tree::LinkTree{Unrooted}) = error("Unrooted trees do not have roots")

import Phylo.API: _getnodes
_getnodes(tree::LinkTree) = skipmissing(tree.nodes)

import Phylo.API: _getnodenames
_getnodenames(tree::LinkTree) = keys(tree.nodedict)

import Phylo.API: _getnodename
_getnodename(::LinkTree{RT, NL, N}, node::N) where {RT, NL, N} = node.name

import Phylo.API: _hasnode
_hasnode(tree::LinkTree{RT, NL}, name::NL) where {RT, NL} =
    haskey(tree.nodedict, name)
_hasnode(tree::LinkTree{RT, NL, N}, node::N) where {RT, NL, N} =
    haskey(tree.nodedict, node.name)

import Phylo.API: _getnode
@traitfn _getnode(tree::T,
                  node::NL) where {RT, NL,
                                   T <: LinkTree{RT, NL};
                                   PreferNodeObjects{T}} =
                                       tree.nodes[tree.nodedict[node]]

import Phylo.API: _getbranches
_getbranches(tree::LinkTree) = skipmissing(tree.branches)

import Phylo.API: _getbranchnames
_getbranchnames(tree::LinkTree) =
    [getbranchname(tree, b) for b in tree.branches if !ismissing(b)]

import Phylo.API: _hasbranch
_hasbranch(tree::LinkTree, id::Int) =
    1 ≤ id ≤ length(tree.branches) && !ismissing(tree.branches[id])
_hasbranch(tree::LinkTree{RT, NL, N, B}, branch::B) where
{RT, NL, N <: LinkNode{RT, NL}, B <: LinkBranch{RT, NL}} =
    branch ∈ skipmissing(tree.branches)

import Phylo.API: _getbranch
_getbranch(tree::LinkTree, id::Int) = tree.branches[id]

import Phylo.API: _createbranch!
function _createbranch!(tree::LinkTree{RT, NL, N, B},
                        from::N, to::N, len::Float64 = NaN,
                        name::Union{Int, Missing} = missing,
                        data::Data = newempty(Data)) where
    {RT <: Rooted, NL, N, Data, B <: LinkBranch{RT, NL, Data}}
    id = ismissing(name) ? length(tree.branches) + 1 : name
    branch = LinkBranch(id, from, to, len, data)
    if ismissing(name)
        push!(tree.branches, branch)
    else
        l = length(tree.branches)
        if id > l
            resize!(tree.branches, id)
            tree.branches[l+1:id-1] .= missing
        else
            ismissing(tree.branches[id]) || error("Branch $name already exists")
        end
        tree.branches[id] = branch
    end
    _addoutbound!(tree, from, branch)
    _addinbound!(tree, to, branch)
    tree.isvalid = missing
    return branch
end
function _createbranch!(tree::LinkTree{Unrooted, NL, N, B},
                        from::N, to::N, len::Float64 = NaN,
                        name::Union{Int, Missing} = missing,
                        data::Data = newempty(Data)) where
    {NL, N, Data, B <: LinkBranch{Unrooted, NL, Data}}
    id = ismissing(name) ? length(tree.branches) + 1 : name
    branch = LinkBranch(id, tree, from, to, len, data)
    push!(tree.branches, branch)
    _addconnection!(tree, from, branch)
    _addconnection!(tree, to, branch)
    tree.isvalid = missing
    return branch
end

import Phylo.API: _deletebranch!
function _deletebranch!(tree::LinkTree{<: Rooted}, id::Int)
    branch = tree.branches[id]
    _removeoutbound!(tree, branch.inout[1], branch)
    _removeinbound!(tree, branch.inout[2], branch)
    tree.branches[id] = missing
    tree.isvalid = missing
    return true
end
function _deletebranch!(tree::LinkTree{Unrooted}, id::Int)
    branch = tree.branches[id]
    _removeconnection!(tree, branch.inout[1], branch)
    _removeconnection!(tree, branch.inout[2], branch)
    tree.branches[id] = missing
    tree.isvalid = missing
    return true
end
function _deletebranch!(tree::LinkTree{<: Rooted, NL, N, B}, branch::B) where
    {NL, N, B}
    _removeoutbound!(tree, branch.inout[1], branch)
    _removeinbound!(tree, branch.inout[2], branch)
    id = _getbranchname(tree, branch)
    tree.branches[id] = missing
    tree.isvalid = missing
    return true
end
function _deletebranch!(tree::LinkTree{Unrooted, NL, N, B}, branch::B) where
    {NL, N, B}
    _removeconnection!(tree, branch.inout[1], branch)
    _removeconnection!(tree, branch.inout[2], branch)
    id = _getbranchname(tree, branch)
    tree.branches[id] = missing
    tree.isvalid = missing
    return true
end

import Phylo.API: _createnode!
function _createnode!(tree::LinkTree{RT, NL, N, B}, name::Union{NL, Missing},
                      data::Data = newempty(Data)) where
    {RT <: Rooted, NL, Data, B, N <: LinkNode{RT, NL, Data, B}}
    nodename = ismissing(name) ? _newnodelabel(tree) : name
    _hasnode(tree, nodename) && error("Node $nodename already exists in tree.")
    node = LinkNode{RT, NL, Data, B}(nodename, data)
    id = length(tree.nodes) + 1
    push!(tree.nodes, node)
    tree.nodedict[nodename] = id
    push!(tree.roots, node)
    return node
end

import Phylo.API: _leafinfotype
leafinfotype(::Type{LinkTree{RT, NL, N, B, TD}}) where {RT, NL, N, B, TD} = TD

import Phylo.API: _getleafinfo
_getleafinfo(tree::LinkTree) = tree.tipdata

import Phylo.API: _setleafinfo!
_setleafinfo!(tree::LinkTree{RT, NL, N, B, TD}, leafinfo::TD) where
{RT, NL, N, B, TD} = tree.tipdata = leafinfo

import Phylo.API: _nodedatatype
_nodedatatype(::Type{LinkTree{RT, NL, N}}) where
{RT, NL, Data, N <: LinkNode{RT, NL, Data}} = _nodedatatype(N)

import Phylo.API: _branchdatatype
_branchdatatype(::Type{LinkTree{RT, NL, N, B, TD}}) where {RT, NL, N, B, TD} =
    _branchdatatype(B)

import Phylo.API: _nnodes
_nnodes(tree::LinkTree) = count((!)∘ismissing, tree.nodes)
import Phylo.API: _nbranches
_nbranches(tree::LinkTree) = count((!)∘ismissing, tree.branches)
