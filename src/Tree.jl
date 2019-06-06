using DataStructures
using Compat
using Compat: @warn
using IterableTables: getiterator
using DataFrames

abstract type AbstractBranchTree{RT, N, B, LI, ND} <:
    AbstractTree{OneTree, RT, String, N, B}
end

import Phylo.API: _hasinbound
function _hasinbound(tree::AbstractBranchTree{<: Rooted}, name::String) 
    return tree.nodes[name].inbound !== nothing
end

import Phylo.API: _getinbound
function _getinbound(tree::AbstractBranchTree{<: Rooted}, name::String)
    return tree.nodes[name].inbound
end

import Phylo.API: _addinbound!
function _addinbound!(tree::AbstractBranchTree{<: Rooted, N, B},
                      name::String,
                      inbound::B) where {N, B <: Branch}
    _hasinbound(tree, name) &&
        error("Node already has an inbound connection")
    tree.nodes[name].inbound = inbound
end

import Phylo.API: _removeinbound!
function _removeinbound!(tree::AbstractBranchTree{<: Rooted, N, B},
                         name::String,
                         inbound::B) where {N, B <: Branch}
    _hasinbound(tree, name) || error("Node has no inbound connection")
    node = tree.nodes[name]
    node.inbound === inbound ||
        error("Node has no inbound connection from branch $inbound")
    node.inbound = nothing
end

import Phylo.API: _leafinfotype
_leafinfotype(::Type{T}) where {RT, N, B, LI,
                                T <: AbstractBranchTree{RT, N, B, LI}} = LI

import Phylo.API: _nodedatatype
_nodedatatype(::Type{<: AbstractBranchTree{RT, N, B, LI, ND}}) where {RT, N, B, LI, ND} = ND

import Phylo.API: _branchdatatype
_branchdatatype(::Type{<: AbstractBranchTree}) = Nothing

import Phylo.API: _getnodes
_getnodes(bt::AbstractBranchTree) = keys(bt.nodes)

import Phylo.API: _getnodenames
_getnodenames(bt::AbstractBranchTree) = keys(bt.nodes)

import Phylo.API: _getbranches
_getbranches(bt::AbstractBranchTree) = values(bt.branches)

import Phylo.API: _getbranchnames
_getbranchnames(bt::AbstractBranchTree) = keys(bt.branches)

import Phylo.API: _getleafinfo
_getleafinfo(bt::AbstractBranchTree) = bt.leafinfos
function _getleafinfo(bt::AbstractBranchTree, leafname)
    return [info for info in getiterator(bt.leafinfos) if info[1] == leafname]
end

import Phylo.API: _setleafinfo!
function _setleafinfo!(bt::AbstractBranchTree, info)
    bt.leafinfos = info
end

import Phylo.API: _resetleaves!
function _resetleaves!(bt::AbstractBranchTree)
    # bt.leafinfos = empty!(bt.leafinfos)
    return true
end

import Phylo.API: _getnodedata
function _getnodedata(bt::AbstractBranchTree, name::String)
    return bt.nodedata[name]
end

import Phylo.API: _setnodedata!
function _setnodedata!(bt::AbstractBranchTree{RT, N, B, LI, ND},
                       name::String, value::ND) where {RT, N, B, LI, ND}
    bt.nodedata[name] = value
end

import Phylo.API: _hasnode
_hasnode(tree::AbstractBranchTree, name::String) = haskey(tree.nodes, name)

import Phylo.API: _getnode
@traitfn _getnode(tree::T,
                  nodename::String) where {T <: AbstractBranchTree;
                                           PreferNodeObjects{T}} =
                                               tree.nodes[name]

import Phylo.API: _createnode!
function _createnode!(tree::AbstractBranchTree{RT, N, B, LI, ND},
                      name::Union{String, Missing},
                      data::ND = ND()) where {RT, N, B, LI, ND}
    nodename = ismissing(name) ? _newnodelabel(tree) : name
    _hasnode(tree, nodename) && error("Node $name already present in tree")
    tree.nodes[nodename] = N(nodename)
    _setnodedata!(tree, nodename, data)
    return nodename
end

import Phylo.API: _deletenode!
function _deletenode!(tree::AbstractBranchTree, node::AbstractNode)
    name = getnodename(tree, node)
    if _hasinbound(tree, node)
        deletebranch!(tree, _getinbound(tree, node))
    end
    for b in _getoutbounds(tree, node)
        deletebranch!(tree, b)
    end
    delete!(tree.nodes, name)
    delete!(tree.nodedata, name)
    return true
end
function _deletenode!(tree::AbstractBranchTree, name::String)
    node = getnode(tree, name)
    if _hasinbound(tree, node)
        deletebranch!(tree, _getinbound(tree, node))
    end
    for b in _getoutbounds(tree, node)
        deletebranch!(tree, b)
    end
    delete!(tree.nodes, name)
    delete!(tree.nodedata, name)
    return true
end

import Phylo.API: _hasbranch
_hasbranch(tree::AbstractBranchTree, name::Int) = haskey(tree.branches, name)

import Phylo.API: _getbranch
_getbranch(tree::AbstractBranchTree, name::Int) = tree.branches[name]

import Phylo.API: _createbranch!
function _createbranch!(tree::AbstractBranchTree{RT}, source, destination,
                        len::Union{Float64, Missing}, name::Union{Int, Missing},
                        data::Nothing = nothing) where RT
    # Add the new branch
    id = ismissing(name) ? _newbranchlabel(tree) : name
    branch = Branch{RT}(id, _getnodename(tree, source),
                        _getnodename(tree, destination), len)
    tree.branches[id] = branch

    # Update the associated source and destination nodes
    _addoutbound!(tree, _getnode(tree, source), branch)
    _addinbound!(tree, _getnode(tree, destination), branch)

    # Return updated tree
    return branch
end

import Phylo.API: _deletebranch!
function _deletebranch!(tree::AbstractBranchTree, name::Int)
    # Find the branch
    branch = _getbranch(tree, name)
    # Remove branch reference from its source node
    _removeoutbound!(tree, _getnode(tree, _src(tree, branch)), branch)
    # Remove branch reference from its destination node
    _removeinbound!(tree, _getnode(tree, _dst(tree, branch)), branch)
    # Remove branch itself
    delete!(tree.branches, name)
    # Return the branch name
    return true
end
function _deletebranch!(tree::AbstractBranchTree, branch::Branch)
    # Find the branch
    name = _getbranchname(tree, branch)
    # Remove branch reference from its source node
    _removeoutbound!(tree, _getnode(tree, _src(tree, branch)), branch)
    # Remove branch reference from its destination node
    _removeinbound!(tree, _getnode(tree, _dst(tree, branch)), branch)
    # Remove branch itself
    delete!(tree.branches, name)
    # Return the branch name
    return true
end

import Phylo.API: _validate!
function _validate!(tree::TREE) where {RT, TREE <: AbstractBranchTree{RT}}
    if _leafinfotype(TREE) != Nothing && length(getiterator(tree.leafinfos)) > 0
        if Set(info[1] for info in getiterator(tree.leafinfos)) !=
            Set(_getleafnames(tree, anyorder))
            @warn "LeafInfo names do not match actual leaves of tree"
            return false
        end
    end

    if Set(keys(tree.nodedata)) != Set(_getnodenames(tree))
        @warn "Node names do not match node records of tree"
        return false
    end

    nr = nroots(tree)
    if RT == OneRoot
        if nr != 1
            @warn "Wrong number of roots for $RT tree ($nr)"
            return false
        end
    elseif RT == ManyRoots
        if nr < 1
            @warn "Wrong number of roots for $RT tree ($nr)"
            return false
        end
    end

    return true
end

import Phylo.API: _hasrootheight
function _hasrootheight(tree::AbstractBranchTree)
    return !ismissing(tree.rootheight)
end

import Phylo.API: _getrootheight
function _getrootheight(tree::AbstractBranchTree)
    return tree.rootheight
end

import Phylo.API: _setrootheight!
function _setrootheight!(tree::AbstractBranchTree, height::Float64)
    tree.rootheight = height
    return height
end

import Phylo.API: _clearrootheight!
function _clearrootheight!(tree::AbstractBranchTree)
    tree.rootheight = missing
end

"""
    BinaryTree

Binary phylogenetic tree object with known leaves and per node data
"""
mutable struct BinaryTree{RT <: Rooted, LI, ND} <:
    AbstractBranchTree{RT, BinaryNode{RT, String, Branch{RT, String}},
                       Branch{RT, String}, LI, ND}
    nodes::OrderedDict{String, BinaryNode{RT, String, Branch{RT, String}}}
    branches::Dict{Int, Branch{RT, String}}
    leafinfos::LI
    nodedata::OrderedDict{String, ND}
    rootheight::Union{Float64, Missing}
end

function BinaryTree(lt::BinaryTree{RT, LI, ND};
                    copyinfo=false) where {RT, LI, ND}
    validate!(lt) || error("Tree to copy is not valid")
    leafnames = getleafnames(lt)
    # Leaf records are conserved across trees, as could be invariant?
    leafinfos = copyinfo ? deepcopy(lt.leafinfos) : lt.leafinfos
    nodes = OrderedDict(leaf => BinaryNode{RT, String, Branch{RT, String}}(leaf)
                        for leaf in leafnames)
    branches = Dict{Int, Branch{RT, String}}()
    nodedata = OrderedDict(leaf => ND() for leaf in leafnames)
    return BinaryTree{RT, LI, ND}(nodes, branches, leafinfos, nodedata,
                                  lt.rootheight)
end

function BinaryTree{RT, LI, ND}(leafnames::Vector{String},
                                treetype::Type{BinaryTree{RT, LI, ND}} =
                                BinaryTree{RT, LI, ND};
                                rootheight::Union{Float64, Missing} =
                                missing) where {RT, LI, ND}
    nodes = OrderedDict(leaf => BinaryNode{RT, String, Branch{RT, String}}(leaf)
                        for leaf in leafnames)
    leafinfos = LI()
    nodedata = OrderedDict(leaf => ND() for leaf in leafnames)
    return BinaryTree{RT, LI, ND}(nodes, Dict{Int, Branch{RT, String}}(),
                                  leafinfos, nodedata, rootheight)
end

function BinaryTree{RT, LI, ND}(numleaves::Int = 0,
                                treetype::Type{BinaryTree{RT, LI, ND}} =
                                BinaryTree{RT, LI, ND};
                                rootheight::Union{Float64, Missing} =
                                missing) where {RT, LI, ND}
    leafnames = ["Leaf $num" for num in Base.OneTo(numleaves)]
    return BinaryTree{RT, LI, ND}(leafnames, treetype; rootheight = rootheight)
end

function BinaryTree{RT, LI, ND}(leafinfos::LI,
                                treetype::Type{BinaryTree{RT, LI, ND}} =
                                BinaryTree{RT, LI, ND};
                                rootheight::Union{Float64, Missing} =
                                missing) where {RT, LI, ND}
    leafnames = unique(info[1] for info in getiterator(leafinfos))
    nodes = OrderedDict(leaf => BinaryNode{RT, String, Branch{RT, String}}(leaf)
                        for leaf in leafnames)
    branches = Dict{Int, Branch{RT, String}}()
    nodedata = OrderedDict(leaf => ND() for leaf in leafnames)
    return BinaryTree{RT, LI, ND}(nodes, branches, leafinfos,
                                  nodedata, rootheight)
end

BinaryTree{RT}(leafinfos::LI;
               rootheight::Union{Float64, Missing} =
               missing) where {RT <: Rootedness, LI} =
               BinaryTree{RT, LI, Dict{String, Any}}(leafinfos;
                                                     rootheight = rootheight)

BinaryTree(leafinfos::LI;
           rootheight::Union{Float64, Missing} = missing) where LI =
    BinaryTree{OneRoot}(leafinfos; rootheight = rootheight)

"""
    NamedBinaryTree

Binary phylogenetic tree object with known leaves
"""
const NamedBinaryTree = BinaryTree{ManyRoots, DataFrame, Dict{String, Any}}

"""
    PolytomousTree

Phylogenetic tree object with polytomous branching, and known leaves
and per node data
"""
mutable struct PolytomousTree{RT <: Rooted, LI, ND} <:
    AbstractBranchTree{RT, Node{RT, String, Branch{RT, String}},
                       Branch{RT, String}, LI, ND}
    nodes::OrderedDict{String, Node{RT, String, Branch{RT, String}}}
    branches::Dict{Int, Branch{RT, String}}
    leafinfos::LI
    nodedata::OrderedDict{String, ND}
    rootheight::Union{Float64, Missing}
end

function PolytomousTree(lt::PolytomousTree{RT, LI, ND};
                        copyinfo=false) where {RT, LI, ND}
    validate!(lt) || error("Tree to copy is not valid")
    leafnames = getleafnames(lt)
    # Leaf records may be conserved across trees, as could be invariant?
    leafinfos = copyinfo ? deepcopy(lt.leafinfos) : lt.leafinfos
    nodes = OrderedDict(leaf => Node{RT, String, Branch{RT, String}}(leaf)
                        for leaf in leafnames)
    branches = Dict{Int, Branch{RT, String}}()
    nodedata = OrderedDict(leaf => ND() for leaf in leafnames)
    return PolytomousTree{RT, LI, ND}(nodes, branches, leafinfos,
                                      nodedata, lt.rootheight)
end

function PolytomousTree{RT, LI, ND}(leafnames::Vector{String},
                                    treetype::Type{PolytomousTree{RT, LI, ND}} =
                                    PolytomousTree{RT, LI, ND};
                                    rootheight::Union{Float64, Missing} =
                                    missing) where {RT, LI, ND}
    nodes = OrderedDict(leaf => Node{RT, String, Branch{RT, String}}(leaf)
                        for leaf in leafnames)
    leafinfos = LI()
    nodedata = OrderedDict(leaf => ND() for leaf in leafnames)
    return PolytomousTree{RT, LI, ND}(nodes, Dict{Int, Branch{RT, String}}(),
                                      leafinfos, nodedata, rootheight)
end

function PolytomousTree{RT, LI, ND}(numleaves::Int = 0,
                                    treetype::Type{PolytomousTree{RT, LI, ND}} =
                                    PolytomousTree{RT, LI, ND};
                                    rootheight::Union{Float64, Missing} =
                                    missing) where {RT, LI, ND}
    leafnames = ["Leaf $num"  for num in Base.OneTo(numleaves)]
    return PolytomousTree{RT, LI, ND}(leafnames, treetype;
                                      rootheight = rootheight)
end

function PolytomousTree{RT, LI, ND}(leafinfos::LI,
                                    treetype::Type{PolytomousTree{RT, LI, ND}} =
                                    PolytomousTree{RT, LI, ND};
                                    rootheight::Union{Float64, Missing} =
                                    missing) where {RT, LI, ND}
    leafnames = unique(info[1] for info in getiterator(leafinfos))
    nodes = OrderedDict(leaf => Node{RT, String, Branch{RT, String}}(leaf)
                        for leaf in leafnames)
    branches = Dict{Int, Branch{RT, String}}()
    nodedata = OrderedDict(leaf => ND() for leaf in leafnames)
    return PolytomousTree{RT, LI, ND}(nodes, branches, leafinfos,
                                      nodedata, rootheight)
end

PolytomousTree{RT}(leafinfos::LI;
                   rootheight::Union{Float64, Missing} = missing) where {RT, LI} =
    PolytomousTree{RT, LI, Dict{String, Any}}(leafinfos; rootheight = rootheight)

PolytomousTree(leafinfos::LI;
               rootheight::Union{Float64, Missing} = missing) where LI =
    PolytomousTree{OneRoot}(leafinfos; rootheight = rootheight)

"""
    NamedPolytomousTree

Polytomous phylogenetic tree object with known leaves
"""
const NamedTree = NamedPolytomousTree =
    PolytomousTree{ManyRoots, DataFrame, Dict{String, Any}}

import Phylo.API: _outdegree
function _outdegree(tree::BinaryTree{<: Rooted}, name::String)
    node = tree.nodes[name]
    return (node.outbounds[1] ≡ nothing ? 0 : 1) +
        (node.outbounds[2] ≡ nothing ? 0 : 1)
end
function _outdegree(tree::PolytomousTree{<: Rooted}, name::String)
    node = tree.nodes[name]
    return length(node.outbounds)
end

import Phylo.API: _hasoutboundspace
_hasoutboundspace(tree::BinaryTree{<: Rooted}, node::String) =
    _outdegree(tree, node) < 2

import Phylo.API: _getoutbounds
function _getoutbounds(tree::BinaryTree{RT}, name::String) where RT <: Rooted
    node = tree.nodes[name]
    return (node.outbounds[1] ≡ nothing ?
            (node.outbounds[2] ≡ nothing ? Branch{RT, String}[] :
             [node.outbounds[2]]) :
            (node.outbounds[2] ≡ nothing ? [node.outbounds[1]] :
             [node.outbounds[1], node.outbounds[2]]))
end
function _getoutbounds(tree::PolytomousTree{<: Rooted}, name::String)
    return tree.nodes[name].outbounds
end

import Phylo.API: _addoutbound!
function _addoutbound!(tree::BinaryTree{RT}, name::String, branch::B) where
    {RT <: Rooted, N, B <: Branch{RT, String}}
    node = tree.nodes[name]
    node.outbounds[1] ≡ nothing ?
        node.outbounds = (branch, node.outbounds[2]) :
        (node.outbounds[2] ≡ nothing ?
         node.outbounds = (node.outbounds[1], branch) :
         error("BinaryNode already has two outbound connections"))
end

function _addoutbound!(tree::PolytomousTree{<: Rooted}, name::String,
                       branch::Branch)
    node = tree.nodes[name]
    push!(node.outbounds, branch)
end

import Phylo.API: _removeoutbound!
function _removeoutbound!(tree::BinaryTree{<: Rooted}, name::String,
                          branch::Branch)
    node = tree.nodes[name]
    node.outbounds[1] ≡ branch ?
        node.outbounds = (node.outbounds[2], nothing) :
        (node.outbounds[2] ≡ branch ?
         node.outbounds = (node.outbounds[1], nothing) :
         error("BinaryNode does not have outbound connection to branch " *
               "$branch"))
end

_removeoutbound!(tree::PolytomousTree{<: Rooted}, name::String,
                 branch::Branch) =
                     filter!(x -> x ≢ branch, tree.nodes[name].outbounds)
