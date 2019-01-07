using DataStructures
using Compat
using Compat: @warn
using IterableTables: getiterator
using DataFrames

abstract type AbstractBranchTree{RT, N, B, LI, ND} <:
    AbstractTree{OneTree, RT, String, N, B}
end

import Phylo.API: _leafinfotype
_leafinfotype(::Type{<: AbstractBranchTree{RT, N, B, LI, ND}}) where
    {RT, N, B, LI, ND} = LI

import Phylo.API: _noderecordtype
_noderecordtype(::Type{<: AbstractBranchTree{RT, N, B, LI, ND}}) where
    {RT, N, B, LI, ND} = ND

import Phylo.API: _branchrecordtype
_branchrecordtype(::Type{<: AbstractBranchTree}) = Nothing

import Phylo.API: _getnodes
_getnodes(bt::AbstractBranchTree) = values(bt.nodes)

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
    bt.leafinfos = empty!(bt.leafinfos)
    return bt
end

import Phylo.API: _getnoderecord
function _getnoderecord(bt::AbstractBranchTree, name)
    return bt.noderecords[name]
end

import Phylo.API: _setnoderecord!
function _setnoderecord!(bt::AbstractBranchTree, name, value)
    bt.noderecords[name] = value
end

import Phylo.API: _hasnode
_hasnode(tree::AbstractBranchTree, name::String) = haskey(tree.nodes, name)

import Phylo.API: _getnode
_getnode(tree::AbstractBranchTree, name::String) = tree.nodes[name]

import Phylo.API: _createnode!
function _createnode!(tree::AbstractBranchTree{RT, N, B, LI, ND},
                      name::String, data::ND = ND()) where {RT, N, B, LI, ND}
    _hasnode(tree, name) && error("Node $name already present in tree")
    tree.nodes[name] = N(name)
    _setnoderecord!(tree, name, data)
    return name
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
    delete!(tree.noderecords, name)
    return name
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
    delete!(tree.noderecords, name)
    return name
end

import Phylo.API: _hasbranch
_hasbranch(tree::AbstractBranchTree, name::Int) = haskey(tree.branches, name)

import Phylo.API: _getbranch
_getbranch(tree::AbstractBranchTree, name::Int) = tree.branches[name]

import Phylo.API: _createbranch!
function _createbranch!(tree::AbstractBranchTree{RT}, source, destination,
                        length::Float64, name::Int,
                        data::Nothing = nothing) where RT
    # Add the new branch
    branch = Branch{RT}(name, _getnodename(tree, source),
                        _getnodename(tree, destination), length)
    tree.branches[name] = branch

    # Update the associated source and destination nodes
    _addoutbound!(tree, _getnode(tree, source), branch)
    _addinbound!(tree, _getnode(tree, destination), branch)

    # Return updated tree
    return name
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
    return name
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
    return branch
end

import Phylo.API: _validate
function _validate(tree::TREE) where {TREE <: AbstractBranchTree}
    if _leafinfotype(TREE) != Nothing && length(getiterator(tree.leafinfos)) > 0
        if Set(info[1] for info in getiterator(tree.leafinfos)) !=
            Set(_getleafnames(tree))
            @warn "LeafInfo names do not match actual leaves of tree"
            return false
        end
    end

    if Set(keys(tree.noderecords)) != Set(_getnodenames(tree))
        @warn "Node names do not match node records of tree"
        return false
    end
    return true
end

import Phylo.API: _hasrootheight
function _hasrootheight(tree::AbstractBranchTree)
    return !isnan(tree.rootheight)
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
    tree.rootheight = NaN
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
    noderecords::OrderedDict{String, ND}
    rootheight::Float64
end

function BinaryTree(lt::BinaryTree{RT, LI, ND};
                    copyinfo=false) where {RT, LI, ND}
    validate(lt) || error("Tree to copy is not valid")
    leafnames = getleafnames(lt)
    # Leaf records are conserved across trees, as could be invariant?
    leafinfos = copyinfo ? deepcopy(lt.leafinfos) : lt.leafinfos
    nodes = OrderedDict(leaf => BinaryNode{RT, String, Branch{RT, String}}(leaf)
                        for leaf in leafnames)
    branches = Dict{Int, Branch{RT, String}}()
    noderecords = OrderedDict(leaf => ND() for leaf in leafnames)
    return BinaryTree{RT, LI, ND}(nodes, branches, leafinfos, noderecords,
                                  lt.rootheight)
end

function BinaryTree{RT, LI, ND}(leafnames::Vector{String},
                                treetype::Type{BinaryTree{RT, LI, ND}} =
                                BinaryTree{RT, LI, ND};
                                rootheight::Float64 = NaN) where {RT, LI, ND}
    nodes = OrderedDict(leaf => BinaryNode{RT, String, Branch{RT, String}}(leaf)
                        for leaf in leafnames)
    leafinfos = LI()
    noderecords = OrderedDict(leaf => ND() for leaf in leafnames)
    return BinaryTree{RT, LI, ND}(nodes, Dict{Int, Branch{RT, String}}(),
                                  leafinfos, noderecords, rootheight)
end

function BinaryTree{RT, LI, ND}(numleaves::Int = 0,
                                treetype::Type{BinaryTree{RT, LI, ND}} =
                                BinaryTree{LI, ND};
                                rootheight::Float64 = NaN) where {RT, LI, ND}
    leafnames = ["Leaf $num" for num in Base.OneTo(numleaves)]
    return BinaryTree{RT, LI, ND}(leafnames, treetype; rootheight = rootheight)
end

function BinaryTree{LI, ND}(leafinfos::LI;
                            rootheight::Float64 = NaN) where {LI, ND}
    leafnames = unique(info[1] for info in getiterator(leafinfos))
    nodes = OrderedDict(leaf => BinaryNode{RT, String, Branch{RT, String}}(leaf)
                        for leaf in leafnames)
    branches = Dict{Int, Branch{RT, String}}()
    noderecords = OrderedDict(leaf => ND() for leaf in leafnames)
    return BinaryTree{RT, LI, ND}(nodes, branches, leafinfos,
                                  noderecords, rootheight)
end

BinaryTree{RT}(leafinfos::LI; rootheight::Float64 = NaN) where {RT, LI} =
    BinaryTree{RT, LI, Dict{String, Any}}(leafinfos; rootheight = rootheight)

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
    noderecords::OrderedDict{String, ND}
    rootheight::Float64
end

function PolytomousTree(lt::PolytomousTree{RT, LI, ND};
                        copyinfo=false) where {RT, LI, ND}
    validate(lt) || error("Tree to copy is not valid")
    leafnames = getleafnames(lt)
    # Leaf records may be conserved across trees, as could be invariant?
    leafinfos = copyinfo ? deepcopy(lt.leafinfos) : lt.leafinfos
    nodes = OrderedDict(leaf => Node{RT, String, Branch{RT, String}}(leaf)
                        for leaf in leafnames)
    branches = Dict{Int, Branch{String}}()
    noderecords = OrderedDict(leaf => ND() for leaf in leafnames)
    return PolytomousTree{RT, LI, ND}(nodes, branches, leafinfos,
                                      noderecords, lt.rootheight)
end

function PolytomousTree{RT, LI, ND}(leafnames::Vector{String},
                                    treetype::Type{PolytomousTree{RT, LI, ND}} =
                                        PolytomousTree{RT, LI, ND};
                                    rootheight::Float64 = NaN) where
    {RT, LI, ND}
    nodes = OrderedDict(leaf => Node{RT, String, Branch{RT, String}}(leaf)
                        for leaf in leafnames)
    leafinfos = LI()
    noderecords = OrderedDict(leaf => ND() for leaf in leafnames)
    return PolytomousTree{RT, LI, ND}(nodes,
                                      OrderedDict{Int,
                                                  Branch{RT, String}}(),
                                      leafinfos, noderecords, rootheight)
end

function PolytomousTree{RT, LI, ND}(numleaves::Int = 0,
                                    treetype::Type{PolytomousTree{RT, LI, ND}} =
                                        PolytomousTree{RT, LI, ND};
                                    rootheight::Float64 = NaN) where
    {RT, LI, ND}
    leafnames = ["Leaf $num"  for num in Base.OneTo(numleaves)]
    return PolytomousTree{RT, LI, ND}(leafnames, treetype;
                                      rootheight = rootheight)
end

function PolytomousTree{RT, LI, ND}(leafinfos::LI,
                                    treetype::Type{PolytomousTree{RT, LI, ND}} =
                                        PolytomousTree{RT, LI, ND};
                                    rootheight::Float64 = NaN) where
    {RT, LI, ND}
    leafnames = unique(info[1] for info in getiterator(leafinfos))
    nodes = OrderedDict(leaf => Node{RT, String, Branch{RT, String}}(leaf)
                        for leaf in leafnames)
    branches = Dict{Int, Branch{String}}()
    noderecords = OrderedDict(leaf => ND() for leaf in leafnames)
    return PolytomousTree{RT, LI, ND}(nodes, branches, leafinfos,
                                      noderecords, rootheight)
end

PolytomousTree{RT}(leafinfos::LI; rootheight::Float64 = NaN) where {RT, LI} =
    PolytomousTree{RT, LI, Dict{String, Any}}(leafinfos;
                                              rootheight = rootheight)

"""
    NamedPolytomousTree

Polytomous phylogenetic tree object with known leaves
"""
const NamedTree = NamedPolytomousTree =
    PolytomousTree{ManyRoots, DataFrame, Dict{String, Any}}
