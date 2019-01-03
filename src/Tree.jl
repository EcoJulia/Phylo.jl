using DataStructures
using Compat
using Compat: @warn
using IterableTables: getiterator
using DataFrames

import Phylo.API: _getnodes, _getbranches
import Phylo.API: _getnodenames, _getbranchnames, _getleafnames
import Phylo.API: _getleafinfo, _setleafinfo!, _leafinfotype
import Phylo.API: _resetleaves!, _getnoderecord, _setnoderecord!
import Phylo.API: _addnode!, _deletenode!, _addbranch!, _deletebranch!, _validate
import Phylo.API: _hasrootheight, _getrootheight, _setrootheight!, _clearrootheight!
import Phylo.API: _getnode, _getbranch, _nleaves

"""
    BinaryTree

Binary phylogenetic tree object with known leaves and per node data
"""
mutable struct BinaryTree{RT <: Rooted, LI, ND} <:
    AbstractTree{OneTree, RT, String,
                 BinaryNode{RT, String, Branch{RT, String}},
                 Branch{RT, String}}
    nodes::OrderedDict{String, BinaryNode{RT, String, Branch{RT, String}}}
    branches::Dict{Int, Branch{RT, String}}
    leafinfos::LI
    noderecords::OrderedDict{String, ND}
    rootheight::Float64
end

_leafinfotype(::Type{BinaryTree{RT, LI, ND}}) where {RT, LI, ND} = LI
_noderecordtype(::Type{BinaryTree{RT, LI, ND}}) where {RT, LI, ND} = ND

function BinaryTree(lt::BinaryTree{RT, LI, ND};
                    copyinfo=false) where {RT, LI, ND}
    validate(lt) || error("Tree to copy is not valid")
    leafnames = getleafnames(lt)
    # Leaf records may be conserved across trees, as could be invariant?
    leafinfos = copyinfo ? deepcopy(lt.leafinfos) : lt.leafinfos
    nodes = OrderedDict(leaf => BinaryNode{RT, String, Branch{RT, String}}() for leaf in leafnames)
    branches = Dict{Int, Branch{RT, String}}()
    noderecords = OrderedDict(leaf => ND() for leaf in leafnames)
    return BinaryTree{RT, LI, ND}(nodes, branches,
                                  leafinfos, noderecords, lt.rootheight)
end

function BinaryTree{RT, LI, ND}(leafnames::Vector{String},
                                treetype::Type{BinaryTree{RT, LI, ND}} =
                                BinaryTree{RT, LI, ND};
                                rootheight::Float64 = NaN) where {RT, LI, ND}
    nodes = OrderedDict(leaf => BinaryNode{RT, String, Branch{RT, String}}() for leaf in leafnames)
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

function BinaryTree{LI, ND}(leafinfos::LI; rootheight::Float64 = NaN) where {LI, ND}
    leafnames = unique(info[1] for info in getiterator(leafinfos))
    nodes = OrderedDict(leaf => BinaryNode{RT, String, Branch{RT, String}}() for leaf in leafnames)
    branches = Dict{Int, Branch{RT, String}}()
    noderecords = OrderedDict(leaf => ND() for leaf in leafnames)
    return BinaryTree{RT, LI, Dict{String, Any}}(nodes, branches, leafinfos,
                                                 noderecords, rootheight)
end

BinaryTree{RT}(leafinfos::LI; rootheight::Float64 = NaN) where
    {RT, LI} =
    BinaryTree{RT, LI, Dict{String, Any}}(leafinfos;
                                          rootheight = rootheight)

function _getnodes(bt::BinaryTree)
    return values(bt.nodes)
end

function _getnodenames(bt::BinaryTree)
    return keys(bt.nodes)
end

function _getbranches(bt::BinaryTree)
    return values(bt.branches)
end

function _getbranchnames(bt::BinaryTree)
    return keys(bt.branches)
end

function _getleafinfo(bt::BinaryTree)
    return bt.leafinfos
end

function _getleafinfo(bt::BinaryTree, leafname)
    return Iterators.filter(info -> info[1] == leafname,
                            getiterator(bt.leafinfos))
end

function _setleafinfo!(bt::BinaryTree, info)
    bt.leafinfos = info
end

function _resetleaves!(bt::BinaryTree)
    bt.leafinfos = empty!(bt.leafinfos)
    return bt
end

function _getnoderecord(bt::BinaryTree, nodename)
    return bt.noderecords[nodename]
end

function _setnoderecord!(bt::BinaryTree, nodename, value)
    bt.noderecords[nodename] = value
end

function _addnode!(tree::BinaryTree{RT, LI, NR}, nodename) where {RT, LI, NR}
    _hasnode(tree, nodename) && error("Node $nodename already present in tree")
    tree.nodes[nodename] = BinaryNode{RT, String, Branch{RT, String}}()
    tree.noderecords[nodename] = NR()
    return nodename
end

function _deletenode!(tree::BinaryTree, nodename)
    node = getnode(tree, nodename)
    if _hasinbound(tree, node)
        deletebranch!(tree, _getinbound(tree, node))
    end
    for b in _getoutbounds(tree, node)
        deletebranch!(tree, b)
    end
    delete!(tree.nodes, nodename)
    delete!(tree.noderecords, nodename)
    return nodename
end

function _validate(tree::BinaryTree)
    if !isempty(tree.leafinfos) && length(getiterator(tree.leafinfos)) > 0
        if Set(map(info -> info[1], getiterator(tree.leafinfos))) !=
            Set(_getleafnames(tree))
            @warn "LeafInfo names do not match actual leaves of tree"
            return false
        end
    end

    if Set(nodenamefilter(_hasoutboundspace, tree)) !=
        Set(nodenamefilter(_isleaf, tree))
        @warn "Nodes must have two or zero outbound connections."
        return false
    end

    if Set(keys(tree.noderecords)) != Set(keys(getnodes(tree)))
        @warn "Leaf records do not match node records of tree"
        return false
    end

    return true
end

function _hasrootheight(tree::BinaryTree)
    return !isnan(tree.rootheight)
end

function _getrootheight(tree::BinaryTree)
    return tree.rootheight
end

function _setrootheight!(tree::BinaryTree, height::Float64)
    tree.rootheight = height
    return height
end

function _clearrootheight!(tree::BinaryTree)
    tree.rootheight = NaN
end

"""
    NamedBinaryTree

Binary phylogenetic tree object with known leaves
"""
const NamedBinaryTree =
    BinaryTree{ManyRoots, DataFrame, Dict{String, Any}}

"""
    PolytomousTree

Phylogenetic tree object with polytomous branching, and known leaves and per node data
"""
mutable struct PolytomousTree{RT <: Rooted, LI, ND} <:
    AbstractTree{OneTree, RT, String,
                 Node{RT, String, Branch{RT, String}},
                 Branch{RT, String}}
    nodes::OrderedDict{String, Node{RT, String, Branch{RT, String}}}
    branches::Dict{Int, Branch{RT, String}}
    leafinfos::LI
    noderecords::OrderedDict{String, ND}
    rootheight::Float64
end

_leafinfotype(::Type{PolytomousTree{RT, LI, ND}}) where {RT, LI, ND} = LI
_noderecordtype(::Type{PolytomousTree{RT, LI, ND}}) where {RT, LI, ND} = ND

function PolytomousTree(lt::PolytomousTree{RT, LI, ND};
                        copyinfo=false) where {RT, LI, ND}
    validate(lt) || error("Tree to copy is not valid")
    leafnames = getleafnames(lt)
    # Leaf records may be conserved across trees, as could be invariant?
    leafinfos = copyinfo ? deepcopy(lt.leafinfos) : lt.leafinfos
    nodes = OrderedDict(map(leaf -> leaf => Node{RT, String, Branch{RT, String}}(), leafnames))
    branches = Dict{Int, Branch{String}}()
    noderecords = OrderedDict(map(leaf -> leaf => ND(), leafnames))
    return PolytomousTree{RT, LI, ND}(nodes, branches, leafinfos,
                                      noderecords, lt.rootheight)
end

function PolytomousTree{RT, LI, ND}(leafnames::Vector{String},
                                treetype::Type{PolytomousTree{RT, LI, ND}} =
                                PolytomousTree{LI, ND};
                                rootheight::Float64 = NaN) where
    {RT, LI, ND}
    nodes = OrderedDict(leaf => Node{RT, String, Branch{RT, String}}()
                        for leaf in leafnames)
    leafinfos = LI()
    noderecords = OrderedDict(leaf => ND() for leaf in leafnames)
    return PolytomousTree{RT, LI, ND}(nodes,
                                      OrderedDict{Int,
                                                  Branch{RT, String}}(),
                                      leafinfos, noderecords,
                                      rootheight)
end

function PolytomousTree{RT, LI, ND}(numleaves::Int = 0,
                                treetype::Type{PolytomousTree{RT, LI, ND}} =
                                PolytomousTree{LI, Dict{String, Any}};
                                rootheight::Float64 = NaN) where
    {RT, LI, ND}
    leafnames = ["Leaf $num"  for num in Base.OneTo(numleaves)]
    return PolytomousTree(leafnames, treetype, rootheight = rootheight)
end

function PolytomousTree{RT, LI, ND}(leafinfos::LI,
                                treetype::Type{PolytomousTree{RT, LI, ND}} =
                                PolytomousTree{RT, LI, Dict{String, Any}};
                                rootheight::Float64 = NaN) where
    {RT, LI, ND}
    leafnames = unique(info[1] for info in getiterator(leafinfos))
    nodes = OrderedDict(leaf => Node{RT, String, Branch{RT, String}}() for leaf in leafnames)
    branches = Dict{Int, Branch{String}}()
    noderecords = OrderedDict(leaf => ND() for leaf in leafnames)
    return PolytomousTree{RT, LI, ND}(nodes, branches, leafinfos,
                                      noderecords, rootheight)
end

PolytomousTree{RT}(leafinfos::LI; rootheight::Float64 = NaN) where {RT, LI} =
    PolytomousTree{RT, LI, Dict{String, Any}}(leafinfos;
                                              rootheight = rootheight)

function _getnodes(pt::PolytomousTree)
    return values(pt.nodes)
end

function _getnodenames(pt::PolytomousTree)
    return keys(pt.nodes)
end

function _getbranches(pt::PolytomousTree)
    return values(pt.branches)
end

function _getbranchnames(pt::PolytomousTree)
    return keys(pt.branches)
end

function _getleafinfo(pt::PolytomousTree)
    return pt.leafinfos
end

function _getleafinfo(pt::PolytomousTree, leaf)
    return Iterators.filter(info -> info[1] == leafname,
        getiterator(pt.leafinfos))
end

function _setleafinfo!(pt::PolytomousTree, info)
    pt.leafinfos = info
end

function _resetleaves!(pt::PolytomousTree)
    pt.leafinfos = empty!(pt.leafinfos)
    return pt
end

function _getnoderecord(pt::PolytomousTree, nodename)
    return pt.noderecords[nodename]
end

function _setnoderecord!(pt::PolytomousTree, nodename, value)
    pt.noderecords[nodename] = value
end

function _addnode!(pt::PolytomousTree{RT, LI, NR}, nodename) where {RT, LI, NR}
    _hasnode(pt, nodename) && error("Node $nodename already present in tree")
    pt.nodes[nodename] = PolytomousNode{RT, String, Branch{RT, String}}()
    pt.noderecords[nodename] = NR()
    return nodename
end

function _deletenode!(pt::PolytomousTree, nodename)
    node = getnode(pt, nodename)
    if _hasinbound(pt, node)
        deletebranch!(pt, _getinbound(pt, node))
    end
    for b in _getoutbounds(pt, node)
        deletebranch!(pt, b)
    end
    delete!(pt.nodes, nodename)
    delete!(pt.noderecords, nodename)
    return nodename
end

function _validate(tree::PolytomousTree)
    if length(getiterator(tree.leafinfos)) > 0
        if Set(info[1] for info in getiterator(tree.leafinfos)) !=
            Set(_getleafnames(tree))
            @warn "LeafInfo names do not match actual leaves of tree"
            return false
        end
    end
    if Set(keys(tree.noderecords)) != Set(keys(getnodes(tree)))
        @warn "Leaf records do not match node records of tree"
        return false
    end

    return true
end

function _hasrootheight(tree::PolytomousTree)
    return !isnan(tree.rootheight)
end

function _getrootheight(tree::PolytomousTree)
    return tree.rootheight
end

function _setrootheight!(tree::PolytomousTree, height::Float64)
    tree.rootheight = height
    return height
end

function _clearrootheight!(tree::PolytomousTree)
    tree.rootheight = NaN
end

"""
    NamedPolytomousTree

Polytomous phylogenetic tree object with known leaves
"""
const NamedTree = NamedPolytomousTree =
    PolytomousTree{ManyRoots, DataFrame, Dict{String, Any}}

_getnodenames(tree::AbstractTree) = collect(keys(_getnodes(tree)))
_getbranchnames(tree::AbstractTree) = collect(keys(_getbranches(tree)))
#  - _hasnode()
_hasnode(tree::AbstractTree, label) = haskey(_getnodes(tree), label)
#  - _hasbranch()
_hasbranch(tree::AbstractTree, label) = haskey(_getbranches(tree), label)
#  - _addbranch!()
function _addbranch!(tree::PolytomousTree, source, destination, length::Float64, label)
    # Add the new branch
    _setbranch!(tree, label, Branch(source, destination, length))

    # Update the associated source and destination nodes
    _addoutbound!(tree, getnode(tree, source), label)
    _addinbound!(tree, getnode(tree, destination), label)

    # Return updated tree
    return label
end
#  - _deletebranch!()
function _deletebranch!(tree::PolytomousTree, label)
    # Find the branch
    branch = _getbranch(tree, label)
    # Remove branch reference from its source node
    _deleteoutbound!(tree, _getnode(tree, _src(tree, branch)), label)
    # Remove branch reference from its destination node
    _deleteinbound!(tree, _getnode(tree, _dst(tree, branch)), label)
    # Remove branch itself
    delete!(tree.branches, label)
    # Return the branch label
    return label
end

"""
    clearrootheight(::AbstractTree)

Clears the tree's root height record.
"""
function clearrootheight!(tree::AbstractTree)
    _clearrootheight!(tree)
end

function _getnode(tree::PolytomousTree, label)
    return tree.nodes[label]
end

function _getbranch(tree::PolytomousTree, label)
    return tree.branches[label]
end

function _setnode!(tree::AbstractTree, label, node)
    return tree.nodes[label] = node
end

function _setbranch!(tree::AbstractTree, label, branch)
    _hasbranch(tree, label) &&
        error("Branch $label already exists")
    return tree.branches[label] = branch
end
