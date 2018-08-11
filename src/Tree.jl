using DataStructures
using Compat
using Compat: @warn
using IterableTables: getiterator
using DataFrames

import Phylo.API: _nodetype, _branchtype, _getnodes, _getbranches
import Phylo.API: _getnodenames, _getbranchnames, _getleafnames
import Phylo.API: _getleafinfo, _setleafinfo!, _leafinfotype
import Phylo.API: _resetleaves!, _getnoderecord, _setnoderecord!
import Phylo.API: _addnode!, _deletenode!, _addbranch!, _deletebranch!, _validate
import Phylo.API: _hasrootheight, _getrootheight, _setrootheight!, _clearrootheight!
import Phylo.API: _getnode, _getbranch, _setnode!, _setbranch!, _nleaves

_branchtype(::AbstractBranchTree{NL, BL}) where {NL, BL} = Branch{NL}

"""
    BinaryTree

Binary phylogenetic tree object with known leaves and per node data
"""
mutable struct BinaryTree{LI, ND} <: AbstractBranchTree{String, Int}
    nodes::OrderedDict{String, BinaryNode{Int}}
    branches::Dict{Int, Branch{String}}
    leafinfos::LI
    noderecords::OrderedDict{String, ND}
    rootheight::Float64
end

_leafinfotype(::BinaryTree{LI, ND}) where {LI, ND} = LI
_nleaves(tree::BinaryTree{LI, ND}) where {LI, ND} =
    length(nodefilter(_isleaf, tree))

function BinaryTree(lt::BinaryTree{LI, ND};
                    copyinfo=false, empty=true) where {LI, ND}
    validate(lt) || error("Tree to copy is not valid")
    leafnames = getleafnames(lt)
    # Leaf records may be conserved across trees, as could be invariant?
    leafinfos = copyinfo ? deepcopy(lt.leafinfos) : lt.leafinfos
    if empty # Empty out everything else
        nodes = OrderedDict(map(leaf -> leaf => BinaryNode{Int}(), leafnames))
        branches = Dict{Int, Branch{String}}()
        noderecords = OrderedDict(map(leaf -> leaf => ND(), leafnames))
    else # Make copies of everything
        nodes = deepcopy(getnodes(lt))
        noderecords = deepcopy(lt.noderecords)
        branches = deepcopy(getbranches(lt))
    end
    return BinaryTree{LI, ND}(nodes, branches,
                              leafinfos, noderecords, lt.rootheight)
end

function BinaryTree{LI, ND}(leaves::Vector{String},
                            treetype::Type{BinaryTree{LI, ND}} =
                            BinaryTree{LI, ND};
                            rootheight::Float64 = NaN) where {LI, ND}
    nodes = OrderedDict(map(leaf -> leaf => BinaryNode{Int}(), leaves))
    leafinfos = LI()
    noderecords = OrderedDict(map(leaf -> leaf => ND(), leaves))
    return BinaryTree{LI, ND}(nodes, OrderedDict{Int, Branch{String}}(),
                              leafinfos, noderecords, rootheight)
end

function BinaryTree{LI, ND}(numleaves::Int = 0,
                            treetype::Type{BinaryTree{LI, ND}} =
                            BinaryTree{LI, ND};
                            rootheight::Float64 = NaN) where {LI, ND}
    leaves = map(num -> "Leaf $num", 1:numleaves)
    nodes = OrderedDict(map(leaf -> leaf => BinaryNode{Int}(), leaves))
    leafinfos = LI()
    noderecords = OrderedDict(map(leaf -> leaf => ND(), leaves))
    return BinaryTree{LI, ND}(nodes, OrderedDict{Int, Branch{String}}(),
                              leafinfos, noderecords, rootheight)
end

function BinaryTree{LI, ND}(leafinfos::LI; rootheight::Float64 = NaN) where {LI, ND}
    leafnames = unique(collect(map(info -> info[1], getiterator(leafinfos))))
    nodes = OrderedDict(map(leaf -> leaf => BinaryNode{Int}(), leafnames))
    branches = Dict{Int, Branch{String}}()
    noderecords = OrderedDict(map(leaf -> leaf => ND(),
                                  leafnames))
    return BinaryTree{LI,
                      Dict{String, Any}}(nodes, branches, leafinfos,
                                         noderecords, rootheight)
end

BinaryTree(leafinfos::LI; rootheight::Float64 = NaN) where LI =
    BinaryTree{LI, Dict{String, Any}}(leafinfos; rootheight = rootheight)

_nodetype(::BinaryTree) = BinaryNode{Int}

function _getnodes(bt::BinaryTree)
    return bt.nodes
end

function _getbranches(bt::BinaryTree)
    return bt.branches
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

function _addnode!(tree::BinaryTree{LI, NR}, nodename) where {LI, NR}
    !_hasnode(tree, nodename) ||
        error("Node $nodename already present in tree")
    _setnode!(tree, nodename, BinaryNode{Int}())
    setnoderecord!(tree, nodename, NR())
    return nodename
end

function _deletenode!(tree::BinaryTree, nodename)
    node = getnode(tree, nodename)
    if _hasinbound(node)
        deletebranch!(tree, _getinbound(node))
    end
    for b in _getoutbounds(node)
        deletebranch!(tree, b)
    end
    delete!(_getnodes(tree), nodename)
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
const NamedTree = NamedBinaryTree = BinaryTree{DataFrame, Dict{String, Any}}



"""
    PolytomousTree

Phylogenetic tree object with polytomous branching, and known leaves and per node data
"""
mutable struct PolytomousTree{LI, ND} <: AbstractBranchTree{String, Int}
    nodes::OrderedDict{String, Node{Int}}
    branches::Dict{Int, Branch{String}}
    leafinfos::LI
    noderecords::OrderedDict{String, ND}
    rootheight::Float64
end

_leafinfotype(::PolytomousTree{LI, ND}) where {LI, ND} = LI
_nleaves(tree::PolytomousTree{LI, ND}) where {LI, ND} =
    length(nodefilter(_isleaf, tree))

function PolytomousTree(lt::PolytomousTree{LI, ND};
                        copyinfo=false, empty=true) where {LI, ND}
    validate(lt) || error("Tree to copy is not valid")
    leafnames = getleafnames(lt)
    # Leaf records may be conserved across trees, as could be invariant?
    leafinfos = copyinfo ? deepcopy(lt.leafinfos) : lt.leafinfos
    if empty # Empty out everything else
        nodes = OrderedDict(map(leaf -> leaf => Node{Int}(), leafnames))
        branches = Dict{Int, Branch{String}}()
        noderecords = OrderedDict(map(leaf -> leaf => ND(), leafnames))
    else # Make copies of everything
        nodes = deepcopy(getnodes(lt))
        noderecords = deepcopy(lt.noderecords)
        branches = deepcopy(getbranches(lt))
    end
    return PolytomousTree{LI, ND}(nodes, branches, leafinfos, noderecords,
                            lt.rootheight)
end

function PolytomousTree{LI, ND}(leafinfos::LI,
                                treetype::Type{PolytomousTree{LI, ND}} =
                                PolytomousTree{LI, Dict{String, Any}};
                                rootheight::Float64 = NaN) where {LI, ND}
    leafnames = unique(collect(map(info -> info[1], getiterator(leafinfos))))
    nodes = OrderedDict(map(leaf -> leaf => Node{Int}(), leafnames))
    branches = Dict{Int, Branch{String}}()
    noderecords = OrderedDict(map(leaf -> leaf => ND(), leafnames))
    return PolytomousTree{LI, ND}(nodes, branches,
                                  leafinfos, noderecords, rootheight)
end

function PolytomousTree{LI, ND}(leaves::Vector{String},
                                treetype::Type{PolytomousTree{LI, ND}} =
                                PolytomousTree{LI, ND};
                                rootheight::Float64 = NaN) where {LI, ND}
    nodes = OrderedDict(map(leaf -> leaf => Node{Int}(), leaves))
    leafinfos = LI()
    noderecords = OrderedDict(map(leaf -> leaf => ND(), leaves))
    return PolytomousTree{LI, ND}(nodes, OrderedDict{Int, Branch{String}}(),
                              leafinfos, noderecords, rootheight)
end

function PolytomousTree{LI, ND}(numleaves::Int = 0,
                                treetype::Type{PolytomousTree{LI, ND}} =
                                PolytomousTree{LI, Dict{String, Any}};
                                rootheight::Float64 = NaN) where {LI, ND}
    leaves = map(num -> "Leaf $num", 1:numleaves)
    nodes = OrderedDict(map(leaf -> leaf => Node{Int}(), leaves))
    leafinfos = LI()
    noderecords = OrderedDict(map(leaf -> leaf => ND(), leaves))
    return PolytomousTree{LI, ND}(nodes, OrderedDict{Int, Branch{String}}(),
                                  leafinfos, noderecords, rootheight)
end

PolytomousTree(leafinfos::LI; rootheight::Float64 = NaN) where LI =
    PolytomousTree{LI, Dict{String, Any}}(leafinfos; rootheight = rootheight)

_nodetype(::PolytomousTree) = Node{Int}

function _getnodes(pt::PolytomousTree)
    return pt.nodes
end

function _getbranches(pt::PolytomousTree)
    return pt.branches
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

function _addnode!(tree::PolytomousTree{LI, NR}, nodename) where {LI, NR}
    !_hasnode(tree, nodename) ||
        error("Node $nodename already present in tree")
    _setnode!(tree, nodename, Node{Int}())
    setnoderecord!(tree, nodename, NR())
    return nodename
end

function _deletenode!(tree::PolytomousTree, nodename)
    node = getnode(tree, nodename)
    if _hasinbound(node)
        deletebranch!(tree, _getinbound(node))
    end
    for b in _getoutbounds(node)
        deletebranch!(tree, b)
    end
    delete!(_getnodes(tree), nodename)
    delete!(tree.noderecords, nodename)
    return nodename
end

function _validate(tree::PolytomousTree)
    if length(getiterator(tree.leafinfos)) > 0
        if Set(map(info -> info[1], getiterator(tree.leafinfos))) !=
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
const NamedPolytomousTree = PolytomousTree{DataFrame, Dict{String, Any}}

_getnodenames(tree::AbstractTree) = collect(keys(_getnodes(tree)))
_getbranchnames(tree::AbstractTree) = collect(keys(_getbranches(tree)))
#  - _hasnode()
_hasnode(tree::AbstractTree, label) = haskey(_getnodes(tree), label)
#  - _hasbranch()
_hasbranch(tree::AbstractTree, label) = haskey(_getbranches(tree), label)
#  - _addbranch!()
function _addbranch!(tree::AbstractBranchTree, source, destination, length::Float64, label)
    # Add the new branch
    _setbranch!(tree, label, Branch(source, destination, length))

    # Update the associated source and destination nodes
    _addoutbound!(getnode(tree, source), label)
    _setinbound!(getnode(tree, destination), label)

    # Return updated tree
    return label
end
#  - _deletebranch!()
function _deletebranch!(tree::AbstractBranchTree, label)
    # Find the branch
    branch = _getbranch(tree, label)
    # Remove branch reference from its source node
    _deleteoutbound!(_getnode(tree, _src(branch)), label)
    # Remove branch reference from its destination node
    _deleteinbound!(_getnode(tree, _dst(branch)), label)
    # Remove branch itself
    delete!(_getbranches(tree), label)
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

function _getnode(tree::AbstractTree, label)
    return _getnodes(tree)[label]
end

function _getbranch(tree::AbstractTree, label)
    return _getbranches(tree)[label]
end

function _setnode!(tree::AbstractTree, label, node)
    return _getnodes(tree)[label] = node
end

function _setbranch!(tree::AbstractTree, label, branch)
    _hasbranch(tree, label) &&
        error("Branch $label already exists")
    return _getbranches(tree)[label] = branch
end
