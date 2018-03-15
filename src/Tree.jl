using DataStructures
using Compat

import Phylo.API: _nodetype, _branchtype, _getnodes, _getbranches
import Phylo.API: _getnodenames, _getbranchnames, _getleafnames
import Phylo.API: _getleafinfo, _setleafinfo!, _resetleaves!, _getnoderecord, _setnoderecord!
import Phylo.API: _addnode!, _deletenode!, _addbranch!, _deletebranch!, _validate
import Phylo.API: _hasrootheight, _getrootheight, _setrootheight!, _clearrootheight!
import Phylo.API: _hasheight, _getheight, _setheight!, _getnode, _getbranch, _setnode!, _setbranch!

_branchtype(::AbstractBranchTree{NL, BL}) where {NL, BL} = Branch{NL}

"""
    BinaryTree

Binary phylogenetic tree object with known leaves and per node data
"""
mutable struct BinaryTree{LI <: AbstractInfo, ND} <: AbstractBranchTree{String, Int}
    nodes::OrderedDict{String, BinaryNode{Int}}
    branches::Dict{Int, Branch{String}}
    leafinfos::OrderedDict{String, LI}
    noderecords::OrderedDict{String, ND}
    rootheight::Float64
end

function BinaryTree(lt::BinaryTree{LI, ND}; copyinfo=true, empty=true) where {LI, ND}
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
    return BinaryTree{LI, ND}(nodes, branches, leafinfos, noderecords,
                            lt.rootheight)
end

function BinaryTree{LI, ND}(leaves::Vector{String},
                            treetype::Type{BinaryTree{LI, ND}} =
                            BinaryTree{LI, ND};
                            rootheight::Float64 = NaN) where {LI <: AbstractInfo, ND}
    nodes = OrderedDict(map(leaf -> leaf => BinaryNode{Int}(), leaves))
    leafinfos = OrderedDict(map(leaf -> leaf => LI(), leaves))
    noderecords = OrderedDict(map(leaf -> leaf => ND(), leaves))
    return BinaryTree{LI, ND}(nodes, OrderedDict{Int, Branch{String}}(),
                              leafinfos, noderecords, rootheight)
end

function BinaryTree{LI, ND}(numleaves::Int = 0,
                            treetype::Type{BinaryTree{LI, ND}} =
                            BinaryTree{LI, ND};
                            rootheight::Float64 = NaN) where {LI <: AbstractInfo, ND}
    leaves = map(num -> "Leaf $num", 1:numleaves)
    nodes = OrderedDict(map(leaf -> leaf => BinaryNode{Int}(), leaves))
    leafinfos = OrderedDict(map(leaf -> leaf => LI(), leaves))
    noderecords = OrderedDict(map(leaf -> leaf => ND(), leaves))
    return BinaryTree{LI, ND}(nodes, OrderedDict{Int, Branch{String}}(),
                              leafinfos, noderecords, rootheight)
end

_nodetype(::BinaryTree) = BinaryNode{Int}

function _getnodes(nt::BinaryTree)
    return nt.nodes
end

function _getbranches(nt::BinaryTree)
    return nt.branches
end

function _getleafnames(nt::BinaryTree)
    return keys(nt.leafinfos)
end

function _getleafinfo(nt::BinaryTree, leaf)
    return nt.leafinfos[leaf]
end

function _setleafinfo!(nt::BinaryTree, leaf, value)
    nt.leafinfos[leaf] = value
end

function _resetleaves!(bt::BinaryTree)
    bt.leafinfos = OrderedDict(map(name -> name => LeafInfo(),
                                   nodenamefilter(isleaf, bt)))
    return bt
end

function _getnoderecord(nt::BinaryTree, nodename)
    return nt.noderecords[nodename]
end

function _setnoderecord!(nt::BinaryTree, nodename, value)
    nt.noderecords[nodename] = value
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
    if Set(nodenamefilter(isleaf, tree)) != Set(getleafnames(tree))
        warn("Leaf names do not match actual leaves of tree")
        return false
    end

    if Set(nodenamefilter(hasoutboundspace, tree)) !=
        Set(nodenamefilter(isleaf, tree))
        warn("Nodes must have two or zero outbound connections.")
        return false
    end
    
    if Set(keys(tree.noderecords)) != Set(keys(getnodes(tree)))
        warn("Leaf records do not match node records of tree")
        return false
    end
    
    rootheight = hasrootheight(tree) ? getrootheight(tree) : NaN
    for leaf in getleafnames(tree)
        if hasheight(tree, leaf)
            if isnan(rootheight)
                rootheight = getheight(tree, leaf) - heighttoroot(tree, leaf)
            end
            if !(getheight(tree, leaf) - rootheight ≈ heighttoroot(tree, leaf))
                warn("Leaf height ($(getheight(tree, leaf))) for $leaf does not match branches")
                return false
            end
        end
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
const NamedTree = NamedBinaryTree = BinaryTree{LeafInfo, Nothing}



"""
    PolytomousTree

Phylogenetic tree object with polytomous branching, and known leaves and per node data
"""
mutable struct PolytomousTree{LI <: AbstractInfo, ND} <: AbstractBranchTree{String, Int}
    nodes::OrderedDict{String, Node{Int}}
    branches::Dict{Int, Branch{String}}
    leafinfos::OrderedDict{String, LI}
    noderecords::OrderedDict{String, ND}
    rootheight::Float64
end

function PolytomousTree(lt::PolytomousTree{LI, ND}; copyinfo=true, empty=true) where {LI, ND}
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

function PolytomousTree{LI, ND}(leaves::Vector{String},
                            treetype::Type{PolytomousTree{LI, ND}} =
                            PolytomousTree{LI, ND};
                            rootheight::Float64 = NaN) where {LI <: AbstractInfo, ND}
    nodes = OrderedDict(map(leaf -> leaf => Node{Int}(), leaves))
    leafinfos = OrderedDict(map(leaf -> leaf => LI(), leaves))
    noderecords = OrderedDict(map(leaf -> leaf => ND(), leaves))
    return PolytomousTree{LI, ND}(nodes, OrderedDict{Int, Branch{String}}(),
                              leafinfos, noderecords, rootheight)
end

function PolytomousTree{LI, ND}(numleaves::Int = 0,
                            treetype::Type{PolytomousTree{LI, ND}} =
                            PolytomousTree{LI, ND};
                            rootheight::Float64 = NaN) where {LI <: AbstractInfo, ND}
    leaves = map(num -> "Leaf $num", 1:numleaves)
    nodes = OrderedDict(map(leaf -> leaf => Node{Int}(), leaves))
    leafinfos = OrderedDict(map(leaf -> leaf => LI(), leaves))
    noderecords = OrderedDict(map(leaf -> leaf => ND(), leaves))
    return PolytomousTree{LI, ND}(nodes, OrderedDict{Int, Branch{String}}(),
                              leafinfos, noderecords, rootheight)
end

_nodetype(::PolytomousTree) = Node{Int}

function _getnodes(nt::PolytomousTree)
    return nt.nodes
end

function _getbranches(nt::PolytomousTree)
    return nt.branches
end

function _getleafnames(nt::PolytomousTree)
    return keys(nt.leafinfos)
end

function _getleafinfo(nt::PolytomousTree, leaf)
    return nt.leafinfos[leaf]
end

function _setleafinfo!(nt::PolytomousTree, leaf, value)
    nt.leafinfos[leaf] = value
end

function _resetleaves!(bt::PolytomousTree)
    bt.leafinfos = OrderedDict(map(name -> name => LeafInfo(),
                                   nodenamefilter(isleaf, bt)))
    return bt
end

function _getnoderecord(nt::PolytomousTree, nodename)
    return nt.noderecords[nodename]
end

function _setnoderecord!(nt::PolytomousTree, nodename, value)
    nt.noderecords[nodename] = value
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
    if Set(nodenamefilter(isleaf, tree)) != Set(getleafnames(tree))
        warn("Leaf names do not match actual leaves of tree")
        return false
    end

    if Set(keys(tree.noderecords)) != Set(keys(getnodes(tree)))
        warn("Leaf records do not match node records of tree")
        return false
    end
    
    rootheight = hasrootheight(tree) ? getrootheight(tree) : NaN
    for leaf in getleafnames(tree)
        if hasheight(tree, leaf)
            if isnan(rootheight)
                rootheight = getheight(tree, leaf) - heighttoroot(tree, leaf)
            end
            if !(getheight(tree, leaf) - rootheight ≈ heighttoroot(tree, leaf))
                warn("Leaf height ($(getheight(tree, leaf))) for $leaf does not match branches")
                return false
            end
        end
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

Binary phylogenetic tree object with known leaves
"""
const NamedPolytomousTree = PolytomousTree{LeafInfo, Nothing}




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


#  - _getleafnames()
function _getleafnames(tree::AbstractTree)
    return keys(OrderedDict(map(leaf -> leaf=>nothing,
                         findleaves(tree) ∪ findunattacheds(tree))))
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

# Interfaces to an info
# ---------------------

function _hasheight(tree::AbstractBranchTree, label)
    return _hasheight(getleafinfo(tree, label))
end

function _getheight(tree::AbstractBranchTree, label)
    return _getheight(getleafinfo(tree, label))
end

function _setheight!(tree::AbstractBranchTree, label, height::Float64)
    ai = getleafinfo(tree, label)
    _setheight!(ai, height)
    setleafinfo!(tree, label, ai)
    return height
end
