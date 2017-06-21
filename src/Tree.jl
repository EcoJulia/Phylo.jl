using DataStructures
importall Phylo.API

"""
    BinaryTree

Binary phylogenetic tree object with known leaves and per node data
"""
type BinaryTree{LI <: AbstractInfo, ND} <: AbstractTree{String, Int}
    nodes::OrderedDict{String, BinaryNode{Int}}
    branches::Dict{Int, Branch{String}}
    leafinfos::OrderedDict{String, LI}
    noderecords::OrderedDict{String, ND}
    rootheight::Nullable{Float64}
end

function BinaryTree{LI, ND}(lt::BinaryTree{LI, ND}; copyinfo=true, empty=true)
    validate(lt) || error("Tree to copy is not valid")
    leafnames = getleafnames(lt)
    # Leaf records may be conserved across trees, as could be invariant?
    leafinfos = copyinfo ? deepcopy(lt.leafinfos) : lt.leafinfos
    if empty # Empty out everything else
        nodes = OrderedDict(map(leaf -> leaf => BinaryNode{Int}(), leafnames))
        branches = OrderedDict{Int, Branch{String}}()
        noderecords = OrderedDict(map(leaf -> leaf => ND(), leafnames))
    else # Make copies of everything
        nodes = deepcopy(nodes)
        noderecords = deepcopy(lt.noderecords)
        branches = deepcopy(getbranches(lt))
    end
    return BinaryTree{LI, ND}(nodes, branches, leafinfos, noderecords,
                            lt.rootheight)
end

function (::Type{BinaryTree{LI, ND}}){LI <: AbstractInfo,
                                      ND}(leaves::Vector{String},
                                          treetype::Type{BinaryTree{LI, ND}} =
                                          BinaryTree{LI, ND};
                                          rootheight::Nullable{Float64} =
                                          Nullable{Float64}())
    nodes = OrderedDict(map(leaf -> leaf => BinaryNode{Int}(), leaves))
    leafinfos = OrderedDict(map(leaf -> leaf => LI(), leaves))
    noderecords = OrderedDict(map(leaf -> leaf => ND(), leaves))
    return BinaryTree{LI, ND}(nodes, OrderedDict{Int, Branch{String}}(),
                              leafinfos, noderecords, rootheight)
end

function (::Type{BinaryTree{LI, ND}}){LI <: AbstractInfo,
                                      ND}(numleaves::Int,
                                          treetype::Type{BinaryTree{LI, ND}} =
                                          BinaryTree{LI, ND};
                                          rootheight::Nullable{Float64} =
                                          Nullable{Float64}())
    leaves = map(num -> "Leaf $num", 1:numleaves)
    nodes = OrderedDict(map(leaf -> leaf => BinaryNode{Int}(), leaves))
    leafinfos = OrderedDict(map(leaf -> leaf => LI(), leaves))
    noderecords = OrderedDict(map(leaf -> leaf => ND(), leaves))
    return BinaryTree{LI, ND}(nodes, OrderedDict{Int, Branch{String}}(),
                              leafinfos, noderecords, rootheight)
end

_nodetype(::BinaryTree) = BinaryNode{Int}

_branchtype(::BinaryTree) = Branch{String}

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

function _getnoderecord(nt::BinaryTree, nodename)
    return nt.noderecords[nodename]
end

function _setnoderecord!(nt::BinaryTree, nodename, value)
    nt.noderecords[nodename] = value
end

function _addnode!{LI, NR}(tree::BinaryTree{LI, NR}, nodename)
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
    if Set(NodeNameIterator(tree, isleaf)) != Set(getleafnames(tree))
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
            if !(getheight(tree, leaf) - rootheight ≈
                 heighttoroot(tree, leaf))
                warn("Leaf height ($(getheight(tree, leaf))) for $leaf does not match branches")
                return false
            end
        end
    end
    return true
end

function _hasrootheight(tree::BinaryTree)
    return !isnull(tree.rootheight)
end

function _getrootheight(tree::BinaryTree)
    return get(tree.rootheight)
end

function _setrootheight!(tree::BinaryTree, height::Float64)
    tree.rootheight = height
    return height
end

function _clearrootheight!(tree::BinaryTree)
    tree.rootheight = Nullable{Float64}()
end

"""
    NamedTree

Binary phylogenetic tree object with known leaves
"""
const NamedTree = BinaryTree{LeafInfo, Void}






_getnodenames(tree::AbstractTree) = collect(keys(_getnodes(tree)))
_getbranchnames(tree::AbstractTree) = collect(keys(_getbranches(tree)))
#  - _hasnode()
_hasnode(tree::AbstractTree, label) = haskey(_getnodes(tree), label)
#  - _hasbranch()
_hasbranch(tree::AbstractTree, label) = haskey(_getbranches(tree), label)
#  - _addbranch!()
function _addbranch!(tree::AbstractTree, source, target, length::Float64, label)
    # Add the new branch
    _setbranch!(tree, label, Branch(source, target, length))
    
    # Update the associated source and target nodes
    _addoutbound!(getnode(tree, source), label)
    _setinbound!(getnode(tree, target), label)
    
    # Return updated tree
    return label
end
#  - _deletebranch!()
function _deletebranch!(tree::AbstractTree, label)
    # Find the branch
    branch = _getbranch(tree, label)
    # Remove branch reference from its source node
    _deleteoutbound!(_getnode(tree, _getsource(branch)), label)
    # Remove branch reference from its target node
    _deleteinbound!(_getnode(tree, _gettarget(branch)), label)
    # Remove branch itself
    delete!(_getbranches(tree), label)
    # Return the branch label
    return label
end

#  - _deletenode!()
function _deletenode!(tree::AbstractTree, label)
    node = _getnode(tree, label)
    if _hasinbound(node)
        _deletebranch!(tree, _getinbound(node))
    end
    for b in _getoutbounds(node)
        _deletebranch!(tree, b)
    end
    delete!(_getnodes(tree), label)
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
function clearrootheight!(tree::BinaryTree)
    _clearrootheight!(tree)
end

function _getnode(tree::BinaryTree, label)
    return _getnodes(tree)[label]
end

function _getbranch(tree::BinaryTree, label)
    return _getbranches(tree)[label]
end

function _setnode!(tree::BinaryTree, label, node)
    return _getnodes(tree)[label] = node
end

function _setbranch!(tree::BinaryTree, label, branch)
    _hasbranch(tree, label) &&
        error("Branch $label already exists")
    return _getbranches(tree)[label] = branch
end

# Interfaces to an info
# ---------------------

function _hasheight(tree::BinaryTree, label)
    return _hasheight(getleafinfo(tree, label))
end

function _getheight(tree::BinaryTree, label)
    return _getheight(getleafinfo(tree, label))
end

function _setheight!(tree::BinaryTree, label, height::Float64)
    ai = getleafinfo(tree, label)
    _setheight!(ai, height)
    setleafinfo!(tree, label, ai)
    return height
end
