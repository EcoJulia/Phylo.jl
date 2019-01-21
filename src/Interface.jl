using Phylo.API
using Compat: mapreduce
import LightGraphs: src, dst, indegree, outdegree, degree

# AbstractTree/Node/Branch type methods
"""
    treetype(::Type{AbstractTree})

Returns tree number (OneTree, ManyTrees) from a tree type.
"""
treetype(::Type{T}) where {TT <: TreeType, RT, NL, N, B,
    T <: AbstractTree{TT, RT, NL, N, B}} = TT

"""
    roottype(::Type{AbstractTree})
    roottype(::Type{AbstractNode})
    roottype(::Type{AbstractBranch})

Returns root type from a tree type.
"""
roottype(::Type{T}) where {TT, RT <: Rootedness, NL, N, B,
                           T <: AbstractTree{TT, RT, NL, N, B}} = RT
roottype(::Type{N}) where {RT <: Rootedness, NL,
                           N <: AbstractNode{RT, NL}} = RT
roottype(::Type{B}) where {RT <: Rootedness, NL,
                           B <: AbstractBranch{RT, NL}} = RT

"""
    nodenametype(::Type{AbstractTree})
    nodenametype(::Type{AbstractNode})
    nodenametype(::Type{AbstractBranch})

Returns type of node names from a tree type.
"""
nodenametype(::Type{T}) where {TT, RT, NL, N, B,
    T <: AbstractTree{TT, RT, NL, N, B}} = NL
nodenametype(::Type{N}) where {RT <: Rootedness, NL,
                               N <: AbstractNode{RT, NL}} = NL
nodenametype(::Type{B}) where {RT <: Rootedness, NL,
                               B <: AbstractBranch{RT, NL}} = NL

"""
    nodetype(::Type{AbstractTree})

Returns type of nodes from a tree type.
"""
nodetype(::Type{T}) where {TT, RT, NL, N <: AbstractNode, B,
    T <: AbstractTree{TT, RT, NL, N, B}} = N

"""
    branchnametype(::AbstractTree)

Returns type of branch names from a branch type.
"""
branchnametype(::Type{<: AbstractTree}) = Int
branchnametype(::Type{<: AbstractNode}) = Int
branchnametype(::Type{<: AbstractBranch}) = Int

"""
    branchtype(::Type{AbstractTree})

Returns type of branches from a tree type.
"""
branchtype(::Type{T}) where {TT, RT, NL, N, B <: AbstractBranch,
    T <: AbstractTree{TT, RT, NL, N, B}} = B

"""
    treenametype(::Type{AbstractTree})

Returns the name type for a tree type.
"""
treenametype(treetype::Type{<: AbstractTree}) = _treenametype(treetype)

"""
    leafinfotype(::Type{<: AbstractTree})

retrieve the leaf info type of a tree.
"""
leafinfotype(::Type{<: AbstractTree}) = _leafinfotype(T)

"""
    noderecordtype(::Type{<: AbstractTree})

retrieve the node info type of a tree.
"""
noderecordtype(::Type{<: AbstractTree}) = _noderecordtype(T)

"""
    branchrecordtype(::Type{<: AbstractTree})

retrieve the branch info type of a tree.
"""
branchrecordtype(::Type{<: AbstractTree}) = _branchrecordtype(T)

# AbstractTree methods
"""
    ntrees(tree::AbstractTree)

Returns the number of trees in a tree object, 1 for a OneTree tree type, and
the count of trees for a ManyTrees type.
"""
ntrees(tree::AbstractTree) = _ntrees(tree)

"""
    gettrees(tree::AbstractTree)

Returns an iterable of trees.
"""
gettrees(tree::AbstractTree) = _gettrees(tree)

"""
    gettree(tree::AbstractTree, label)

Returns a single OneTree object tree corresponding to the label given.
"""
function gettree end
gettree(tree::AbstractTree{OneTree}) = tree
function gettree(tree::AbstractTree{ManyTrees})
    @assert _ntrees(tree) == 1 "Must be only one tree (found $(_ntrees(tree)))"
    return first(_gettrees(tree))
end
gettree(tree::AbstractTree, label) = _gettree(tree, label)

"""
    gettreenames(tree::AbstractTree)

Returns the names of the trees.
"""
function gettreenames end
gettreenames(tree::AbstractTree{OneTree}) = [_gettreename(tree)]
gettreenames(tree::AbstractTree{ManyTrees}) = _gettreenames(tree)

"""
    gettreename(tree::AbstractTree)

Returns the name of the single tree.
"""
gettreename(tree::AbstractTree{OneTree}) = _gettreename(tree)
function gettreename(tree::AbstractTree{ManyTrees})
    @assert _ntrees(tree) == 1 "Must be only one tree (found $(_ntrees(tree)))"
    return first(_gettreenames(tree))
end

"""
    nleaves(::AbstractTree)

Returns the number of leaves (tips) in a tree.
"""
nleaves(tree::AbstractTree) = _nleaves(tree)

"""
    nroots(::AbstractTree)

Returns the number of roots in a tree. For OneTree types, Unrooted trees will
return 0, OneRoot trees should return 1, and manyroots tree (ones with multiple
subtrees) will return the number of subtrees. Manytrees types will return a
Dict of counts of the number of roots for each tree in the set.
"""
function nroots end
nroots(::AbstractTree{OneTree, Unrooted}) = 0
nroots(tree::AbstractTree{OneTree, <: Rooted}) = _nroots(tree)
nroots(tree::AbstractTree{ManyTrees}) =
    Dict(name => nroots(gettree(tree, name)) for name in _gettreenames(tree))

"""
    getroots(::AbstractTree)
    getroots(::AbstractTree, id)

Returns a vector containing the root(s) of a single (OneTree) tree
or a set of (ManyTrees) trees.
"""
function getroots end
getroots(tree::AbstractTree{OneTree}) = _getroots(tree)
getroots(tree::AbstractTree{OneTree, Unrooted, NL, N, B}) where {NL, N, B} = N[]
getroots(tree::AbstractTree{ManyTrees}) =
    Dict(name => getroots(gettree(tree, name)) for name in _gettreenames(tree))

"""
    getroot(::AbstractTree)

Returns the root of a single tree (must be only one tree for a ManyTrees tree).
"""
function getroot end
getroot(tree::AbstractTree{OneTree, <: Rooted}) = _getroot(tree)
getroot(tree::AbstractTree{ManyTrees, <: Rooted}) = getroot(gettree(tree))

"""
    getnodes(::AbstractTree)

Returns the vector of nodes of a single tree, or a Dict of vectors of nodes
for multiple trees.
"""
function getnodes end
getnodes(tree::AbstractTree{OneTree}) = _getnodes(tree)
getnodes(tree::AbstractTree{ManyTrees}) =
    Dict(name => getnodes(gettree(tree, name)) for name in _gettreenames(tree))

"""
    nnodes(::AbstractTree)

Returns the number of nodes of a single tree, or a Dict of numbers of nodes
for multiple trees.
"""
function nnodes end
nnodes(tree::AbstractTree{OneTree}) = _nnodes(tree)
nnodes(tree::AbstractTree{ManyTrees}) =
    Dict(name => nnodes(gettree(tree, name)) for name in _gettreenames(tree))

"""
    nbranches(::AbstractTree)

Returns the number of branches of a single tree, or a Dict of numbers of
branches for multiple trees.
"""
function nbranches end
nbranches(tree::AbstractTree{OneTree}) = _nbranches(tree)
nbranches(tree::AbstractTree{ManyTrees}) =
    Dict(name => nbranches(gettree(tree, name)) for name in _gettreenames(tree))

"""
    getnodenames(tree::AbstractTree)

Return a vector of node names of a single tree (identified by id for a
ManyTrees tree), or a Dict of vectors of node names for multiple trees.
"""
function getnodenames end
getnodenames(tree::AbstractTree{OneTree}) = _getnodenames(tree)
getnodenames(tree::AbstractTree{ManyTrees}) =
    Dict(name => getnodenames(gettree(tree, name))
         for name in _gettreenames(tree))

"""
    getbranches(::AbstractTree)

Returns the vector of branches of a single tree, or a Dict of vectors of
branches for multiple trees.
"""
function getbranches end
getbranches(tree::AbstractTree{OneTree}) = _getbranches(tree)
getbranches(tree::AbstractTree{ManyTrees}) =
    Dict(name => getbranches(gettree(tree, name))
         for name in _gettreenames(tree))

"""
    getbranchnames(tree::AbstractTree)

Return a vector of branch names of a single tree, or a Dict of vectors of
branch names for multiple trees.
"""
function getbranchnames end
getbranchnames(tree::AbstractTree{OneTree}) = _getbranchnames(tree)
getbranchnames(tree::AbstractTree{ManyTrees}) =
    Dict(name => getbranchnames(gettree(tree, name))
         for name in _gettreenames(tree))

"""
    createbranch!(tree::AbstractTree, source, destination[, length::Float64];
                  data)

Add a branch from `source` to `destination` on `tree` with optional length and
data. `source` and `destination` can be either nodes or nodenames.
"""
function createbranch!(tree::AbstractTree{OneTree, <: Rooted},
                       source, destination, length::Float64 = NaN;
                       name::Int = _newbranchlabel(tree), data = missing)
    sn = _getnode(tree, source)
    dn = _getnode(tree, destination)
    _hasnode(tree, sn) ||
        error("Tree does not have an available source node called " *
              _getnodename(tree, source))
    _hasoutboundspace(tree, sn) ||
        error(_getnodename(tree, source) *
              " already has maximum number of outbound connections ($(_outdegree(tree, sn)))")
    _hasnode(tree, dn) ||
        error("Tree does not have a destination node called " *
              _getnodename(tree, destination))
    _hasinboundspace(tree, dn) ||
        error("Tree does not have an available destination node called " *
              _getnodename(tree, destination))
    sn !== dn || error("Branch must connect different nodes")

    return ismissing(data) ? _createbranch!(tree, sn, dn, length, name) :
        _createbranch!(tree, sn, dn, length, name, data)
end
function createbranch!(tree::AbstractTree{OneTree, Unrooted},
                       source, destination, length::Float64 = NaN;
                       name::Int = _newbranchlabel(tree), data = missing)
    sn = _getnode(tree, source)
    dn = _getnode(tree, destination)
    _hasnode(tree, sn) ||
        error("Tree does not have an available source node called " *
              _getnodename(tree, source))
    _hasspace(tree, sn) ||
        error(_getnodename(tree, source) *
              " already has maximum number of connections " *
              "$(_degree(tree, sn)))")
    _hasnode(tree, dn) ||
        error("Tree does not have a destination node called " *
              _getnodename(tree, destination))
    _hasspace(tree, dn) ||
        error(_getnodename(tree, destination) *
              " already has maximum number of connections " *
              "$(_degree(tree, dn)))")
    sn !== dn || error("Branch must connect different nodes")

    return ismissing(data) ? _createbranch!(tree, sn, dn, length, name) :
        _createbranch!(tree, sn, dn, length, name, data)
end

"""
    deletebranch!(tree::AbstractTree, branch)
    deletebranch!(tree::AbstractTree, src, dst)

Delete the branch `branch` from `tree`, or branch connecting `src` node to
`dst` node.
"""
function deletebranch! end
function deletebranch!(tree::AbstractTree{OneTree}, branch)
    _hasbranch(tree, branch) ||
        error("Tree does not have a branch called " *
              "$(_getbranchname(tree, branch))")
    return _deletebranch!(tree, branch)
end
function deletebranch!(tree::AbstractTree{OneTree, RT, NL, N, B},
                       source::Union{N, NL}, dest::Union{N, NL}) where
    {RT, NL, N, B}
    _hasbranch(tree, source, dest) ||
        error("Tree does not have a branch between " *
              "$(_getnodename(tree, source)) and $(_getnodename(tree, dest))")
    return _deletebranch!(tree, _getbranch(tree, source, dest))
end

"""
    createnode!(tree::AbstractTree[, nodename]; data)

Create a node on a tree with optional node info.
"""
function createnode!(tree::AbstractTree{OneTree, RT, NL, N, B},
                     nodename::NL = _newnodelabel(tree);
                     data = missing) where {RT, NL, N, B}
    !_hasnode(tree, nodename) ||
        error("Node $nodename already present in tree")
    return ismissing(data) ? _createnode!(tree, nodename) :
        _createnode!(tree, nodename, data)
end

"""
    createnodes!(tree::AbstractTree, count::Integer)
    createnodes!(tree::AbstractTree, nodenames)
    createnodes!(tree::AbstractTree, nodedict)

Add a number of nodes, a vector with given names, or a Dict with node names
and associated node info to a tree.
"""
function createnodes! end
createnodes!(tree::AbstractTree{OneTree}, count::Int) =
    [createnode!(tree) for count in Base.OneTo(count)]
function createnodes!(tree::T, nodenames::V) where
    {RT, NL, N, B,
     T <: AbstractTree{OneTree, RT, NL, N, B},
     V <: AbstractVector{NL}}
    here = [name for name in nodenames if _hasnode(tree, name)]
    isempty(here) ||
        error("Nodes $here already present in tree")
    return [_createnode!(tree, name) for name in nodenames]
end
function createnodes!(tree::T, nodedata::NDS) where
    {RT, NL, N, B, T <: AbstractTree{OneTree, RT, NL, N, B},
     ND, NDS <: Dict{NL, ND}}
    here = [name for name in keys(nodedata) if _hasnode(tree, name)]
    isempty(here) || error("Nodes $here already present in tree")
    noderecordtype(T) == ND || error("Node info types does not match")
    return [_createnode!(tree, info.first, info.second) for info in nodedata]
end

"""
    deletenode!(tree::AbstractTree, node)

Delete a node (or a name) from a tree
"""
function deletenode!(tree::AbstractTree{OneTree}, node)
    _hasnode(tree, node) || error("Tree does not contain node $node")
    return _deletenode!(tree, _getnode(tree, node))
end

"""
    hasnode(tree::AbstractTree, node)

Returns whether a tree has a given node (or node name) or not.
"""
hasnode(tree::AbstractTree{OneTree}, node) = _hasnode(tree, node)

"""
    getnode(tree::AbstractTree, nodename)

Returns a node from a tree.
"""
function getnode(tree::AbstractTree{OneTree}, node)
    _hasnode(tree, node) || error("Node $node does not exist")
    return _getnode(tree, node)
end

"""
    getnodename(::AbstractTree, node)

Returns the node name associated with a node from a tree. For some
node types, it will be able to extract the node name without reference to
the tree.
"""
function getnodename end
function getnodename(tree::AbstractTree{OneTree}, node)
    _hasnode(tree, node) || error("Node $node does not exist")
    return _getnodename(tree, node)
end

"""
    hasbranch(tree::AbstractTree, branch)
    hasbranch(tree::AbstractTree, source, dest)

Does `tree` have a branch `branch` or a branch from `source` to `dest`?
"""
function hasbranch end
hasbranch(tree::AbstractTree{OneTree}, branch) = _hasbranch(tree, branch)
function hasbranch(tree::AbstractTree{OneTree, RT, NL, N, B},
          source::Union{N, NL}, dest::Union{N, NL}) where {RT, NL, N, B}
    _hasnode(tree, source) || error("Node $source does not exist")
    _hasnode(tree, dest) || error("Node $dest does not exist")
    return _hasbranch(tree, source, dest)
end

"""
    getbranch(tree::AbstractTree, branch)
    getbranch(tree::AbstractTree, source, dest)

Returns a branch from a tree by name or by source and destination node.
"""
function getbranch(tree::AbstractTree{OneTree, RT, NL, N, B},
                   branch::Union{Int, B}) where {RT, NL, N, B}
    _hasbranch(tree, branch) || error("Branch $branch does not exist")
    return _getbranch(tree, branch)
end
function getbranch(tree::AbstractTree{OneTree, RT, NL, N, B},
                   source::Union{N, NL}, dest::Union{N, NL}) where {RT, NL, N, B}
    _hasbranch(tree, source, dest) || error("Branch $branch does not exist")
    return _getbranch(tree, source, dest)
end

"""
    getbranchname(::AbstractTree, branch)
    getbranchname(branch)

Returns the branch name associated with a branch from a tree. For some
branch types, it will be able to extract the branch name without reference to
the tree.
"""
function getbranchname end
getbranchname(tree::T, branch::Union{Int, B}) where
    {TT, RT, NL, N, B <: AbstractBranch{RT, NL},
     T <: AbstractTree{TT, RT, NL, N, B}} = _getbranchname(tree, branch)

"""
    hasrootheight(tree::AbstractTree)


"""
function hasrootheight end
hasrootheight(tree::AbstractTree{OneTree, <: Rooted}) = _hasrootheight(tree)
hasrootheight(tree::AbstractTree{OneTree, Unrooted}) = false
hasrootheight(tree::AbstractTree{ManyTrees}) =
    all(hasrootheight(tree) for tree in gettrees(tree))

"""
    getrootheight(tree::AbstractTree)


"""
function getrootheight(tree::AbstractTree{OneTree, <: Rooted})
    return _getrootheight(tree)
end

"""
    setrootheight!(tree::AbstractTree, height)


"""
function setrootheight!(tree::AbstractTree{OneTree, <: Rooted}, height)
    return _setrootheight!(tree, height)
end

# AbstractNode methods
"""
    isleaf(tree::AbstractTree, node)

Is the node (referenced by name or node object) a leaf of the tree?
"""
function isleaf end
isleaf(tree::AbstractTree{OneTree}, node) = _isleaf(tree, node)

"""
    isroot(tree::AbstractTree, node)

Is the node (referenced by name or node object) a root of the tree?
"""
function isroot end
isroot(tree::AbstractTree{OneTree}, node) = _isroot(tree, node)

"""
    isinternal(tree::AbstractTree, node)

Is the node (referenced by name or node object) internal to the tree (neither
root nor leaf)?
"""
function isinternal end
isinternal(tree::AbstractTree{OneTree}, node) = _isinternal(tree, node)

"""
    isunattached(tree::AbstractTree, node)

Is the node (referenced by name or node object) unattached (i.e. not connected
to other nodes)?
"""
function isunattached end
isunattached(tree::AbstractTree{OneTree}, node) = _isunattached(tree, node)

"""
    indegree(tree::AbstractTree, node)

Return in degree of node.
"""
function indegree end
indegree(tree::AbstractTree{OneTree}, node) = _indegree(tree, node)

"""
    outdegree(tree::AbstractTree, node)

Return out degree of node.
"""
function outdegree end
outdegree(tree::AbstractTree{OneTree}, node) = _outdegree(tree, node)

"""
    degree(tree::AbstractTree, node)

Return the degree of a node including all connections.
"""
function degree end
degree(tree::AbstractTree{OneTree}, node) = _degree(tree, node)

"""
    hasoutboundspace(tree::AbstractTree, node)

Does the node have space for an[other] outbound connection?
"""
function hasoutboundspace end
hasoutboundspace(tree::AbstractTree{OneTree}, node) =
    _hasoutboundspace(tree, node)

"""
    hasinbound(tree::AbstractTree, node)

Does the node have an inbound connection?
"""
function hasinbound end
hasinbound(tree::AbstractTree{OneTree}, node) = _hasinbound(tree, node)

"""
    hasinboundspace(tree::AbstractTree, node)

Does the node have space for an inbound connection?
"""
function hasinboundspace end
hasinboundspace(tree::AbstractTree{OneTree}, node) =
    _hasinboundspace(tree, node)

"""
    getinbound(tree::AbstractTree, node)

return the inbound branch to this node (returns name for node name, branch for
node).
"""
function getinbound end
getinbound(tree::AbstractTree{OneTree, <:Rooted}, node) =
    _getinbound(tree, node)

"""
    getparent(tree::AbstractTree, node)

Return [the name of] the parent node for this node [name]. Second method may
not be implemented for some node types.
"""
function getparent end
getparent(tree::AbstractTree{OneTree, <:Rooted}, node) =
    _getparent(tree, node)

"""
    getoutbounds(tree::AbstractTree, nodename)

Return the names of the outbound branches from this node.
"""
function getoutbounds end
getoutbounds(tree::AbstractTree{OneTree, <:Rooted}, node) =
    _getoutbounds(tree, node)

"""
    getchildren(tree::AbstractTree, node)

Return the [name(s) of] the child node(s) for this node [name].
"""
function getchildren end
getchildren(tree::AbstractTree{OneTree, <:Rooted}, node) =
    _getchildren(tree, node)

"""
    getconnections(tree::AbstractTree, nodee)

Returns all of the branches connected to a node.
"""
function getconnections end
getconnections(tree::AbstractTree{OneTree}, node) =
    _getconnections(tree, node)

"""
    getsiblings(tree::AbstractTree, node)

Returns all of the siblings of a node. Must be implemented for any unrooted
AbstractNode subtype, can be inferred from _getparent and _getchildren for
a rooted node.
"""
function getsiblings end
getsiblings(tree::AbstractTree{OneTree}, node) = _getsiblings(tree, node)

"""
    hasheight(tree::AbstractTree, node)

Does the node have a height defined?
"""
function hasheight end
hasheight(tree::AbstractTree, node) = false
function hasheight(tree::AbstractTree{OneTree, <:Rooted}, node)
    return _hasheight(tree, node) ||
        (_hasrootheight(tree) &&
         mapreduce(b -> haslength(tree, b), &, branchhistory(tree, node);
         init = _hasrootheight(tree)))
end

"""
    getheight(tree::AbstractTree, node)

Return the height of the node.
"""
function getheight(tree::AbstractTree{OneTree, <:Rooted}, node)
    return _hasheight(tree, node) ? _getheight(tree, node) :
        mapreduce(b -> getlength(tree, b), +, branchhistory(tree, node);
                  init = getrootheight(tree))
end

"""
    setheight!(tree::AbstractTree, nodename, height)

Set the height of the node.
"""
function setheight!(tree::AbstractTree{OneTree, <:Rooted}, nodename, height)
    return _setheight!(tree, nodename, height)
end

# AbstractBranch methods
"""
    src(tree::AbstractTree, branch)

Return the source node for this branch.
"""
function src end
src(tree::AbstractTree{OneTree, <:Rooted}, branch) = _src(tree, branch)

"""
    dst(tree::AbstractTree, branch)

Return the destination node for this branch.
"""
function dst end
dst(tree::AbstractTree{OneTree, <:Rooted}, branch) = _dst(tree, branch)

"""
    conns(tree::AbstractTree, branch)

Return the nodes connected to `branch`.
"""
function conns end
conns(tree::AbstractTree{OneTree}, branch) = _conns(tree, branch)

"""
    conn(tree::AbstractTree, branch, exclude)

Return the other node connected to `branch` that is not `exclude`.
"""
function conn end
conn(tree::AbstractTree{OneTree}, branch, exclude) =
    _conn(tree, branch, exclude)

"""
    getlength(tree::AbstractTree, branch)

Return the length of this branch.
"""
function getlength end
getlength(tree::AbstractTree{OneTree}, branch) = _getlength(tree, branch)

"""
    getleafnames(::AbstractTree)

Retrieve the leaf names from the tree.
"""
getleafnames(tree::AbstractTree) = _getleafnames(tree)

"""
    getleaves(::AbstractTree)

Retrieve the leaves from the tree.
"""
getleaves(tree::AbstractTree) = _getleaves(tree)

"""
    getleafinfo(::AbstractTree, label)

retrieve the leaf info for a leaf of the tree.
"""
getleafinfo(tree::AbstractTree) = _getleafinfo(tree)

"""
    setleafinfo!(::AbstractTree, table)

Set the leaf info for the leaves of the tree.
"""
setleafinfo!(tree::AbstractTree, table) = _setleafinfo!(tree, table)

"""
    getnoderecord(::AbstractTree, label)

retrieve the node info for a leaf of the tree.
"""
getnoderecord(tree::AbstractTree{OneTree}, label) = _getnoderecord(tree, label)

"""
    setnoderecord!(::AbstractTree, label, value)

Set the node info for a node of the tree.
"""
setnoderecord!(tree::AbstractTree{OneTree}, label, value) =
    _setnoderecord!(tree, label, value)

"""
    getbranchinfo(::AbstractTree, label)

retrieve the branch info for a leaf of the tree.
"""
getbranchinfo(tree::AbstractTree{OneTree}, label) = _getbranchinfo(tree, label)

"""
    setbranchinfo!(::AbstractTree, label, value)

Set the branch info for a branch of the tree.
"""
setbranchinfo!(tree::AbstractTree{OneTree}, label, value) =
    _setbranchinfo!(tree, label, value)

"""
    resetleaves!(::AbstractTree)

Reset the leaf records to the current leaves, deleting all leaf records.
"""
resetleaves!(tree::AbstractTree{OneTree}) = _resetleaves!(tree)

"""
    validate(tree::AbstractTree)

Validate the tree by making sure that it is connected up correctly.
"""
function validate(tree::T) where
    {TT, RT, NL, N, B, T <: AbstractTree{TT, RT, NL, N, B}}
    nodes = _getnodes(tree)
    nodenames = _getnodenames(tree)
    branches = _getbranches(tree)
    branchnames = _getbranchnames(tree)
    if !isempty(nodes) || !isempty(branches)
        # We need to validate the connections
        if Set(_getinbound(tree, node) for node in nodes
               if _hasinbound(tree, node)) != Set(branches)
            warn("Inbound branches must exactly match Branch labels")
            return false
        end

        if Set(mapreduce(node -> _getoutbounds(tree, node), append!,
                         nodes; init = B[])) != Set(branches)
            warn("Node outbound branches must exactly match Branch labels")
            return false
        end

        if !(Set(_src(tree, branch) for branch in branches) ⊆ Set(nodenames))
            warn("Branch sources must be node labels")
            return false
        end

        if !(Set(_dst(tree, branch) for branch in branches) ⊆ Set(nodenames))
            warn("Branch destinations must be node labels")
            return false
        end
    end

    return _validate(tree)
end

"""
    traversal(::AbstractTree, ::TraversalOrder)
    traversal(::AbstractTree, ::TraversalOrder, init)

Return an iterable object for a tree containing nodes in given order -
preorder, inorder, postorder or breadthfirst - optionally starting from init.
"""
function traversal end
traversal(tree::AbstractTree, order::TraversalOrder = preorder) =
    _traversal(tree, order)
traversal(tree::AbstractTree{TT, RT, NL, N},
          order::TraversalOrder, init::Union{N, NL}) where {TT, RT, NL, N} =
    _traversal(tree, order, [init])

"""
    getancestors(tree::AbstractTree, nodename)

Return the name of all of the nodes that are ancestral to this node.
"""
function getancestors(tree::AbstractTree{OneTree, <:Rooted}, nodename)
    return _treepast(tree, nodename)[2][2:end]
end


"""
    getdescendants(tree::AbstractTree, nodename)

Return the names of all of the nodes that descend from this node.
"""
function getdescendants(tree::AbstractTree{OneTree, <:Rooted}, nodename)
    return _treefuture(tree, nodename)[2][2:end]
end
