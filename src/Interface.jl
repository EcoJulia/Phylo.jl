using Phylo.API
using Compat: mapreduce

# AbstractTree type methods
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
    nodenametype(::AbstractTree)

Returns type of node names from a tree type.
"""
nodenametype(::Type{T}) where {TT, RT, NL, N, B,
    T <: AbstractTree{TT, RT, NL, N, B}} = NL

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
branchnametype(::Type{T}) where {TT, RT, NL, N, B,
    T <: AbstractTree{TT, RT, NL, N, B}} = Int

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
treenametype(treetype::Type{T}) where {TT, RT, NL, N, B,
    T <: AbstractTree{TT, RT, NL, N, B}} = _treenametype(treetype)

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
    getonetree(tree::AbstractTree, label)

Returns a OneTree object tree corresponding to the label given.
"""
function getonetree end
getonetree(tree::AbstractTree) = _getonetree(tree)
getonetree(tree::AbstractTree, label) = _getonetree(tree, label)

"""
    gettreenames(tree::AbstractTree)

Returns an iterable of trees.
"""
gettreenames(tree::AbstractTree) = _gettreenames(tree)

"""
    gettreename(tree::AbstractTree)

Returns an iterable of trees.
"""
gettreename(tree::AbstractTree{OneTree, RT, NL, N, B}) where {RT, NL, N, B} =
    _gettreename(tree)

"""
    nleaves(::AbstractTree)

Returns the number of leaves (tips) in a tree.
"""
nleaves(tree::AbstractTree) = _nleaves(tree)

"""
    nroots(::AbstractTree)
    nroots(::AbstractTree, label)

Returns the number of roots in a tree. For OneTree types, Unrooted trees will
return 0, OneRoot trees will return 1, and manyroots tree (ones with multiple
subtrees) will return the number of subtrees. manytree types will return a
vector of counts of the number of roots for each tree in the set.
"""
function nroots end
nroots(tree::AbstractTree, label) = _nroots(_getonetree(tree, label))
nroots(::AbstractTree{OneTree, Unrooted, NL, N, B}) where {NL, N, B} = 0
nroots(tree::AbstractTree{OneTree, RT, NL, N, B}) where
    {RT <: Rooted, NL, N, B} = _nroots(tree)
nroots(tree::AbstractTree{ManyTrees, RT, NL, N, B}) where
    {RT, NL, N, B} = [nroots(t) for t in _gettrees(tree)]

"""
    getroots(::AbstractTree)
    getroots(::AbstractTree, id)

Returns a vector containing the root(s) of a single (OneTree) tree
or a set of (ManyTrees) trees.
"""
function getroots end
getroots(tree::AbstractTree{OneTree, RT, NL, N, B}) where
    {RT, NL, N, B} = _getroots(tree)
getroots(tree::AbstractTree{OneTree, Unrooted, NL, N, B}) where
    {NL, N, B} = N[]
getroots(tree::AbstractTree{ManyTrees, RT, NL, N, B}) where
    {RT, NL, N, B} = [getroots(t) for t in _gettrees(tree)]


"""
    getroot(::AbstractTree)
    getroot(::AbstractTree, id)

Returns the root of a single tree (identified by id for a ManyTrees tree),
or a vector of roots for multiple trees.
"""
function getroot end
getroot(tree::AbstractTree{OneTree, OneRoot, NL, N, B}) where
    {NL, N, B} = _getroot(tree)
getroot(tree::AbstractTree{ManyTrees, OneRoot, NL, N, B}, id) where
    {NL, N, B} = getroot(_getonetree(tree, id))
getroot(tree::AbstractTree{ManyTrees, OneRoot, NL, N, B}) where
    {NL, N, B} = [getroot(t) for t in _gettrees(tree)]

"""
    getnodes(::AbstractTree)
    getnodes(::AbstractTree, id)

Returns the vector of nodes of a single tree (identified by id for a
ManyTrees tree), or a vector of vectors of nodes for multiple trees.
"""
function getnodes end
getnodes(tree::AbstractTree{OneTree, RT, NL, N, B}) where
    {RT, NL, N, B} = _getnodes(tree)
getnodes(tree::AbstractTree{ManyTrees, RT, NL, N, B}, id) where
    {RT, NL, N, B} = _getnodes(_getonetree(tree, id))
getnodes(tree::AbstractTree{ManyTrees, RT, NL, N, B}) where
    {RT, NL, N, B} = [getnodes(t) for t in _gettrees(tree)]

"""
    nnodes(::AbstractTree)
    nnodes(::AbstractTree, id)

Returns the number of nodes of a single tree (identified by id for a
ManyTrees tree), or a vector of numbers of nodes for multiple trees.
"""
function nnodes end
nnodes(tree::AbstractTree{OneTree, RT, NL, N, B}) where
    {RT, NL, N, B} = _nnodes(tree)
nnodes(tree::AbstractTree{ManyTrees, RT, NL, N, B}, id) where
    {RT, NL, N, B} = _nnodes(_getonetree(tree, id))
nnodes(tree::AbstractTree{ManyTrees, RT, NL, N, B}) where
    {RT, NL, N, B} = [nnodes(t) for t in _gettrees(tree)]

"""
    getnodenames(tree::AbstractTree)
    getnodenames(tree::AbstractTree, id)

Return a vector of node names of a single tree (identified by id for a
ManyTrees tree), or a vector of vectors of node names for multiple trees.
"""
function getnodenames end
getnodenames(tree::AbstractTree{OneTree, RT, NL, N, B}) where
    {RT, NL, N, B} = _getnodenames(tree)
getnodenames(tree::AbstractTree{ManyTrees, RT, NL, N, B}, id) where
    {RT, NL, N, B} = _getnodenames(_getonetree(tree, id))
getnodenames(tree::AbstractTree{ManyTrees, RT, NL, N, B}) where
    {RT, NL, N, B} = [_getnodenames(t) for t in _gettrees(tree)]

"""
    getbranches(::AbstractTree)
    getbranches(::AbstractTree, id)

Returns the vector of branches of a single tree (identified by id for a
ManyTrees tree), or a vector of vectors of branches for multiple trees.
"""
function getbranches end
getbranches(tree::AbstractTree{OneTree, RT, NL, N, B}) where
    {RT, NL, N, B} = _getbranches(tree)
getbranches(tree::AbstractTree{ManyTrees, RT, NL, N, B}, id) where
    {RT, NL, N, B} = _getbranches(_getonetree(tree, id))
getbranches(tree::AbstractTree{ManyTrees, RT, NL, N, B}) where
    {RT, NL, N, B} = [getbranches(t) for t in _gettrees(tree)]

"""
    getbranchnames(tree::AbstractTree)
    getbranchnames(tree::AbstractTree, id)

Return a vector of branch names of a single tree (identified by id for a
ManyTrees tree), or a vector of vectors of branch names for multiple trees.
"""
function getbranchnames end
getbranchnames(tree::AbstractTree{OneTree, RT, NL, N, B}) where
    {RT, NL, N, B} = _getbranchnames(tree)
getbranchnames(tree::AbstractTree{ManyTrees, RT, NL, N, B}, id) where
    {RT, NL, N, B} = _getbranchnames(_getonetree(tree, id))
getbranchnames(tree::AbstractTree{ManyTrees, RT, NL, N, B}) where
    {RT, NL, N, B} = [_getbranchnames(t) for t in _gettrees(tree)]

"""
    createbranch!(tree::AbstractTree, source, destination[, length::Float64];
                  data)

Add a branch from `source` to `destination` on `tree` with optional length and
data. `source` and `destination` can be either nodes or nodenames.
"""
function createbranch!(tree::AbstractTree{OneTree, RT, NL, N, B},
                       source::Union{N, NL}, destination::Union{N, NL},
                       length::Float64 = NaN; data = Dict{String, Any}()) where
    {RT, NL, N, B}
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

    return _createbranch!(tree, sn, dn, length, data)
end

"""
    deletebranch!(tree::AbstractTree, branch)

Delete the branch `branch` from `tree`.
"""
function deletebranch!(tree::AbstractTree{OneTree, RT, NL, N, B},
                       branch::Union{B, Int}) where {RT, NL, N, B}
    _hasbranch(tree, branch) ||
        error("Tree does not have a branch called " *
              _getbranchname(tree, branch))
    return _deletebranch!(tree, branch)
end

"""
    createnode!(tree::AbstractTree[, nodename]; data)

Create a node on a tree with optional node info.
"""
function createnode!(tree::AbstractTree{OneTree, RT, NL, N, B},
                     nodename::NL = _newnodelabel(tree);
                     data = Dict{String, Any}()) where {RT, NL, N, B}
    !_hasnode(tree, nodename) ||
        error("Node $nodename already present in tree")
    return _createnode!(tree, nodename, data)
end

"""
    createnodes!(tree::AbstractTree, count::Integer)
    createnodes!(tree::AbstractTree, nodenames)
    createnodes!(tree::AbstractTree, nodedict)

Add a number of nodes, a vector with given names, or a Dict with node names
and associated node info to a tree.
"""
function createnodes! end
createnodes!(tree::AbstractTree{OneTree, RT, NL, N, B}, count::Int) where
    {RT, NL, N, B} = map(name -> createnode!(tree), 1:count)
function createnodes!(tree::T, nodenames::V) where
    {RT, NL, N, B,
     T <: AbstractTree{OneTree, RT, NL, N, B},
     V <: AbstractVector{NL}}
    names = filter(name -> _hasnode(tree, name), nodenames)
    isempty(names) ||
        error("Nodes $names already present in tree")
    return map(name -> _createnode!(tree, name), nodenames)
end
function createnodes!(tree::T, nodedata::NDS) where
    {RT, NL, N, B, T <: AbstractTree{OneTree, RT, NL, N, B},
     ND, NDS <: Dict{NL, ND}}
    names = filter(name -> _hasnode(tree, name), keys(nodedata))
    isempty(names) ||
        error("Nodes $names already present in tree")
    return map(node -> _createnode!(tree, node.first, data = node.second),
               nodedata)
end

"""
    deletenode!(tree::AbstractTree, node)

Delete a node (or a name) from a tree
"""
function deletenode!(tree::AbstractTree{OneTree, RT, NL, N, B},
                     node::Union{NL, N}) where {RT, NL, N, B}
    return _deletenode!(tree, node)
end

"""
    hasnode(tree::AbstractTree, node)

Returns whether a tree has a given node (or node name) or not.
"""
function hasnode(tree::AbstractTree{OneTree, RT, NL, N, B},
                 node::Union{NL, N}) where {RT, NL, N, B}
    return _hasnode(tree, node)
end

"""
    getnode(tree::AbstractTree, nodename)

Returns a node from a tree.
"""
function getnode(tree::AbstractTree{OneTree, RT, NL, N, B},
                 node::Union{NL, N}) where {RT, NL, N, B}
    _hasnode(tree, node) ||
        error("Node $(_getnodename(tree, node)) does not exist")
    return _getnode(tree, node)
end

"""
    getnodename(::AbstractTree, id)
    getnodename(id)

Returns the node name associated with a node from a tree. For some id and
node types, it will be able to extract the node name without reference to
the tree.
"""
function getnodename end
getnodename(tree::T, node::Union{NL, N}) where
    {TT, RT, NL, N <: AbstractNode{RT, NL}, B,
     T <: AbstractTree{TT, RT, NL, N, B}} = _getnodename(tree, node)
getnodename(node::N) where
    {RT, NL, N <: AbstractNode{RT, NL}} = _getnodename(node)

"""
    hasbranch(tree::AbstractTree, branchname)


"""
function hasbranch(tree::AbstractTree{OneTree, RT, NL, N, B},
                   branch::Union{Int, B}) where {RT, NL, N, B}
    return _hasbranch(tree, branch)
end

"""
    getbranch(tree::AbstractTree, branchname)

Returns a branch from a tree
"""
function getbranch(tree::AbstractTree{OneTree, RT, NL, N, B},
                   branch::Union{Int, B}) where {RT, NL, N, B}
    _hasbranch(tree, branch) ||
        error("Branch $branch does not exist")
    return _getbranch(tree, branch)
end

"""
    getbranchname(::AbstractTree, id)
    getbranchname(id)

Returns the branch name associated with a branch from a tree. For some id and
branch types, it will be able to extract the branch name without reference to
the tree.
"""
function getbranchname end
getbranchname(tree::T, branch::Union{Int, B}) where
    {TT, RT, NL, N, B <: AbstractBranch{RT, NL},
     T <: AbstractTree{TT, RT, NL, N, B}} = _getbranchname(tree, branch)
getbranchname(branch::B) where
    {RT, NL, B <: AbstractBranch{RT, NL}} = _getbranchname(branch)

"""
    hasrootheight(tree::AbstractTree)


"""
function hasrootheight(tree::AbstractTree)
    return _hasrootheight(tree)
end

"""
    getrootheight(tree::AbstractTree)


"""
function getrootheight(tree::AbstractTree)
    return _getrootheight(tree)
end

"""
    setrootheight!(tree::AbstractTree, height)


"""
function setrootheight!(tree::AbstractTree, height)
    return _setrootheight!(tree, height)
end

# AbstractNode methods
"""
    isleaf(tree::AbstractTree, id)


"""
function isleaf end
isleaf(tree::AbstractTree, id) = _isleaf(tree, _getnode(tree, id))

"""
    isroot(node::AbstractNode)
    isroot(tree::AbstractTree, nodename)


"""
function isroot end
isroot(node::AbstractNode) = _isroot(node)
isroot(tree::AbstractTree, nodename) = _isroot(_getnode(tree, nodename))

"""
    isinternal(node::AbstractNode)
    isinternal(tree::AbstractTree, node)

Is the node internal to the tree (neither root nor leaf)?
"""
function isinternal end

function isinternal(node::AbstractNode)
    return _isinternal(node)
end

function isinternal(tree::AbstractTree, node)
    return _isinternal(_getnode(tree, node))
end

"""
    isunattached(node::AbstractNode)
    isunattached(tree::AbstractTree, node)

Is the node unattached (i.e. not connected to other nodes)?
"""
function isunattached end
isunattached(node::AbstractNode) = _isunattached(node)
isunattached(tree::AbstractTree, node) = _isunattached(_getnode(tree, node))

"""
    indegree(node::AbstractNode)
    indegree(tree::AbstractTree, node)

Return in degree of node
"""
function indegree end
indegree(node::AbstractNode) = _indegree(node)
indegree(tree::AbstractTree, node) = _indegree(_getnode(tree, node))

"""
    outdegree(node::AbstractNode)
    outdegree(tree::AbstractTree, node)

Return out degree of node
"""
function outdegree end
outdegree(node::AbstractNode) = _outdegree(node)
outdegree(tree::AbstractTree, node) = _outdegree(_getnode(tree, node))

"""
    degree(node::AbstractNode)
    degree(tree::AbstractTree, node)

Return the degree of a node including all connections
"""
function degree end
degree(node::AbstractNode) = _degree(node)
degree(tree::AbstractTree, node) = _degree(_getnode(tree, node))

"""
    hasoutboundspace(node::AbstractNode)
    hasoutboundspace(tree::AbstractTree, node)

Does the node have space for an[other] outbound connection?
"""
function hasoutboundspace end
hasoutboundspace(node::AbstractNode) = _hasoutboundspace(node)
hasoutboundspace(tree::AbstractTree, node) =
    _hasoutboundspace(_getnode(tree, node))

"""
    hasinbound(node::AbstractNode)
    hasinbound(tree::AbstractTree, node)

Does the node have an inbound connection?
"""
function hasinbound end
hasinbound(node::AbstractNode) = _hasinbound(node)
hasinbound(tree::AbstractTree, node) = _hasinbound(_getnode(tree, node))

"""
    hasinboundspace(node::AbstractNode)
    hasinboundspace(tree::AbstractTree, node)

Does the node have space for an inbound connection?
"""
function hasinboundspace end
hasinboundspace(node::AbstractNode) = _hasinboundspace(node)
hasinboundspace(tree::AbstractTree, node) =
    _hasinboundspace(_getnode(tree, node))

"""
    getinbound(node::AbstractNode)
    getinbound(tree::AbstractTree, node)

return the name of the inbound branch to this node.
"""
function getinbound end
getinbound(node::AbstractNode) = _getinbound(node)
getinbound(tree::AbstractTree, node) = _getinbound(_getnode(tree, node))

"""
    getparent(tree::AbstractTree, nodename)
    getparent(tree::AbstractTree, node)

Return [the name of] the parent node for this node [name].
"""
function getparent end
getparent(tree::AbstractTree{OneTree, RT, NL, N, B}, nodename::NL) where
    {RT <: Rooted, NL, N, B} =
    _getnodename(tree, _getparent(tree, _getnode(tree, nodename)))
getparent(tree::AbstractTree{OneTree, RT, NL, N, B}, node::N) where
    {RT, NL, N, B} = _getparent(tree, node)

"""
    getancestors(tree::AbstractTree, nodename)

Return the name of all of the nodes that are ancestral to this node.
"""
function getancestors(tree::AbstractTree, nodename)
    return _treepast(tree, nodename)[2][2:end]
end

"""
    getoutbounds(node::AbstractNode)
    getoutbounds(tree::AbstractTree, nodename)

Return the names of the outbound branches from this node.
"""
function getoutbounds end
getoutbounds(node::AbstractNode) = _getoutbounds(node)
getoutbounds(tree::AbstractTree, node) =
    _getoutbounds(_getnode(tree, nodename))

"""
    getchildren(tree::AbstractTree, node)
    getchildren(tree::AbstractTree, nodename)

Return the [name(s) of] the child node(s) for this node [name].
"""
function getchildren end
getchildren(tree::AbstractTree{OneTree, RT, NL, N, B}, nodename::NL) where
    {RT <: Rooted, NL, N, B} = map(branch -> _getnodename(dst(tree, branch)),
                                   getoutbounds(tree, nodename))
getchildren(tree::AbstractTree{OneTree, RT, NL, N, B}, node::N) where
    {RT <: Rooted, NL, N, B} = map(branch -> dst(tree, branch),
                                   getoutbounds(tree, node))

"""
    getdescendants(tree::AbstractTree, nodename)

Return the names of all of the nodes that descend from this node.
"""
function getdescendants(tree::AbstractTree, nodename)
    return _treefuture(tree, nodename)[2][2:end]
end

"""
    getconnections(tree::AbstractTree, node::AbstractNode)

Returns all of the branches connected to a node.
"""
function getconnections end
getconnections(tree::AbstractTree{OneTree, RT, NL, N, B}, nodename::NL) where
    {RT <: Rooted, NL, N, B} =
    map(branch -> _getbranchname(branch),
        _getconnections(tree, _getnode(tree, nodename)))
getconnections(tree::AbstractTree{OneTree, RT, NL, N, B}, node::N) where
    {RT <: Rooted, NL, N, B} = _getconnections(tree, node)

"""
    getsiblings(tree::AbstractTree, node::AbstractNode)

Returns all of the siblings of a node. Must be implemented for any unrooted
AbstractNode subtype, can be inferred from _getparent and _getchildren for
a rooted node.
"""
function getsiblings end
getsiblings(tree::AbstractTree{OneTree, RT, NL, N, B}, nodename::NL) where
    {RT <: Rooted, NL, N, B} = map(node -> _getnodename(node),
                                   _getsiblings(tree, _getnode(tree, nodename)))
getsiblings(tree::AbstractTree{OneTree, RT, NL, N, B}, node::N) where
    {RT <: Rooted, NL, N, B} = _getsiblings(tree, node)



"""
    hasheight(tree::AbstractTree, nodename)

Does the node have a height defined?
"""
function hasheight end

function hasheight(tree::AbstractTree, nodename)
    return _hasheight(tree, nodename) ||
        (_hasrootheight(tree) &&
         mapreduce(b -> haslength(tree, b), &, branchhistory(tree, nodename);
         init = _hasrootheight(tree)))
end

"""
    getheight(tree::AbstractTree, nodename)

Return the height of the node.
"""
function getheight(tree::AbstractTree, nodename)
    return _hasheight(tree, nodename) ? _getheight(tree, nodename) :
        mapreduce(b -> getlength(tree, b), +, branchhistory(tree, nodename);
                  init = getrootheight(tree))
end

"""
    setheight!(tree::AbstractTree, nodename, height)

Set the height of the node.
"""
function setheight!(tree::AbstractTree, nodename, height)
    return _setheight!(tree, nodename, height)
end

# AbstractBranch methods
"""
    src(branch::AbstractBranch)
    src(tree::AbstractTree, branchname)

Return the source node for this branch.
"""
function src end

function src(branch::AbstractBranch)
    return _src(branch)
end

function src(tree::AbstractTree, branchname)
    return _src(_getbranch(tree, branchname))
end

"""
    dst(branch::AbstractBranch)
    dst(tree::AbstractTree, branchname)

Return the destination node for this branch.
"""
function dst end

function dst(branch::AbstractBranch)
    return _dst(branch)
end

function dst(tree::AbstractTree, branchname)
    return _dst(_getbranch(tree, branchname))
end

"""
    Pair(branch::AbstractBranch)
    Pair(tree::AbstractTree, branchname)

Return a Pair containing the source and destination for this branch.
"""
function Pair end

function Pair(branch::AbstractBranch)
    return Pair(src(branch), dst(branch))
end

function Pair(tree::AbstractTree, branchname)
    return Pair(_getbranch(tree, branchname))
end

"""
    Tuple(branch::AbstractBranch)
    Tuple(tree::AbstractTree, branchname)

Return a Tuple containing the source and destination for this branch.
"""
function Tuple end

function Tuple(branch::AbstractBranch)
    return (src(branch), dst(branch))
end

function Tuple(tree::AbstractTree, branchname)
    return Tuple(_getbranch(tree, branchname))
end

"""
    getlength(branch::AbstractBranch)
    getlength(tree::AbstractTree, branchname)

Return the length of this branch.
"""
function getlength end

function getlength(branch::AbstractBranch)
    return _getlength(branch)
end

function getlength(tree::AbstractTree, branchname)
    return _getlength(_getbranch(tree, branchname))
end

"""
    resetleaves!(::AbstractTree)

Reset the leaf records to the current leaves, deleting all leaf records.
"""
resetleaves!(tree::AbstractTree) = _resetleaves!(tree)

"""
    getleafnames(::AbstractTree)

Retrieve the leaf names from the tree.
"""
getleafnames(tree::AbstractTree) = collect(_getleafnames(tree))

"""
    leafinfotype(::Type{<: AbstractTree})

retrieve the leaf info type of a tree.
"""
function leafinfotype(::Type{T}) where T <: AbstractTree
    return _leafinfotype(T)
end

"""
    nodeinfotype(::Type{<: AbstractTree})

retrieve the node info type of a tree.
"""
function nodeinfotype(::Type{T}) where T <: AbstractTree
    return _nodeinfotype(T)
end

"""
    branchinfotype(::Type{<: AbstractTree})

retrieve the branch info type of a tree.
"""
function branchinfotype(::Type{T}) where T <: AbstractTree
    return _branchinfotype(T)
end

"""
    getleafinfo(::AbstractTree, label)

retrieve the leaf info for a leaf of the tree.
"""
function getleafinfo(tree::AbstractTree)
    return _getleafinfo(tree)
end

"""
    setleafinfo!(::AbstractTree, table)

Set the leaf info for the leaves of the tree.
"""
function setleafinfo!(tree::AbstractTree, table)
    return _setleafinfo!(tree, table)
end

"""
    getnodeinfo(::AbstractTree, label)

retrieve the node info for a leaf of the tree.
"""
function getnodeinfo(tree::AbstractTree, label)
    return _getnodeinfo(tree, label)
end

"""
    setnodeinfo!(::AbstractTree, label, value)

Set the node info for a node of the tree.
"""
function setnodeinfo!(tree::AbstractTree, label, value)
    return _setnodeinfo!(tree, label, value)
end

"""
    getbranchinfo(::AbstractTree, label)

retrieve the branch info for a leaf of the tree.
"""
function getbranchinfo(tree::AbstractTree, label)
    return _getbranchinfo(tree, label)
end

"""
    setbranchinfo!(::AbstractTree, label, value)

Set the branch info for a branch of the tree.
"""
function setbranchinfo!(tree::AbstractTree, label, value)
    return _setbranchinfo!(tree, label, value)
end

"""
    validate(tree::AbstractTree)

Validate the tree by making sure that it is connected up correctly.
"""
function validate(tree::T) where
    {TT, RT, NL, N, B, T <: AbstractTree{TT, RT, NL, N, B}}
    nodes = _getnodes(tree)
    branches = _getbranches(tree)
    if !isempty(nodes) || !isempty(branches)
        # We need to validate the connections
        if Set(mapreduce(_getinbound, push!, nodefilter(_hasinbound, tree);
                         init = Int[])) != Set(keys(branches))
            warn("Inbound branches must exactly match Branch labels")
            return false
        end

        if Set(mapreduce(_getoutbounds, append!, nodeiter(tree);
                         init = Int[])) != Set(keys(branches))
            warn("Node outbound branches must exactly match Branch labels")
            return false
        end

        if !(mapreduce(_src, push!, branchiter(tree); init = NL[]) ⊆
             Set(keys(nodes)))
            warn("Branch sources must be node labels")
            return false
        end

        if !(mapreduce(_dst, push!, branchiter(tree); init = NL[]) ⊆
             Set(keys(nodes)))
            warn("Branch destinations must be node labels")
            return false
        end
    end

    return _validate(tree)
end
