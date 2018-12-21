using Phylo.API
using Compat: mapreduce
import LightGraphs: src, dst, indegree, outdegree, degree

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

"""
    leafinfotype(::Type{<: AbstractTree})

retrieve the leaf info type of a tree.
"""
leafinfotype(::Type{T}) where T <: AbstractTree = _leafinfotype(T)

"""
    nodeinfotype(::Type{<: AbstractTree})

retrieve the node info type of a tree.
"""
nodeinfotype(::Type{T}) where T <: AbstractTree = _nodeinfotype(T)

"""
    branchinfotype(::Type{<: AbstractTree})

retrieve the branch info type of a tree.
"""
branchinfotype(::Type{T}) where T <: AbstractTree = _branchinfotype(T)

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

"""
    nleaves(::AbstractTree)

Returns the number of leaves (tips) in a tree.
"""
nleaves(tree::AbstractTree) = _nleaves(tree)

"""
    nroots(::AbstractTree)
    nroots(::AbstractTree, label)

Returns the number of roots in a tree. For OneTree types, Unrooted trees will
return 0, OneRoot trees should return 1, and manyroots tree (ones with multiple
subtrees) will return the number of subtrees. manytree types will return a
vector of counts of the number of roots for each tree in the set.
"""
function nroots end
nroots(tree::AbstractTree, label) = _nroots(_getonetree(tree, label))
nroots(::AbstractTree{OneTree, Unrooted}) = 0
nroots(tree::AbstractTree{OneTree, RT}) where {RT <: Rooted} = _nroots(tree)
nroots(tree::AbstractTree{ManyTrees}) = [nroots(t) for t in _gettrees(tree)]

"""
    getroots(::AbstractTree)
    getroots(::AbstractTree, id)

Returns a vector containing the root(s) of a single (OneTree) tree
or a set of (ManyTrees) trees.
"""
function getroots end
getroots(tree::AbstractTree{OneTree}) = _getroots(tree)
getroots(tree::AbstractTree{OneTree, Unrooted, NL, N, B}) where {NL, N, B} = N[]
getroots(tree::AbstractTree{ManyTrees}) = [getroots(t) for t in _gettrees(tree)]


"""
    getroot(::AbstractTree)
    getroot(::AbstractTree, id)

Returns the root of a single tree (identified by id for a ManyTrees tree),
or a vector of roots for multiple trees.
"""
function getroot end
getroot(tree::AbstractTree{OneTree, OneRoot}) = _getroot(tree)
getroot(tree::AbstractTree{ManyTrees, OneRoot}, id) =
    getroot(_getonetree(tree, id))
getroot(tree::AbstractTree{ManyTrees, OneRoot}) =
    [getroot(t) for t in _gettrees(tree)]

"""
    getnodes(::AbstractTree)
    getnodes(::AbstractTree, id)

Returns the vector of nodes of a single tree (identified by id for a
ManyTrees tree), or a vector of vectors of nodes for multiple trees.
"""
function getnodes end
getnodes(tree::AbstractTree{OneTree}) = _getnodes(tree)
getnodes(tree::AbstractTree{ManyTrees}, id) = _getnodes(_getonetree(tree, id))
getnodes(tree::AbstractTree{ManyTrees}) = [getnodes(t) for t in _gettrees(tree)]

"""
    nnodes(::AbstractTree)
    nnodes(::AbstractTree, id)

Returns the number of nodes of a single tree (identified by id for a
ManyTrees tree), or a vector of numbers of nodes for multiple trees.
"""
function nnodes end
nnodes(tree::AbstractTree{OneTree}) = _nnodes(tree)
nnodes(tree::AbstractTree{ManyTrees}, id) = _nnodes(_getonetree(tree, id))
nnodes(tree::AbstractTree{ManyTrees}) = [nnodes(t) for t in _gettrees(tree)]

"""
    nbranches(::AbstractTree)
    nbranches(::AbstractTree, id)

Returns the number of branches of a single tree (identified by id for a
ManyTrees tree), or a vector of numbers of branches for multiple trees.
"""
function nbranches end
nbranches(tree::AbstractTree{OneTree}) = _nbranches(tree)
nbranches(tree::AbstractTree{ManyTrees}, id) = _nbranches(_getonetree(tree, id))
nbranches(tree::AbstractTree{ManyTrees}) =
    [nbranches(t) for t in _gettrees(tree)]

"""
    getnodenames(tree::AbstractTree)
    getnodenames(tree::AbstractTree, id)

Return a vector of node names of a single tree (identified by id for a
ManyTrees tree), or a vector of vectors of node names for multiple trees.
"""
function getnodenames end
getnodenames(tree::AbstractTree{OneTree}) = _getnodenames(tree)
getnodenames(tree::AbstractTree{ManyTrees}, id) =
    _getnodenames(_getonetree(tree, id))
getnodenames(tree::AbstractTree{ManyTrees}) =
    [_getnodenames(t) for t in _gettrees(tree)]

"""
    getbranches(::AbstractTree)
    getbranches(::AbstractTree, id)

Returns the vector of branches of a single tree (identified by id for a
ManyTrees tree), or a vector of vectors of branches for multiple trees.
"""
function getbranches end
getbranches(tree::AbstractTree{OneTree}) = _getbranches(tree)
getbranches(tree::AbstractTree{ManyTrees}, id) =
    _getbranches(_getonetree(tree, id))
getbranches(tree::AbstractTree{ManyTrees}) =
    [getbranches(t) for t in _gettrees(tree)]

"""
    getbranchnames(tree::AbstractTree)
    getbranchnames(tree::AbstractTree, id)

Return a vector of branch names of a single tree (identified by id for a
ManyTrees tree), or a vector of vectors of branch names for multiple trees.
"""
function getbranchnames end
getbranchnames(tree::AbstractTree{OneTree}) = _getbranchnames(tree)
getbranchnames(tree::AbstractTree{ManyTrees}, id) =
    _getbranchnames(_getonetree(tree, id))
getbranchnames(tree::AbstractTree{ManyTrees}) =
    [_getbranchnames(t) for t in _gettrees(tree)]

"""
    createbranch!(tree::AbstractTree, source, destination[, length::Float64];
                  data)

Add a branch from `source` to `destination` on `tree` with optional length and
data. `source` and `destination` can be either nodes or nodenames.
"""
function createbranch!(tree::AbstractTree{OneTree}, source, destination,
                       length::Float64 = NaN; data = missing)
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
    deletebranch!(tree::AbstractTree, src, dst)

Delete the branch `branch` from `tree`, or branch connecting `src` node to
`dst` node.
"""
function deletebranch! end
function deletebranch!(tree::AbstractTree{OneTree}, branch)
    _hasbranch(tree, branch) ||
        error("Tree does not have a branch called " *
              _getbranchname(tree, branch))
    return _deletebranch!(tree, branch)
end
function deletebranch!(tree::AbstractTree{OneTree, RT, NL, N, B},
                       src::Union{N, NL}, dst::Union{N, NL}) where
    {RT, NL, N, B}
    _hasbranch(tree, src, dst) ||
        error("Tree does not have a branch called " *
              _getbranchname(tree, branch))
    return _deletebranch!(tree, _getbranch(src, dst))
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
createnodes!(tree::AbstractTree{OneTree}, count::Int) =
    map(name -> createnode!(tree), Base.OneTo(count))
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
function deletenode!(tree::AbstractTree{OneTree}, node)
    return _deletenode!(tree, node)
end

"""
    hasnode(tree::AbstractTree, node)

Returns whether a tree has a given node (or node name) or not.
"""
function hasnode(tree::AbstractTree{OneTree}, node)
    return _hasnode(tree, node)
end

"""
    getnode(tree::AbstractTree, nodename)

Returns a node from a tree.
"""
function getnode(tree::AbstractTree{OneTree}, node)
    _hasnode(tree, node) ||
        error("Node $(_getnodename(tree, node)) does not exist")
    return _getnode(tree, node)
end

"""
    getnodename(::AbstractTree, node)
    getnodename(node)

Returns the node name associated with a node from a tree. For some
node types, it will be able to extract the node name without reference to
the tree.
"""
function getnodename end
getnodename(tree::AbstractTree{OneTree}, node) = _getnodename(tree, node)
getnodename(node::AbstractNode) = _getnodename(node)

"""
    hasbranch(tree::AbstractTree, branch)
    hasbranch(tree::AbstractTree, source, dest)

Does `tree` have a branch `branch` or a branch from `source` to `dest`?
"""
function hasbranch end
hasbranch(tree::AbstractTree{OneTree, RT, NL, N, B},
          branch::Union{Int, B}) where {RT, NL, N, B} = _hasbranch(tree, branch)
hasbranch(tree::AbstractTree{OneTree, RT, NL, N, B},
          source::Union{N, NL}, dest::Union{N, NL}) where {RT, NL, N, B} =
    _hasbranch(tree, source, dest)

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
    _hasbranch(tree, source, dest) ||
        error("Branch $branch does not exist")
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
getbranchname(branch::B) where
    {RT, NL, B <: AbstractBranch{RT, NL}} = _getbranchname(branch)

"""
    hasrootheight(tree::AbstractTree)


"""
function hasrootheight end
hasrootheight(::AbstractTree) = false
hasrootheight(tree::AbstractTree{OneTree, <:Rooted}) = _hasrootheight(tree)

"""
    getrootheight(tree::AbstractTree)


"""
function getrootheight(tree::AbstractTree{OneTree, <:Rooted})
    return _getrootheight(tree)
end

"""
    setrootheight!(tree::AbstractTree, height)


"""
function setrootheight!(tree::AbstractTree{OneTree, <:Rooted}, height)
    return _setrootheight!(tree, height)
end

# AbstractNode methods
"""
    isleaf(tree::AbstractTree, node)
    isleaf(node)

Is the node (referenced by name or node object) a leaf of the tree? Second
method may not be possible for some node types.
"""
function isleaf end
isleaf(tree::AbstractTree{OneTree}, node) = _isleaf(tree, node)
isleaf(node::AbstractNode) = _isleaf(node)

"""
    isroot(tree::AbstractTree, node)
    isroot(node::AbstractNode)

Is the node (referenced by name or node object) a root of the tree? Second
method may not be possible for some node types.
"""
function isroot end
isroot(tree::AbstractTree{OneTree}, node) = _isroot(tree, node)
isroot(node::AbstractNode) = _isroot(node)

"""
    isinternal(tree::AbstractTree, node)
    isinternal(node::AbstractNode)

Is the node (referenced by name or node object) internal to the tree (neither
root nor leaf)? Second method may not be implemented for some node types.
"""
function isinternal end
isinternal(tree::AbstractTree{OneTree}, node) = _isinternal(tree, node)
isinternal(node::AbstractNode) = _isinternal(node)

"""
    isunattached(tree::AbstractTree, node)
    isunattached(node::AbstractNode)

Is the node (referenced by name or node object) unattached (i.e. not connected
to other nodes)? Second method may not be implemented for some node types.
"""
function isunattached end
isunattached(tree::AbstractTree{OneTree}, node) = _isunattached(tree, node)
isunattached(node::AbstractNode) = _isunattached(node)

"""
    indegree(tree::AbstractTree, node)
    indegree(node::AbstractNode)

Return in degree of node. Second method may not be implemented for some node
types.
"""
function indegree end
indegree(tree::AbstractTree{OneTree}, node) = _indegree(tree, node)
indegree(node::AbstractNode) = _indegree(node)

"""
    outdegree(tree::AbstractTree, node)
    outdegree(node::AbstractNode)

Return out degree of node. Second method may not be implemented for some node
types.
"""
function outdegree end
outdegree(tree::AbstractTree{OneTree}, node) = _outdegree(tree, node)
outdegree(node::AbstractNode) = _outdegree(node)

"""
    degree(tree::AbstractTree, node)
    degree(node::AbstractNode)

Return the degree of a node including all connections. Second method may not
be implemented for some node types.
"""
function degree end
degree(tree::AbstractTree{OneTree}, node) = _degree(tree, node)
degree(node::AbstractNode) = _degree(node)

"""
    hasoutboundspace(tree::AbstractTree, node)
    hasoutboundspace(node::AbstractNode)

Does the node have space for an[other] outbound connection? Second method may
not be implemented for some node types.
"""
function hasoutboundspace end
hasoutboundspace(tree::AbstractTree{OneTree}, node) =
    _hasoutboundspace(tree, node)
hasoutboundspace(node::AbstractNode) = _hasoutboundspace(node)

"""
    hasinbound(tree::AbstractTree, node)
    hasinbound(node::AbstractNode)

Does the node have an inbound connection? Second method may not be implemented
for some node types.
"""
function hasinbound end
hasinbound(tree::AbstractTree{OneTree}, node) = _hasinbound(tree, node)
hasinbound(node::AbstractNode) = _hasinbound(node)

"""
    hasinboundspace(tree::AbstractTree, node)
    hasinboundspace(node::AbstractNode)

Does the node have space for an inbound connection? Second method may not be
implemented for some node types.
"""
function hasinboundspace end
hasinboundspace(tree::AbstractTree{OneTree}, node) =
    _hasinboundspace(tree, node)
hasinboundspace(node::AbstractNode) = _hasinboundspace(node)

"""
    getinbound(tree::AbstractTree, node)
    getinbound(node::AbstractNode)

return the inbound branch to this node (returns name for node name, branch for
node). Second method may not be implemented for some node types.
"""
function getinbound end
getinbound(tree::AbstractTree{OneTree, <:Rooted}, node) =
    _getinbound(tree, node)
getinbound(node::AbstractNode{<:Rooted}) = _getinbound(node)

"""
    getparent(tree::AbstractTree, node)
    getparent(node::AbstractNode)

Return [the name of] the parent node for this node [name]. Second method may
not be implemented for some node types.
"""
function getparent end
getparent(tree::AbstractTree{OneTree, <:Rooted}, node) =
    _getparent(tree, nodename)
getparent(node::AbstractNode{<:Rooted}) = _getparent(node)

"""
    getoutbounds(tree::AbstractTree, nodename)
    getoutbounds(node::AbstractNode)

Return the names of the outbound branches from this node. Second method may
not be implemented for some node types.
"""
function getoutbounds end
getoutbounds(tree::AbstractTree{OneTree, <:Rooted}, node) =
    _getoutbounds(tree, nodename)
getoutbounds(node::AbstractNode{<:Rooted}) = _getoutbounds(node)

"""
    getchildren(tree::AbstractTree, node)
    getchildren(node::AbstractNode)

Return the [name(s) of] the child node(s) for this node [name]. Second method
may not be implemented for some node types.
"""
function getchildren end
getchildren(tree::AbstractTree{OneTree, <:Rooted}, node) =
    _getchildren(tree, node)
getchildren(node::AbstractNode{<:Rooted}) = _getchildren(node)

"""
    getconnections(tree::AbstractTree, nodee)
    getconnections(node::AbstractNode)

Returns all of the branches connected to a node. Second method may
not be implemented for some node types.
"""
function getconnections end
getconnections(tree::AbstractTree{OneTree}, node) =
    _getconnections(tree, node)
getconnections(node::AbstractNode) = _getconnections(node)

"""
    getsiblings(tree::AbstractTree, node)
    getsiblings(node::AbstractNode)

Returns all of the siblings of a node. Must be implemented for any unrooted
AbstractNode subtype, can be inferred from _getparent and _getchildren for
a rooted node.
"""
function getsiblings end
getsiblings(tree::AbstractTree{OneTree}, node) =
    _getsiblings(tree, node)
getsiblings(node::AbstractNode) = _getsiblings(node)

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
    src(branch::AbstractBranch)

Return the source node for this branch.
"""
function src end
src(tree::AbstractTree{OneTree, <:Rooted}, branch) = _src(tree, branch)
src(branch::AbstractBranch{<:Rooted}) = _src(branch)

"""
    dst(tree::AbstractTree, branch)
    dst(branch::AbstractBranch)

Return the destination node for this branch.
"""
function dst end
dst(tree::AbstractTree{OneTree, <:Rooted}, branch) = _dst(tree, branch)
dst(branch::AbstractBranch{<:Rooted}) = _dst(branch)

"""
    link(tree::AbstractTree, branch, exclude)
    link(branch::AbstractBranch, exclude)

Return the node linked to `branch` that is not `exclude`.
"""
function link end
link(tree::AbstractTree{OneTree}, branch, exclude) =
    _link(tree, branch, exclude)
link(branch::AbstractBranch, exclude) = _link(branch, exclude)

"""
    getlength(tree::AbstractTree, branch)
    getlength(branch::AbstractBranch)

Return the length of this branch.
"""
function getlength end
getlength(tree::AbstractTree{OneTree}, branch) = _getlength(tree, branch)
getlength(branch::AbstractBranch) = _getlength(branch)

"""
    getleafnames(::AbstractTree)

Retrieve the leaf names from the tree.
"""
getleafnames(tree::AbstractTree) = collect(_getleafnames(tree))

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
    getnodeinfo(::AbstractTree, label)

retrieve the node info for a leaf of the tree.
"""
getnodeinfo(tree::AbstractTree{OneTree}, label) = _getnodeinfo(tree, label)

"""
    setnodeinfo!(::AbstractTree, label, value)

Set the node info for a node of the tree.
"""
setnodeinfo!(tree::AbstractTree{OneTree}, label, value) =
    _setnodeinfo!(tree, label, value)

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
