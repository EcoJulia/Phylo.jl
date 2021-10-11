using Phylo.API
import LightGraphs: src, dst, indegree, outdegree, degree
using SimpleTraits

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
leafinfotype(t::Type{<: AbstractTree}) = _leafinfotype(t)

"""
    nodedatatype(::Type{<: AbstractTree})

retrieve the node info type of a tree.
"""
nodedatatype(t::Type{<: AbstractTree}) = _nodedatatype(t)

"""
    branchdatatype(::Type{<: AbstractTree})

retrieve the branch info type of a tree.
"""
branchdatatype(t::Type{<: AbstractTree}) = _branchdatatype(t)

"""
    branchdims(::Type{<: AbstractTree})

retrieve the dimensions of the branch lengths for the tree.
"""
branchdims(t::Type{<: AbstractTree}) = _branchdims(t)

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
function gettree(trees::AbstractTree{ManyTrees})
    @assert _ntrees(trees) == 1 "Must be only one tree, found $(_ntrees(trees))"
    return first(_gettrees(tree))
end
gettree(tree::AbstractTree, label) = _gettree(tree, label)

"""
    gettreenames(tree::AbstractTree)

Returns the names of the trees.
"""
function gettreenames end
gettreenames(tree::AbstractTree{OneTree}) = [_gettreename(tree)]
gettreenames(trees::AbstractTree{ManyTrees}) = _gettreenames(trees)

"""
    gettreeinfo(tree::AbstractTree)
    gettreeinfo(tree::AbstractTree, treename)

Returns the info data associated with the tree(s).
"""
function gettreeinfo end
gettreeinfo(tree::AbstractTree) = _gettreeinfo(tree)
gettreeinfo(tree::AbstractTree, treename) = _gettreeinfo(tree, treename)

"""
    gettreename(tree::AbstractTree)

Returns the name of the single tree.
"""
gettreename(tree::AbstractTree{OneTree}) = _gettreename(tree)
function gettreename(trees::AbstractTree{ManyTrees})
    @assert _ntrees(trees) == 1 "Must be only one tree. found $(_ntrees(tree))"
    return first(_gettreenames(trees))
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
subtrees) will return the number of subtrees. ManyTrees types will return a
Dict of counts of the number of roots for each tree in the set.
"""
function nroots end
nroots(::AbstractTree{OneTree, Unrooted}) = 0
nroots(tree::AbstractTree{OneTree, <: Rooted}) = _nroots(tree)
nroots(trees::AbstractTree{ManyTrees}, name) = nroots(gettree(trees, name))
nroots(trees::AbstractTree{ManyTrees}) =
    Dict(name => nroots(gettree(trees, name)) for name in _gettreenames(trees))

"""
    getroots(::AbstractTree)
    getroots(::AbstractTree, id)

Returns a vector containing the root(s) of a single (OneTree) tree
or a set of (ManyTrees) trees.
"""
function getroots end
getroots(tree::AbstractTree{OneTree}) = _getroots(tree)
getroots(tree::AbstractTree{OneTree, Unrooted, NL, N}) where {NL, N} = N[]
getroots(trees::AbstractTree{ManyTrees}, name) = getroots(gettree(trees, name))
getroots(trees::AbstractTree{ManyTrees}) =
    Dict(name => getroots(gettree(trees, name))
         for name in _gettreenames(trees))

"""
    getroot(::AbstractTree)

Returns the root of a single tree (must be only one tree for a ManyTrees tree).
"""
function getroot end
getroot(tree::AbstractTree{OneTree, <: Rooted}) = _getroot(tree)
getroot(trees::AbstractTree{ManyTrees}, name) = getroot(gettree(trees, name))
getroot(trees::AbstractTree{ManyTrees}) =
    Dict(name => getroot(gettree(trees, name))
         for name in _gettreenames(trees))

"""
    getnodes(::AbstractTree[, ::TraversalOrder])

Returns the vector of nodes of a single tree, or a Dict of vectors of nodes
for multiple trees.
"""
function getnodes end
getnodes(tree::AbstractTree{OneTree}, order::TraversalOrder = preorder) =
    _getnodes(tree, order)
getnodes(trees::AbstractTree{ManyTrees}, name,
         order::TraversalOrder = preorder) =
             getnodes(gettree(trees, name), order)
getnodes(trees::AbstractTree{ManyTrees}, order::TraversalOrder = preorder) =
    Dict(name => getnodes(gettree(trees, name), order)
         for name in _gettreenames(trees))

"""
    nnodes(::AbstractTree)

Returns the number of nodes of a single tree, or a Dict of numbers of nodes
for multiple trees.
"""
function nnodes end
nnodes(tree::AbstractTree{OneTree}) = _nnodes(tree)
nnodes(trees::AbstractTree{ManyTrees}, name) = nnodes(gettree(trees, name))
nnodes(trees::AbstractTree{ManyTrees}) =
    Dict(name => nnodes(gettree(trees, name))
         for name in _gettreenames(trees))

"""
    ninternal(::AbstractTree)

Returns the number of internal nodes of a single tree, or a Dict of numbers of nodes
for multiple trees.
"""
ninternal(tree::AbstractTree) = nnodes(tree) - nleaves(tree)

"""
    nbranches(::AbstractTree)

Returns the number of branches of a single tree, or a Dict of numbers of
branches for multiple trees.
"""
function nbranches end
nbranches(tree::AbstractTree{OneTree}) = _nbranches(tree)
nbranches(trees::AbstractTree{ManyTrees}, name) =
    nbranches(gettree(trees, name))
nbranches(trees::AbstractTree{ManyTrees}) =
    Dict(name => nbranches(gettree(trees, name))
         for name in _gettreenames(trees))

"""
    getnodenames(::AbstractTree[, ::TraversalOrder])

Return a vector of node names of a single tree (identified by id for a
ManyTrees tree), or a Dict of vectors of node names for multiple trees.
"""
function getnodenames end
getnodenames(tree::AbstractTree{OneTree}, order::TraversalOrder = preorder) =
    _getnodenames(tree, order)
getnodenames(trees::AbstractTree{ManyTrees}, name,
             order::TraversalOrder = preorder) =
    getnodenames(gettree(trees, name), order)
getnodenames(trees::AbstractTree{ManyTrees}, order::TraversalOrder = preorder) =
    Dict(name => getnodenames(gettree(trees, name), order)
         for name in _gettreenames(trees))

"""
    getbranches(::AbstractTree)

Returns the vector of branches of a single tree, or a Dict of vectors of
branches for multiple trees.
"""
function getbranches end
getbranches(tree::AbstractTree{OneTree}) = _getbranches(tree)
getbranches(trees::AbstractTree{ManyTrees}, name) =
    getbranches(gettree(trees, name))
getbranches(trees::AbstractTree{ManyTrees}) =
    Dict(name => getbranches(gettree(trees, name))
         for name in _gettreenames(trees))

"""
    getbranchnames(tree::AbstractTree)

Return a vector of branch names of a single tree, or a Dict of vectors of
branch names for multiple trees.
"""
function getbranchnames end
getbranchnames(tree::AbstractTree{OneTree}) = _getbranchnames(tree)
getbranchnames(trees::AbstractTree{ManyTrees}, name) =
    getbranchnames(gettree(trees, name))
getbranchnames(trees::AbstractTree{ManyTrees}) =
    Dict(name => getbranchnames(gettree(trees, name))
         for name in _gettreenames(trees))

"""
    createbranch!(tree::AbstractTree, src, dst[, len::Number];
                  data)

Add a branch from `src` to `dst` on `tree` with optional length and
data. `source` and `destination` can be either nodes or nodenames.
"""
function createbranch! end
@traitfn function createbranch!(tree::T, src::N1, dst::N2, len = missing;
                                name = missing, data = missing) where
    {T <: AbstractTree{OneTree, <: Rooted}, N1, N2;
     !MatchNodeTypes{T, N1, N2}}
    hasnode(tree, src) ||
        error("Tree does not have an available source node - $src")
    hasnode(tree, dst) ||
        error("Tree does not have a destination node - $dst")
    sn = getnode(tree, src)
    dn = getnode(tree, dst)
    _hasoutboundspace(tree, sn) ||
        error(_getnodename(tree, sn) *
              " already has maximum number of outbound connections " *
              "($(_outdegree(tree, sn)))")
    _hasinboundspace(tree, dn) ||
        error("Tree does not have an available destination node called " *
              _getnodename(tree, dn))
    sn ≢ dn || error("Branch must connect different nodes")
    
    return ismissing(data) ?
        _createbranch!(tree, sn, dn, len, name) :
        _createbranch!(tree, sn, dn, len, name, data)
end
@traitfn function createbranch!(tree::T, src::N1, dst::N2, len = missing;
                                name = missing, data = missing) where
    {T <: AbstractTree{OneTree, <: Rooted}, N1, N2; MatchNodeTypes{T, N1, N2}}
    _hasnode(tree, src) ||
        error("Tree does not have an available source node - $src")
    _hasnode(tree, dst) ||
        error("Tree does not have a destination node - $dst")
    _hasoutboundspace(tree, src) ||
        error(_getnodename(tree, src) *
              " already has maximum number of outbound connections " *
              "($(_outdegree(tree, src)))")
    _hasinboundspace(tree, dst) ||
        error("Tree does not have an available destination node called " *
              _getnodename(tree, dst))
    src ≢ dst || error("Branch must connect different nodes")
    
    return ismissing(data) ?
        _createbranch!(tree, src, dst, len, name) :
        _createbranch!(tree, src, dst, len, name, data)
end
@traitfn function createbranch!(tree::T, src::N1, dst::N2, len = missing;
                                name = missing, data = missing) where
    {T <: AbstractTree{OneTree, Unrooted}, N1, N2;
     !MatchNodeTypes{T, N1, N2}}
    hasnode(tree, src) ||
        error("Tree does not have an available source node - $src")
    hasnode(tree, dst) ||
        error("Tree does not have a destination node - $dst")
    sn = getnode(tree, src)
    dn = getnode(tree, dst)
    _hasspace(tree, sn) ||
        error(_getnodename(tree, sn) *
              " already has maximum number of outbound connections " *
              "($(_outdegree(tree, sn)))")
    _hasspace(tree, dn) ||
        error("Tree does not have an available destination node called " *
              _getnodename(tree, dn))
    sn ≢ dn || error("Branch must connect different nodes")
    
    return ismissing(data) ?
        _createbranch!(tree, sn, dn, len, name) :
        _createbranch!(tree, sn, dn, len, name, data)
end
@traitfn function createbranch!(tree::T,
                                src::N1, dst::N2, len = missing;
                                name = missing, data = missing) where
    {T <: AbstractTree{OneTree, Unrooted}, N1, N2; MatchNodeTypes{T, N1, N2}}
    _hasnode(tree, src) ||
        error("Tree does not have an available source node - $src")
    _hasnode(tree, dst) ||
        error("Tree does not have a destination node - $dst")
    _hasspace(tree, src) ||
        error(_getnodename(tree, src) *
              " already has maximum number of outbound connections " *
              "($(_outdegree(tree, src)))")
    _hasspace(tree, dst) ||
        error("Tree does not have an available destination node called " *
              _getnodename(tree, dst))
    src ≢ dst || error("Branch must connect different nodes")
    
    return ismissing(data) ?
        _createbranch!(tree, src, dst, len, name) :
        _createbranch!(tree, src, dst, len, name, data)
end

"""
    deletebranch!(tree::AbstractTree, branch)
    deletebranch!(tree::AbstractTree, src, dst)

Delete the branch `branch` from `tree`, or branch connecting `src` node to
`dst` node.
"""
function deletebranch! end
@traitfn function deletebranch!(tree::T, branch::B) where
    {T <: AbstractTree{OneTree}, B; !MatchBranchType{T, B}}
    hasbranch(tree, branch) ||
        error("Tree does not have a branch $branch")
    return _deletebranch!(tree, getbranch(tree, branch))
end
@traitfn function deletebranch!(tree::T, branch::B) where
    {T <: AbstractTree{OneTree}, B; MatchBranchType{T, B}}
    _hasbranch(tree, branch) ||
        error("Tree does not have a branch called " *
              "$(_getbranchname(tree, branch))")
    return _deletebranch!(tree, branch)
end
@traitfn function deletebranch!(tree::T, src::N1, dst::N2) where
    {T <: AbstractTree{OneTree}, N1, N2; !MatchNodeTypes{T, N1, N2}}
    hasbranch(tree, src, dst) ||
        error("Tree does not have a branch between " *
              "$(_getnodename(tree, src)) and $(_getnodename(tree, dst))")
    return _deletebranch!(tree, getbranch(tree, src, dst))
end
@traitfn function deletebranch!(tree::T, src::N1, dst::N2) where
    {T <: AbstractTree{OneTree}, N1, N2; MatchNodeTypes{T, N1, N2}}
    hasbranch(tree, src, dst) ||
        error("Tree does not have a branch between " *
              "$(_getnodename(tree, src)) and $(_getnodename(tree, dst))")
    return _deletebranch!(tree, _getbranch(tree, src, dst))
end

"""
    createnode!(tree::AbstractTree[, nodename]; data)

Create a node on a tree with optional node info.
"""
function createnode!(tree::AbstractTree{OneTree, RT, NL, N, B},
                     nodename = missing; data = missing) where {RT, NL, N, B}
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
function createnodes!(tree::AbstractTree{OneTree, RT, NL, N, B},
                      nodenames::AbstractVector{NL}) where {RT, NL, N, B}
    here = [name for name in nodenames if _hasnode(tree, name)]
    isempty(here) || error("Nodes $here already present in tree")
    return [_createnode!(tree, name) for name in nodenames]
end
@traitfn function createnodes!(tree::T, nodedata::NDS) where
    {RT, NL, N, B, T <: AbstractTree{OneTree, RT, NL, N, B},
     ND, NDS <: Dict{NL, ND}; HoldsNodeData{T, ND}}
    here = [name for name in keys(nodedata) if _hasnode(tree, name)]
    isempty(here) || error("Nodes $here already present in tree")
    return [_createnode!(tree, info.first, info.second) for info in nodedata]
end

"""
    deletenode!(tree::AbstractTree, node)

Delete a node (or a name) from a tree
"""
function deletenode! end
@traitfn deletenode!(tree::T, node::NL) where
{NL, RT, T <: AbstractTree{OneTree, RT, NL}; !MatchNodeType{T, NL}} =
    deletenode!(tree, getnode(tree, node))
@traitfn function deletenode!(tree::T, node::N) where
    {T <: AbstractTree{OneTree}, N; MatchNodeType{T, N}}
    _hasnode(tree, node) || error("Tree does not have this node")
    return _deletenode!(tree, node)
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
    hasnode(tree, node) || error("Node $node does not exist")
    return _getnode(tree, node)
end

"""
    getnodename(::AbstractTree, node)

Returns the node name associated with a node from a tree. For some
node types, it will be able to extract the node name without reference to
the tree.
"""
function getnodename end
@traitfn function getnodename(tree::T, node::NL) where
    {NL, RT, T <: AbstractTree{OneTree, RT, NL}; !MatchNodeType{T, NL}}
    hasnode(tree, node) || error("Node $node does not exist")
    return _getnodename(tree, node)
end
@traitfn function getnodename(tree::T, node::N) where
    {T <: AbstractTree{OneTree}, N; MatchNodeType{T, N}}
    _hasnode(tree, node) || error("Node $node does not exist")
    return _getnodename(tree, node)
end

"""
    hasbranch(tree::AbstractTree, branch)
    hasbranch(tree::AbstractTree, source, dest)

Does `tree` have a branch `branch` or a branch from `source` to `dest`?
"""
function hasbranch end
@traitfn hasbranch(tree::T, branch::B) where
{T <: AbstractTree{OneTree}, B; MatchBranchType{T, B}} =
    _hasbranch(tree, branch)
@traitfn hasbranch(tree::T, branch::B) where
{T <: AbstractTree{OneTree}, B; !MatchBranchType{T, B}} =
    _hasbranch(tree, _getbranch(tree, branch))
@traitfn hasbranch(tree::T, src::N1, dst::N2) where
{T <: AbstractTree{OneTree}, N1, N2; !MatchNodeTypes{T, N1, N2}} =
    hasnode(tree, src) && hasnode(tree, dst) &&
    _hasbranch(tree, _getnode(tree, src), _getnode(tree, dst))
@traitfn hasbranch(tree::T, src::N1, dst::N2) where
{T <: AbstractTree{OneTree}, N1, N2; MatchNodeTypes{T, N1, N2}} =
    _hasnode(tree, src) && _hasnode(tree, dst) && _hasbranch(tree, src, dst)

"""
    getbranch(tree::AbstractTree, branch)
    getbranch(tree::AbstractTree, source, dest)

Returns a branch from a tree by name or by source and destination node.
"""
function getbranch end
@traitfn function getbranch(tree::T, branch::B) where
    {T <: AbstractTree{OneTree}, B; !MatchBranchType{T, B}}
    hasbranch(tree, branch) || error("Branch $branch does not exist")
    return _getbranch(tree, branch)
end
@traitfn function getbranch(tree::T, branch::B) where
    {T <: AbstractTree{OneTree}, B; MatchBranchType{T, B}}
    _hasbranch(tree, branch) || error("Branch $branch does not exist")
    return branch
end
@traitfn getbranch(tree::T, src::N1, dst::N2) where
{T <: AbstractTree{OneTree}, N1, N2; !MatchNodeTypes{T, N1, N2}} =
    getbranch(tree, getnode(tree, src), getnode(tree, dst))
@traitfn getbranch(tree::T, src::N1, dst::N2) where
{T <: AbstractTree{OneTree}, N1, N2; MatchNodeTypes{T, N1, N2}} =
    hasnode(tree, src) && hasnode(tree, dst) && _getbranch(tree, src, dst)

"""
    getbranchname(::AbstractTree, branch)
    getbranchname(branch)

Returns the branch name associated with a branch from a tree. For some
branch types, it will be able to extract the branch name without reference to
the tree.
"""
function getbranchname end
@traitfn function getbranchname(tree::T, branch::B) where
    {T <: AbstractTree{OneTree}, B; !MatchBranchType{T, B}}
    hasbranch(tree, branch) || error("Branch $branch does not exist")
    return _getbranchname(tree, branch)
end
@traitfn function getbranchname(tree::T, branch::B) where
    {T <: AbstractTree{OneTree}, B; MatchBranchType{T, B}}
    _hasbranch(tree, branch) || error("Branch $branch does not exist")
    return _getbranchname(tree, branch)
end
@traitfn function getbranchname(tree::T, branch::Int) where
    {T <: AbstractTree{OneTree}; MatchBranchType{T, Int}}
    hasbranch(tree, branch) || error("Branch $branch does not exist")
    return branch
end

"""
    hasrootheight(tree::AbstractTree)

Does the tree have an explicit root height
"""
function hasrootheight end
hasrootheight(tree::AbstractTree{OneTree, <: Rooted}) = _hasrootheight(tree)
hasrootheight(tree::AbstractTree{OneTree, Unrooted}) = false
hasrootheight(trees::AbstractTree{ManyTrees}) =
    all(hasrootheight(tree) for tree in gettrees(trees))

"""
    getrootheight(tree::AbstractTree)

Get the tree's root height.
"""
getrootheight(tree::AbstractTree{OneTree, <: Rooted}) = _getrootheight(tree)
getrootheight(trees::AbstractTree{ManyTrees}, name) =
    getrootheight(gettree(trees, name))
getrootheight(trees::AbstractTree{ManyTrees}) =
    Dict(name => getrootheight(gettree(trees, name))
         for name in _gettreenames(trees))

"""
    setrootheight!(tree::AbstractTree, height)

Set the tree's root height.
"""
function setrootheight! end
function setrootheight!(tree::AbstractTree{OneTree, <: Rooted}, height)
    return _setrootheight!(tree, height)
end
function setrootheight!(trees::AbstractTree{ManyTrees},
                        heights::AbstractVector)
    for (tree, height) in zip(gettrees(tree), heights)
        setrootheight!(tree, height)
    end
end

# AbstractNode methods
"""
    isleaf(tree::AbstractTree, node)

Is the node (referenced by name or node object) a leaf of the tree?
"""
function isleaf end
@traitfn function isleaf(tree::T, node::NL) where
    {NL, RT, T <: AbstractTree{OneTree, RT, NL}; !MatchNodeType{T, NL}}
    hasnode(tree, node) || error("Node $node does not exist")
    return _isleaf(tree, _getnode(tree, node))
end
@traitfn function isleaf(tree::T, node::N) where
    {T <: AbstractTree{OneTree}, N; MatchNodeType{T, N}}
    _hasnode(tree, node) || error("Node $node does not exist")    
    return _isleaf(tree, node)
end

"""
    isroot(tree::AbstractTree, node)

Is the node (referenced by name or node object) a root of the tree?
"""
function isroot end
@traitfn function isroot(tree::T, node::NL) where
    {NL, RT, T <: AbstractTree{OneTree, RT, NL}; !MatchNodeType{T, NL}}
    hasnode(tree, node) || error("Node $node does not exist")
    return _isroot(tree, _getnode(tree, node))
end
@traitfn function isroot(tree::T, node::N) where
    {T <: AbstractTree{OneTree}, N; MatchNodeType{T, N}}
    _hasnode(tree, node) || error("Node $node does not exist")
    return _isroot(tree, node)
end

"""
    isinternal(tree::AbstractTree, node)

Is the node (referenced by name or node object) internal to the tree (neither
root nor leaf)?
"""
function isinternal end
@traitfn function isinternal(tree::T, node::NL) where
    {NL, RT, T <: AbstractTree{OneTree, RT, NL}; !MatchNodeType{T, NL}}
    hasnode(tree, node) || error("Node $node does not exist")
    return _isinternal(tree, _getnode(tree, node))
end
@traitfn function isinternal(tree::T, node::N) where
    {T <: AbstractTree{OneTree}, N; MatchNodeType{T, N}}
    _hasnode(tree, node) || error("Node $node does not exist")
    return _isinternal(tree, node)
end

"""
    isunattached(tree::AbstractTree, node)

Is the node (referenced by name or node object) unattached (i.e. not connected
to other nodes)?
"""
function isunattached end
@traitfn function isunattached(tree::T, node::NL) where
    {NL, RT, T <: AbstractTree{OneTree, RT, NL}; !MatchNodeType{T, NL}}
    hasnode(tree, node) || error("Node $node does not exist")
    return _isunattached(tree, _getnode(tree, node))
end
@traitfn function isunattached(tree::T, node::N) where
    {T <: AbstractTree{OneTree}, N; MatchNodeType{T, N}}
    _hasnode(tree, node) || error("Node $node does not exist")
    return _isunattached(tree, node)
end

"""
    indegree(tree::AbstractTree, node)

Return in degree of node.
"""
function indegree end
@traitfn function indegree(tree::T, node::NL) where
    {NL, T <: AbstractTree{OneTree, <: Rooted, NL}; !MatchNodeType{T, NL}}
    hasnode(tree, node) || error("Node $node does not exist")
    return _indegree(tree, _getnode(tree, node))
end
@traitfn function indegree(tree::T, node::N) where
    {T <: AbstractTree{OneTree, <: Rooted}, N; MatchNodeType{T, N}}
    _hasnode(tree, node) || error("Node $node does not exist")
    return _indegree(tree, node)
end

"""
    outdegree(tree::AbstractTree, node)

Return out degree of node.
"""
function outdegree end
@traitfn function outdegree(tree::T, node::NL) where
    {NL, T <: AbstractTree{OneTree, <: Rooted, NL}; !MatchNodeType{T, NL}}
    hasnode(tree, node) || error("Node $node does not exist")
    return _outdegree(tree, _getnode(tree, node))
end
@traitfn function outdegree(tree::T, node::N) where
    {T <: AbstractTree{OneTree, <: Rooted}, N; MatchNodeType{T, N}}
    _hasnode(tree, node) || error("Node $node does not exist")
    return _outdegree(tree, node)
end

"""
    degree(tree::AbstractTree, node)

Return the degree of a node including all connections.
"""
function degree end
@traitfn function degree(tree::T, node::NL) where
    {NL, RT, T <: AbstractTree{OneTree, RT, NL}; !MatchNodeType{T, NL}}
    hasnode(tree, node) || error("Node $node does not exist")
    return _degree(tree, _getnode(tree, node))
end
@traitfn function degree(tree::T, node::N) where
    {T <: AbstractTree{OneTree}, N; MatchNodeType{T, N}}
    _hasnode(tree, node) || error("Node $node does not exist")
    return _degree(tree, node)
end

"""
    hasoutboundspace(tree::AbstractTree, node)

Does the node have space for an[other] outbound connection?
"""
function hasoutboundspace end
@traitfn function hasoutboundspace(tree::T, node::NL) where
    {NL, T <: AbstractTree{OneTree, <: Rooted, NL}; !MatchNodeType{T, NL}}
    hasnode(tree, node) || error("Node $node does not exist")
    return _hasoutboundspace(tree, _getnode(tree, node))
end
@traitfn function hasoutboundspace(tree::T, node::N) where
    {T <: AbstractTree{OneTree, <: Rooted}, N; MatchNodeType{T, N}}
    _hasnode(tree, node) || error("Node $node does not exist")
    return _hasoutboundspace(tree, node)
end

"""
    hasinbound(tree::AbstractTree, node)

Does the node have an inbound connection?
"""
function hasinbound end
@traitfn function hasinbound(tree::T, node::NL) where
    {NL, T <: AbstractTree{OneTree, <: Rooted, NL}; !MatchNodeType{T, NL}}
    hasnode(tree, node) || error("Node $node does not exist")
    return _hasinbound(tree, _getnode(tree, node))
end
@traitfn function hasinbound(tree::T, node::N) where
    {T <: AbstractTree{OneTree, <: Rooted}, N; MatchNodeType{T, N}}
    _hasnode(tree, node) || error("Node $node does not exist")
    return _hasinbound(tree, node)
end

"""
    hasinboundspace(tree::AbstractTree, node)

Does the node have space for an inbound connection?
"""
function hasinboundspace end
@traitfn function hasinboundspace(tree::T, node::NL) where
    {NL, T <: AbstractTree{OneTree, <: Rooted, NL}; !MatchNodeType{T, NL}}
    hasnode(tree, node) || error("Node $node does not exist")
    return _hasinboundspace(tree, _getnode(tree, node))
end
@traitfn function hasinboundspace(tree::T, node::N) where
    {T <: AbstractTree{OneTree, <: Rooted}, N; MatchNodeType{T, N}}
    _hasnode(tree, node) || error("Node $node does not exist")
    return _hasinboundspace(tree, node)
end

"""
    getinbound(tree::AbstractTree, node)

return the inbound branch to this node (returns name for node name, branch for
node).
"""
function getinbound end
@traitfn function getinbound(tree::T, node::NL) where
    {NL, T <: AbstractTree{OneTree, <: Rooted, NL}; !MatchNodeType{T, NL}}
    hasnode(tree, node) || error("Node $node does not exist")
    return _getinbound(tree, _getnode(tree, node))
end
@traitfn function getinbound(tree::T, node::N) where
    {T <: AbstractTree{OneTree, <: Rooted}, N; MatchNodeType{T, N}}
    _hasnode(tree, node) || error("Node $node does not exist")
    return _getinbound(tree, node)
end

"""
    getparent(tree::AbstractTree, node)

Return [the name of] the parent node for this node [name]. Second method may
not be implemented for some node types.
"""
function getparent end
@traitfn function getparent(tree::T, node::NL) where
    {NL, T <: AbstractTree{OneTree, <: Rooted, NL}; !MatchNodeType{T, NL}}
    hasnode(tree, node) || error("Node $node does not exist")
    return _getnodename(tree, _getparent(tree, _getnode(tree, node)))
end
@traitfn function getparent(tree::T, node::N) where
    {T <: AbstractTree{OneTree, <: Rooted}, N; MatchNodeType{T, N}}
    _hasnode(tree, node) || error("Node $node does not exist")
    return _getparent(tree, node)
end

"""
    getoutbounds(tree::AbstractTree, nodename)

Return the names of the outbound branches from this node.
"""
function getoutbounds end
@traitfn function getoutbounds(tree::T, node::NL) where
    {NL, T <: AbstractTree{OneTree, <: Rooted, NL}; !MatchNodeType{T, NL}}
    hasnode(tree, node) || error("Node $node does not exist")
    return _getoutbounds(tree, _getnode(tree, node))
end
@traitfn function getoutbounds(tree::T, node::N) where
    {T <: AbstractTree{OneTree, <: Rooted}, N; MatchNodeType{T, N}}
    _hasnode(tree, node) || error("Node $node does not exist")
    return _getoutbounds(tree, node)
end

"""
    getchildren(tree::AbstractTree, node)

Return the [name(s) of] the child node(s) for this node [name].
"""
function getchildren end
@traitfn function getchildren(tree::T, node::NL) where
    {NL, T <: AbstractTree{OneTree, <: Rooted, NL}; !MatchNodeType{T, NL}}
    hasnode(tree, node) || error("Node $node does not exist")
    return _getnodename.(tree, _getchildren(tree, _getnode(tree, node)))
end
@traitfn function getchildren(tree::T, node::N) where
    {T <: AbstractTree{OneTree, <: Rooted}, N; MatchNodeType{T, N}}
    _hasnode(tree, node) || error("Node $node does not exist")
    return _getchildren(tree, node)
end

"""
    getconnections(tree::AbstractTree, nodee)

Returns all of the branches connected to a node.
"""
function getconnections end
@traitfn function getconnections(tree::T, node::NL) where
    {NL, RT, T <: AbstractTree{OneTree, RT, NL}; !MatchNodeType{T, NL}}
    hasnode(tree, node) || error("Node $node does not exist")
    return _getconnections(tree, _getnode(tree, node))
end
@traitfn function getconnections(tree::T, node::N) where
    {T <: AbstractTree{OneTree}, N; MatchNodeType{T, N}}
    _hasnode(tree, node) || error("Node $node does not exist")
    return _getconnections(tree, node)
end

"""
    getsiblings(tree::AbstractTree, node)

Returns all of the siblings of a node. Must be implemented for any unrooted
AbstractNode subtype, can be inferred from _getparent and _getchildren for
a rooted node.
"""
function getsiblings end
@traitfn getsiblings(tree::T, node::NL) where
{NL, RT, T <: AbstractTree{OneTree, RT, NL}; !MatchNodeType{T, NL}} =
    _getsiblings(tree, getnode(tree, node))
@traitfn function getsiblings(tree::T, node::N) where
    {T <: AbstractTree{OneTree}, N; MatchNodeType{T, N}}
    _hasnode(tree, node) || error("Node $node does not exist")
    return _getsiblings(tree, node)
end

"""
    hasheight(tree::AbstractTree, node)

Does the node have a height defined?
"""
function hasheight end
@traitfn hasheight(tree::T, node::NL) where
{NL, T <: AbstractTree{OneTree, <: Rooted, NL}; !MatchNodeType{T, NL}} =
    hasheight(tree,  getnode(tree, node))
@traitfn function hasheight(tree::AbstractTree{OneTree, <: Rooted},
                            node::N) where
    {T <: AbstractTree{OneTree, <: Rooted}, N; MatchNodeType{T, N}}
    hasnode(tree, node) || error("Node $node does not exist")
    return _hasheight(tree, node) ||
        (_hasrootheight(tree) &&
         mapreduce(b -> haslength(tree, b), &, branchhistory(tree, node);
                   init = _hasrootheight(tree)))
end

"""
    getheight(tree::AbstractTree, node)

Return the height of the node.
"""
function getheight end
@traitfn getheight(tree::T, node::NL) where
{NL, T <: AbstractTree{OneTree, <: Rooted, NL}; !MatchNodeType{T, NL}} =
    getheight(tree,  getnode(tree, node))
@traitfn function getheight(tree::T, node::N) where
    {T <: AbstractTree{OneTree, <: Rooted}, N; MatchNodeType{T, N}}
    _hasnode(tree, node) || error("Node $node does not exist")
    return _hasheight(tree, node) ? _getheight(tree, node) :
        mapreduce(b -> getlength(tree, b), +, branchhistory(tree, node);
                  init = hasrootheight(tree) ?
                  getrootheight(tree) :
                  0.0 * upreferred(branchdims(T)))
end

"""
    setheight!(tree::AbstractTree, nodename, height)

Set the height of the node.
"""
function setheight! end
@traitfn setheight!(tree::T, node::NL, height::Number) where
{NL, T <: AbstractTree{OneTree, <: Rooted, NL}; !MatchNodeType{T, NL}} =
    _setheight!(tree, getnode(tree, node), height)
@traitfn function setheight!(tree::T, node::N, height::Number) where
    {T <: AbstractTree{OneTree, <: Rooted}, N; MatchNodeType{T, N}}
    _hasnode(tree, node) || error("Node $node does not exist")
    return _setheight!(tree, node, height)
end

# AbstractBranch methods
"""
    src(tree::AbstractTree, branch)

Return the source node for this branch.
"""
function src end
@traitfn src(tree::T, branch::B) where
{T <: AbstractTree{OneTree, <: Rooted}, B; !MatchBranchType{T, B}} =
    _src(tree, getbranch(tree, branch))
@traitfn function src(tree::T, branch::B) where
    {T <: AbstractTree{OneTree, <: Rooted}, B; MatchBranchType{T, B}}
    _hasbranch(tree, branch) || error("Branch $branch does not exist")
    return _src(tree, branch)
end

"""
    dst(tree::AbstractTree, branch)

Return the destination node for this branch.
"""
function dst end
@traitfn function dst(tree::T, branch::B) where
    {T <: AbstractTree{OneTree, <: Rooted}, B; !MatchBranchType{T, B}}
    hasbranch(tree, branch) || error("Branch $branch does not exist")
    return _dst(tree, getbranch(tree, branch))
end
@traitfn function dst(tree::T, branch::B) where
    {T <: AbstractTree{OneTree, <: Rooted}, B; MatchBranchType{T, B}}
    _hasbranch(tree, branch) || error("Branch $branch does not exist")
    return _dst(tree, branch)
end

"""
    conns(tree::AbstractTree, branch)

Return the nodes connected to `branch`.
"""
function conns end
@traitfn function conns(tree::T, branch::B) where
    {T <: AbstractTree{OneTree}, B; !MatchBranchType{T, B}}
    hasbranch(tree, branch) || error("Branch $branch does not exist")
    return _conns(tree, getbranch(tree, branch))
end
@traitfn function conns(tree::T, branch::B) where
    {T <: AbstractTree{OneTree}, B; MatchBranchType{T, B}}
    _hasbranch(tree, branch) || error("Branch $branch does not exist")
    return _conns(tree, branch)
end

"""
    conn(tree::AbstractTree, branch, exclude)

Return the other node connected to `branch` that is not `exclude`.
"""
function conn end
@traitfn function conn(tree::T, branch::B, exclude::N) where
    {T <: AbstractTree{OneTree}, B, N; !MatchBranchNodeType{T, B, N}}
    hasbranch(tree, branch) || error("Branch $branch does not exist")
    hasnode(tree, exclude) || error("Node $exclude does not exist")
    return _conn(tree, _getbranch(tree, branch), _getnode(tree, exclude))
end
@traitfn function conn(tree::T, branch::B, exclude::N) where
    {T <: AbstractTree{OneTree}, B, N; MatchBranchNodeType{T, B, N}}
    _hasbranch(tree, branch) || error("Branch $branch does not exist")
    _hasnode(tree, exclude) || error("Node $exclude does not exist")
    return _conn(tree, branch, exclude)
end

"""
    getlength(tree::AbstractTree, branch)

Return the length of this branch.
"""
function getlength end
@traitfn function getlength(tree::T, branch::B) where
    {T <: AbstractTree{OneTree}, B, N; !MatchBranchType{T, B}}
    hasbranch(tree, branch) || error("Branch $branch does not exist")
    return _getlength(tree, _getbranch(tree, branch))
end
@traitfn function getlength(tree::T, branch::B) where
    {T <: AbstractTree{OneTree}, B, N; MatchBranchType{T, B}}
    _hasbranch(tree, branch) || error("Branch $branch does not exist")
    return _getlength(tree, branch)
end

"""
    getleafnames(::AbstractTree[, ::TraversalOrder])

Retrieve the leaf names from the tree (in some specific order).
"""
getleafnames(tree::AbstractTree, order::TraversalOrder = preorder) =
    _getleafnames(tree, order)

"""
    getleaves(::AbstractTree[, ::TraversalOrder])

Retrieve the leaves from the tree.
"""
getleaves(tree::AbstractTree{OneTree}, order::TraversalOrder = preorder) =
    _getleaves(tree, order)
getleaves(trees::AbstractTree{ManyTrees}, name,
          order::TraversalOrder = preorder) =
    _getleaves(gettree(trees, name), order)
getleaves(trees::AbstractTree{ManyTrees}, order::TraversalOrder = preorder) =
    Dict(name => _getleaves(gettree(trees, name), order)
         for name in _gettreenames(trees))

"""
    getleafinfo(::AbstractTree[, label])

retrieve the leaf info for a leaf of the tree.
"""
function getleafinfo end
getleafinfo(tree::AbstractTree) = _getleafinfo(tree)
@traitfn getleafinfo(tree::T, leaf::NL) where
{RT, NL, T <: AbstractTree{OneTree, RT, NL}; !MatchNodeType{T, NL}} =
    _getleafinfo(tree, getnode(tree, leaf))
@traitfn function getleafinfo(tree::T, leaf::N) where
    {T <: AbstractTree{OneTree}, N; MatchNodeType{T, N}}
    _hasnode(tree, leaf) || error("Node $node does not exist")
    return _getleafinfo(tree, leaf)
end

"""
    setleafinfo!(::AbstractTree, table)

Set the leaf info for the leaves of the tree.
"""
function setleafinfo!(tree::AbstractTree, table)
    _setleafinfo!(tree, table)
    validate!(tree) || error("Leaf info not consistent with tree")
end

"""
    getnodedata(::AbstractTree, node)

retrieve the node data for a node of the tree.
"""
function getnodedata end
@traitfn getnodedata(tree::T, node::NL) where
{NL, RT, T <: AbstractTree{OneTree, RT, NL}; !MatchNodeType{T, NL}} =
    _getnodedata(tree, getnode(tree, node))
@traitfn function getnodedata(tree::T, node::N) where
    {T <: AbstractTree{OneTree}, N; MatchNodeType{T, N}}
    _hasnode(tree, node) || error("Node $node does not exist")
    return _getnodedata(tree, node)
end
@traitfn getnodedata(tree::T, node::NL, label) where
{NL, RT, T <: AbstractTree{OneTree, RT, NL}; !MatchNodeType{T, NL}} =
    _getnodedata(tree, getnode(tree, node), label)
@traitfn function getnodedata(tree::T, node::N, label) where
    {T <: AbstractTree{OneTree}, N; MatchNodeType{T, N}}
    _hasnode(tree, node) || error("Node $node does not exist")
    return _getnodedata(tree, node, label)
end

"""
    setnodedata!(::AbstractTree, node, label, value)
    setnodedata!(::AbstractTree, node, data)

Set the node data for a node of the tree.
"""
function setnodedata! end
@traitfn setnodedata!(tree::T, node::NL, label, value) where
{NL, RT, T <: AbstractTree{OneTree, RT, NL}; !MatchNodeType{T, NL}} =
    _setnodedata!(tree, getnode(tree, node), label, value)
@traitfn function setnodedata!(tree::T, node::N, label, value) where
    {T <: AbstractTree{OneTree}, N; MatchNodeType{T, N}}
    _hasnode(tree, node) || error("Node $node does not exist")
    return _setnodedata!(tree, node, label, value)
end
@traitfn setnodedata!(tree::T, node::NL, data) where
{NL, RT, T <: AbstractTree{OneTree, RT, NL}; !MatchNodeType{T, NL}} =
    _setnodedata!(tree, getnode(tree, node), data)
@traitfn function setnodedata!(tree::T, node::N, data) where
    {T <: AbstractTree{OneTree}, N; MatchNodeType{T, N}}
    _hasnode(tree, node) || error("Node $node does not exist")
    return _setnodedata!(tree, node, data)
end

"""
    getbranchdata(::AbstractTree, label)

retrieve the branch data for a leaf of the tree.
"""
function getbranchdata end
@traitfn function getbranchdata(tree::T, branch::B) where
    {T <: AbstractTree{OneTree}, B; !MatchBranchType{T, B}}
    hasbranch(tree, branch) || error("Branch $branch does not exist")
    return _getbranchdata(tree, _getbranch(tree, branch))
end
@traitfn function getbranchdata(tree::T, branch::B) where
    {T <: AbstractTree{OneTree}, B; MatchBranchType{T, B}}
    _hasbranch(tree, branch) || error("Branch $branch does not exist")
    return _getbranchdata(tree, branch)
end

"""
    setbranchdata!(::AbstractTree, branch, label, value)
    setbranchdata!(::AbstractTree, branch, data)

Set the branch data for a branch of the tree.
"""
function setbranchdata! end
@traitfn function setbranchdata!(tree::T, branch::B, label, value) where
    {T <: AbstractTree{OneTree}, B; !MatchBranchType{T, B}}
    hasbranch(tree, branch) || error("Branch $branch does not exist")
    return _setbranchdata!(tree, _getbranch(tree, branch), label, value)
end
@traitfn function setbranchdata!(tree::T, branch::B, label, value) where
    {T <: AbstractTree{OneTree}, B; MatchBranchType{T, B}}
    _hasbranch(tree, branch) || error("Branch $branch does not exist")
    return _setbranchdata!(tree, branch, label, value)
end
@traitfn function setbranchdata!(tree::T, branch::B, data) where
    {T <: AbstractTree{OneTree}, B; !MatchBranchType{T, B}}
    hasbranch(tree, branch) || error("Branch $branch does not exist")
    return _setbranchdata!(tree, _getbranch(tree, branch), data)
end
@traitfn function setbranchdata!(tree::T, branch::B, data) where
    {T <: AbstractTree{OneTree}, B; MatchBranchType{T, B}}
    _hasbranch(tree, branch) || error("Branch $branch does not exist")
    return _setbranchdata!(tree, branch, data)
end

"""
    validate!(tree::AbstractTree)

Validate the tree by making sure that it is connected up correctly.
"""
function validate!(tree::T) where
    {TT, RT, NL, N, B, T <: AbstractTree{TT, RT, NL, N, B}}
    _resetleaves!(tree)
    nodes = _getnodes(tree)
    nodenames = _getnodenames(tree)
    branches = _getbranches(tree)
    branchnames = _getbranchnames(tree)
    if !isempty(nodes) || !isempty(branches)
        # We need to validate the connections
        if Set(_getinbound(tree, node) for node in nodes
               if _hasinbound(tree, node)) != Set(branches)
            @warn "Inbound branches must exactly match Branch labels"
            return false
        end

        if Set(mapreduce(node -> _getoutbounds(tree, node), append!,
                         nodes; init = B[])) != Set(branches)
            @warn "Node outbound branches must exactly match Branch labels"
            return false
        end

        if !(Set(_getnodename(tree, _src(tree, branch))
                 for branch in branches) ⊆ Set(nodenames))
            @warn "Branch sources must be nodes"
            return false
        end

        if !(Set(_getnodename(tree, _dst(tree, branch))
                for branch in branches) ⊆ Set(nodenames))
            @warn "Branch destinations must be nodes"
            return false
        end
    end

    return _validate!(tree)
end

"""
    traversal(::AbstractTree, ::TraversalOrder)
    traversal(::AbstractTree, ::TraversalOrder, init)

Return an iterable object for a tree containing nodes in given order -
preorder, inorder, postorder or breadthfirst - optionally starting from init.
"""
function traversal end
traversal(tree::AbstractTree{OneTree}, order::TraversalOrder = preorder) =
    _traversal(tree, order)
@traitfn function traversal(tree::T, order::TraversalOrder, init::NL) where
    {NL, RT, T <: AbstractTree{OneTree, RT, NL}; !MatchNodeType{T, NL}}
    hasnode(tree, init) || error("Node $init does not exist")
    return _getnodename.(tree, _traversal(tree, order, [getnode(tree, init)]))
end
@traitfn function traversal(tree::T, order::TraversalOrder, init::N) where
    {T <: AbstractTree{OneTree}, N; MatchNodeType{T, N}}
    _hasnode(tree, init) || error("Node $init does not exist")
    return _traversal(tree, order, [init])
end

"""
    getancestors(tree::AbstractTree, node)

Return the name of all of the nodes that are ancestral to this node.
"""
function getancestors end
@traitfn function getancestors(tree::T, init::NL) where
    {NL, T <: AbstractTree{OneTree, <: Rooted, NL}; !MatchNodeType{T, NL}}
    hasnode(tree, init) || error("Node $init does not exist")
    return getnodename.(tree, _treehistory(tree, getnode(tree, init))[2][2:end])
end
@traitfn function getancestors(tree::T, init::N) where
    {T <: AbstractTree{OneTree, <: Rooted}, N; MatchNodeType{T, N}}
    hasnode(tree, init) || error("Node $init does not exist")
    return _treehistory(tree, init)[2][2:end]
end

"""
    getdescendants(tree::AbstractTree, node)

Return the names of all of the nodes that descend from this node.
"""
function getdescendants end
@traitfn function getdescendants(tree::T, init::NL) where
    {NL, T <: AbstractTree{OneTree, <: Rooted, NL}; !MatchNodeType{T, NL}}
    hasnode(tree, init) || error("Node $init does not exist")
    return getnodename.(tree, _treefuture(tree, getnode(tree, init))[2][2:end])
end
@traitfn function getdescendants(tree::T, init::N) where
    {T <: AbstractTree{OneTree, <: Rooted}, N; MatchNodeType{T, N}}
    hasnode(tree, init) || error("Node $init does not exist")
    return _treefuture(tree, init)[2][2:end]
end
