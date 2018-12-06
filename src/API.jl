using Phylo
using Phylo: Rootedness, Rooted, TreeType
using Phylo: AbstractNode, AbstractBranch, AbstractTree

# Need to be able to create new node and branch labels
function _newlabel(ids::Vector{Label}) where Label <: Integer
    return isempty(ids) ? 1 : maximum(ids) + 1
end
function _newlabel(ids::Vector{Vector{Label}}) where Label <: Integer
    return mapreduce(_newlabel, max, 1, ids)
end

function _newlabelfree(names::Vector{String}, prefix)
    names = collect(Iterators.filter(n -> length(n) > length(prefix) &&
                                     n[1:length(prefix)]==prefix, names))
    start = length(names) + 1
    name = prefix * "$start"
    while (name ∈ names)
        start += 1
        name = prefix * "$start"
    end
    return start
end
_newlabel(names::Vector{String}, prefix) =
    prefix * "$(_newlabelfree(names, prefix))"
_newlabel(names::Vector{Vector{String}}, prefix) =
    prefix * "$(mapreduce(_newlabelfree, max, 1, names))"

"""
    _newbranchlabel(tree::AbstractTree)

Returns a new unique branch name for a tree.
"""
function _newbranchlabel end
_newbranchlabel(tree::AbstractTree) = _newlabel(getbranchnames(tree))

"""
    _newnodelabel(tree::AbstractTree)

Returns a new unique node name for a tree.
"""
function _newnodelabel end
_newnodelabel(tree::AbstractTree{TT, RT, String, N, B}) where
    {TT, RT, N, B} = _newlabel(collect(getnodenames(tree)), "Node ")
_newnodelabel(tree::AbstractTree{TT, RT, I, N, B}) where
    {TT, RT, I <: Integer, N, B} = _newlabel(collect(getnodenames(tree)))

# AbstractTree type methods
"""
    _treenametype(::Type{AbstractTree})

Returns the label type for a tree type. Must be implemented for any tree type.
"""
function _treenametype end
_treenametype(::Type{T}) where
    {RT, NL, N, B, T <: AbstractTree{OneTree, RT, NL, N, B}} = Int

# AbstractTree methods
"""
    _ntrees(::AbstractTree)

Returns the number of trees in an object. Must be implemented for any ManyTrees
type.
"""
function _ntrees end
_ntrees(tree::AbstractTree{OneTree, RT, NL, N, B}) where
    {RT, NL, N, B} = 1
_ntrees(tree::AbstractTree{ManyTrees, RT, NL, N, B}) where
    {RT, NL, N, B} = length(_gettreenames(tree))

"""
    _gettreenames(::AbstractTree)

Returns the names for the trees. Must be implemented for any ManyTrees type,
should be implemented for any OneTree type where they have names.
"""
function _gettreenames end
_gettreenames(tree::AbstractTree{OneTree, RT, NL, N, B}) where
    {RT, NL, N, B} = [_gettreename(tree)]

"""
    _gettreenames(::AbstractTree)

Returns the names for the trees. Must be implemented for any ManyTrees type,
should be implemented for any OneTree type where they have names.
"""
function _gettreename end
_gettreename(tree::AbstractTree{OneTree, RT, NL, N, B}) where
    {RT, NL, N, B} = 1

"""
    _getonetree(::AbstractTree)
    _getonetree(::Pair{Label, AbstractTree})
    _getonetree(::AbstractTree, id)

Returns a tree - either itself if it is a single tree, or the single tree
in a set with label id. Must be implemented for any ManyTrees type.
"""
function _getonetree end
_getonetree(tree::AbstractTree{OneTree, RT, NL, N, B}) where
    {RT, NL, N, B} = tree
_getonetree(treepair::Pair{Label, T}) where
    {Label, RT, NL, N, B, T <: AbstractTree{OneTree, RT, NL, N, B}} = treepair[2]
function _getonetree(tree::T, id) where
    {RT, NL, N, B, T <: AbstractTree{OneTree, RT, NL, N, B}}
    id isa treenametype(T) || throw(TypeError(:_getonetree, "index argument",
                                               treenametype(T), id))
    id ∈ _gettreenames(tree) || throw(BoundsError(tree, id))
    return tree
end

"""
    _getnodes(tree::AbstractTree)

Returns a vector of nodes for a OneTree tree. Either _getnodes() or
_getnodenames() must be implemented for any OneTree tree type.
"""
function _getnodes end
_getnodes(tree::AbstractTree{OneTree, RT, NL, N, B}) where {RT, NL, N, B} =
    [_getnode(tree, name) for name in _getnodenames(tree)]

"""
    _getnodenames(tree::AbstractTree)

Returns a vector of node names for a OneTree tree. Either _getnodes() or
_getnodenames() must be implemented for any OneTree tree type.
"""
function _getnodenames end
_getnodenames(tree::AbstractTree{OneTree, RT, NL, N, B}) where {RT, NL, N, B} =
    [_getnodename(tree, node) for node in _getnodes(tree)]

"""
    _nnodes(::AbstractTree)

Returns the number of nodes (internal nodes and leaves) in a single tree. May
be implemented for any OneTree tree type (otherwise infers from _getnodes()).
"""
_nnodes(tree::AbstractTree{OneTree, RT, NL, N, B}) where {RT, NL, N, B} =
    length(_getnodes(tree))

"""
    _getleaves(::AbstractTree)

Returns the leaves (tips) of a single tree. May be implemented for any
OneTree type (otherwise determined from _getnodes() and _isleaf() functions).
"""
function _getleaves end
_getleaves(tree::AbstractTree{OneTree, RT, NL, N, B}) where {RT, NL, N, B} =
    filter(_isleaf, _getnodes(tree))

"""
    _nleaves(::AbstractTree)

Returns the number of leaves (tips) in a tree. May be implemented for any
OneTree type (otherwise determined from the _getleaves() function).
"""
function _nleaves end
_nleaves(tree::AbstractTree{ManyTrees, RT, NL, N, B}) where {RT, NL, N, B} =
    _nleaves(_getonetree(tree, _gettreenames(tree)[1]))
_nleaves(tree::AbstractTree{OneTree, RT, NL, N, B}) where {RT, NL, N, B} =
    length(_getleaves(tree))

"""
    _getroots(::AbstractTree)
    _getroots(::Rooted)

Returns the root(s) of a OneTree tree. May be implemented for any
OneTree type (otherwise determined from _getnodes() and _isroot() functions).
"""
function _getroots end
_getroots(tree::AbstractTree{ManyTrees, RT, NL, N, B}) where
    {RT, NL, N, B} = mapreduce(_getroots, append!, _gettrees(tree); init = N[])

"""
    _getroot(::AbstractTree)
    _getroot(::Rooted)

Returns the unique root of a rooted tree. May be implemented for any
OneTree type (otherwise determined from _getroots).
"""
function _getroot end
function _getroot(tree::AbstractTree{OneTree, RT, NL, N, B}) where
    {RT <: Rooted, NL, N, B}
    roots = _getroots(tree)
    @assert length(roots) == 1 "More than one root for tree ($(length(roots)))"
    return first(roots)
end

"""
    _nroots(::AbstractTree)

Returns the number of roots (subtrees) in a OneTree tree. May be implemented
for any ManyRoots type (otherwise infers from _getroots()).
"""
function _nroots end
_nroots(tree::AbstractTree{OneTree, RT, NL, N, B}) where
    {RT <: Rooted, NL, N, B} = length(_getroots(tree))
_nroots(tree::AbstractTree{OneTree, Unrooted, NL, N, B}) where {NL, N, B} = 0

"""
    _getnode(::AbstractTree, id)
    _getnode(id)

Returns the node associated with id (which could be a name, a node or a pair)
from a tree. For some id types, it will be possible to extract the node
without reference to the tree. Must be implemented for any tree and node label
type.
"""
function _getnode end
_getnode(::T, node::N) where
    {TT, RT, NL, N <: AbstractNode{RT, NL}, B,
     T <: AbstractTree{TT, RT, NL, N, B}} = node
_getnode(::T, pair::Pair{NL, N}) where
    {TT, RT, NL, N <: AbstractNode{RT, NL}, B,
     T <: AbstractTree{TT, RT, NL, N, B}} = pair[2]
_getnode(pair::Pair{NL, N}) where {RT, NL, N <: AbstractNode{RT, NL}} = pair[2]

"""
    _getnodename(::AbstractTree, id)
    _getnodename(id)

Returns the name of a node associated with id (which could be a name,
a node or a pair) from a tree. For some id and node types, it will be
able to extract the node name without reference to the tree.
"""
function _getnodename end
_getnodename(::T, pair::Pair{NL, N}) where
    {RT, NL, N <: AbstractNode{RT, NL}, B,
     T <: AbstractTree{OneTree, RT, NL, N, B}} = pair[1]
_getnodename(::T, nodename::NL) where
    {RT, NL, N, B,
     T <: AbstractTree{OneTree, RT, NL, N, B}} = nodename
_getnodename(pair::Pair{NL, N}) where
    {RT, NL, N <: AbstractNode{RT, NL}} = pair[1]


"""
    _getbranch(::AbstractTree, id)
    _getbranch(id)

Returns the branch associated with id (which could be a name, a branch or a
pair) from a tree. For some id types, it will be possible to extract the branch
without reference to the tree. Must be implemented for any OneTree tree type.
"""
function _getbranch end
_getbranch(::T, branch::B) where
    {TT, RT, NL, N, B <: AbstractBranch{RT, NL},
     T <: AbstractTree{TT, RT, NL, N, B}} = branch
_getbranch(::T, pair::Pair{B}) where
    {TT, RT, NL, N, B <: AbstractBranch{NL},
     T <: AbstractTree{TT, RT, NL, N, B}} = pair[2]
_getbranch(pair::Pair{B}) where
    {NL, B <: AbstractBranch{NL}} = pair[2]

"""
    _getbranchname(::AbstractTree, id)
    _getbranchname(id)

Returns the branch associated with id (which could be a name, a branch or a
pair) from a tree. For some id and branch types, it will be able to extract
the branch name without reference to the tree.
"""
function _getbranchname end
_getbranchname(::T, pair::Pair{B}) where
    {TT, RT, NL, N, B <: AbstractBranch{RT, NL},
     T <: AbstractTree{TT, RT, NL, N, B}} = pair[1]
_getbranchname(::T, branchname::Int) where
    {TT, RT, NL, N, B,
     T <: AbstractTree{TT, RT, NL, N, B}} = branchname
_getbranchname(pair::Pair{B}) where
    {RT, NL, B <: AbstractBranch{RT, NL}} = pair[1]

"""
    _hasnode(tree::AbstractTree, nodename)

May be implemented for any OneTree tree type.
"""
function _hasnode end
_hasnode(tree::AbstractTree{OneTree, RT, NL, N, B}, nodename::NL) where
    {RT, NL, N, B} = nodename ∈ _getnodenames(tree)

"""
    _getbranches(tree::AbstractTree)

Returns a vector of branches for a OneTree tree. Either _getbranches() or
_getbranchnames() must be implemented for any OneTree tree type.
"""
function _getbranches end
_getbranches(tree::AbstractTree{OneTree, RT, NL, N, B}) where
    {RT, NL, N, B} =
    [_getbranch(tree, name) for name in _getbranchnames(tree)]

"""
    _getbranchnames(tree::AbstractTree)

Returns a vector of branch names for a OneTree tree. Either _getbranches() or
_getbranchnames() must be implemented for any OneTree tree type.
"""
function _getbranchnames end
_getbranchnames(tree::AbstractTree{OneTree, RT, NL, N, B}) where
    {RT, NL, N, B} =
    [_getbranchname(tree, branch) for branch in _getbranches(tree)]

"""
    _hasbranch(tree::AbstractTree, branchname)

Tests whether a branch name is present in a tree. May be implemented for any
OneTree tree type.
"""
function _hasbranch end
_hasbranch(tree::AbstractTree{OneTree, RT, NL, N, B}, branchname::NL) where
    {RT, NL, N, B} = branchname ∈ _getbranchnames(tree)

"""
    _createbranch!(tree::AbstractTree, source, destination,
                   length::Float64, data)

Create a new branch and add it to a tree. Must be implemented for any
AbstractTree subtype.
"""
function _createbranch! end

"""
    _deletebranch!(tree::AbstractTree, branch)

Delete a branch, reoving it from a tree. Must be implemented for any
AbstractTree subtype.
"""
function _deletebranch! end

"""
    _createnode!(tree::AbstractTree, nodename)

Must be implemented for any AbstractTree subtype.
"""
function _createnode! end

"""
    _deletenode!(tree::AbstractTree, nodename)

Must be implemented for any AbstractTree subtype.
"""
function _deletenode! end

"""
    _hasrootheight(::AbstractTree)


"""
function _hasrootheight(::AbstractTree)
    return false
end

"""
    _getrootheight(::AbstractTree)


"""
function _getrootheight(::AbstractTree)
    throw(NullException())
    return NaN
end

"""
    _setrootheight!(::AbstractTree, value)


"""
function _setrootheight!(::AbstractTree, value)
    throw(NullException())
    return value
end
"""
    _validate(::AbstractTree)


"""
function _validate(::AbstractTree)
    return true
end


# AbstractNode methods
"""
    _isleaf(tree::AbstractTree, node::AbstractNode)

Is the node a leaf? Does not need to be implemented for any node type –
inferred from _outdegree or _degree - unless tree knows which nodes are
leaves and not nodes.
"""
function _isleaf end
_isleaf(::AbstractTree, node::AbstractNode) = _isleaf(node)
_isleaf(node::AbstractNode) = _outdegree(node) == 0
function _isleaf(node::AbstractNode{Unrooted, NL}) where NL
    _degree(node) == 0 && return true
    _degree(node) > 1 && return false
    return missing
end

"""
    _isroot(tree::AbstractTree, node::AbstractNode)

Is the node a root node of the tree?
"""
function _isroot end
_isroot(::AbstractTree, node::AbstractNode) = _isroot(node)
_isroot(::AbstractTree, node::AbstractNode{Unrooted, NL}) where NL = false
_isroot(tree::AbstractTree{TT, RT, NL, N, B}, node::N) where
    {TT, RT <: Rooted, NL, N, B} = node ∈ _getroots(tree)

"""
    _isinternal(node::AbstractNode)

Is the node internal to the tree?
"""
function _isinternal end
_isinternal(::AbstractTree, node::AbstractNode) = _isinternal(node)
_isinternal(node::AbstractNode{RT, NL}) where {RT <: Rooted, NL} =
    _outdegree(node) > 0 && _hasinbound(node)
_isinternal(node::AbstractNode{Unrooted, NL}) where NL =
    _degree(node) < 2 ? false : missing

"""
    _isunattached(node::AbstractNode)

Does the node currently form its own (sub)tree?
"""
_isunattached(node::AbstractNode) = _degree(node) == 0

"""
    _indegree(node::AbstractNode)

In degree of node. Can be implemented for rooted nodes, otherwise inferred
from _hasinbound.
"""
function _indegree end
_indegree(::AbstractTree, node::AbstractNode) = _indegree(node)
_indegree(node::AbstractNode) = Int(_hasinbound(node))
_indegree(node::AbstractNode{Unrooted, NL}) where NL =
    _degree(node) == 0 ? 0 : missing

"""
    _hasinboundspace(node::AbstractNode)

Is there space for a new inbound connection on a node?
"""
function _hasinboundspace end
_hasinboundspace(::AbstractTree, node::AbstractNode) = !_hasinbound(node)
_hasinboundspace(node::AbstractNode) = !_hasinbound(node)

"""
    _outdegree(node::AbstractNode)

Out degree of node.
"""
function _outdegree end
_outdegree(::AbstractTree, node::AbstractNode) = _outdegree(node)
_outdegree(node::AbstractNode) = length(_getoutbounds(node))
_outdegree(node::AbstractNode{Unrooted, NL}) where NL =
    _degree(node) == 0 ? 0 : missing

"""
    _hasoutboundspace(node::AbstractNode)

Is there space for a new outbound connection on a node? Must be implemented if
a node has a limit on the number of outbound connections (eg for a binary tree)
"""
function _hasoutboundspace end
_hasoutboundspace(::AbstractTree, node::AbstractNode) = _hasoutboundspace(node)
_hasoutboundspace(::AbstractNode) = true

"""
    _degree(node::AbstractNode)

Degree of node. Must be implemented for Unrooted nodes, otherwise can be
inferred from indegree and outdegree.
"""
function _degree end
_degree(::AbstractTree, node::AbstractNode) = _degree(node)
_degree(node::AbstractNode) = _indegree(node) + _outdegree(node)

"""
    _hasinbound(node::AbstractNode)

Must be implemented for any AbstractNode subtype.
"""
function _hasinbound end
_hasinbound(::AbstractTree, node::AbstractNode) = _hasinbound(node)
_hasinbound(node::AbstractNode{Unrooted, NL}) where NL =
    _degree(node) == 0 ? false : missing

"""
    _getinbound(node::AbstractNode)

Get the inbound connection. Must be implemented for any rooted AbstractNode subtype.
"""
function _getinbound end
_getinbound(::AbstractTree, node::AbstractNode) = _getinbound(node)
_getinbound(node::AbstractNode{Unrooted, NL}) where NL =
    error("Can't ask for inbound connection on Unrooted trees")

"""
    _getparent(tree::AbstractTree, node)

Return the parent node for this node. Can be implemented for Rooted node types.
"""
function _getparent end
_getparent(tree::AbstractTree{OneTree, Unrooted, NL, N, B}, nodename::NL) where
    {NL, N, B} = error("Can't ask for parent of Unrooted tree")
_getparent(tree::AbstractTree{OneTree, RT, NL, N, B}, node::N) where
    {RT <: Rooted, NL, N, B} = _src(tree, _getinbound(tree, node))

"""
    _addinbound!(node::AbstractNode, inbound)

Adds a branch to the input of a rooted node. Must be implemented for any rooted
AbstractNode subtype.
"""
function _addinbound! end
_addinbound!(node::AbstractNode{Unrooted, NL}, branch) where NL =
    error("Unrooted trees don't have inbound connections")

"""
    _removeinbound!(node::AbstractNode, inbound)

Removes a branch from the input of a rooted node. Must be implemented for any
rooted AbstractNode subtype.
"""
function _removeinbound! end
_removeinbound!(node::AbstractNode{Unrooted, NL}, branch) where NL =
    error("Unrooted trees don't have inbound connections")

"""
    _getoutbounds(node::AbstractNode)

Returns the outbound connections of a rooted node. Must be implemented for any
rooted AbstractNode subtype.
"""
function _getoutbounds end
_getoutbounds(node::AbstractNode{Unrooted, NL}) where NL =
    error("Unrooted trees don't have outbound connections")

"""
    _addoutbound!(node::AbstractNode, branch)

Add an outbound branch to a rooted node. Must be implemented for any Rooted
AbstractNode subtype.
"""
function _addoutbound! end
_addoutbound!(node::AbstractNode{Unrooted, NL}, branch) where NL =
    error("Unrooted trees don't have inbound and outbound connections")

"""
    _removeoutbound!(node::AbstractNode, branch)

Remove an outbound branch from a rooted node. Must be implemented for any
AbstractNode subtype.
"""
function _removeoutbound! end
_removeoutbound!(node::AbstractNode{Unrooted, NL}, branch) where NL =
    error("Unrooted trees don't have outbound connections")

"""
    _getchildren(tree::AbstractTree, node)
    _getchildren(tree::AbstractTree, nodename)

Return the child node(s) for this node. May be implemented for any rooted
AbstractNode subtype.
"""
function _getchildren end
_getchildren(tree::AbstractTree{OneTree, Unrooted, NL, N, B},
             nodename::NL) where {NL, N, B} =
    error("Can't ask for children of Unrooted tree")
_getchildren(tree::AbstractTree{OneTree, RT, NL, N, B}, node::N) where
    {RT <: Rooted, NL, N, B} = map(branch -> _dst(tree, branch),
                                   _getoutbounds(tree, node))

"""
    _getconnections(tree::AbstractTree, node::AbstractNode)

Returns all of the connections of a node. Must be implemented for any unrooted
AbstractNode subtype, can be inferred from _getinbound and _getoutbounds for
a rooted node.
"""
function _getconnections end
_getconnections(tree::AbstractTree, node::AbstractNode) =
    _hasinbound(node) ?
    append!([_getinbound(node)], _getoutbounds(node)) : _getoutbounds(node)

"""
    _getsiblings(tree::AbstractTree, node::AbstractNode)

Returns all of the siblings of a node. Must be implemented for any unrooted
AbstractNode subtype, can be inferred from _getparent and _getchildren for
a rooted node.
"""
function _getsiblings end
_getsiblings(node::AbstractNode) =
    _hasinbound(node) ? append!([_getparent(node)], _getchildren(node)) :
                        _getchildren(node)

"""
    _addconnection!(node::AbstractNode, branch)

Add a connection to an unrooted node. Must be implemented for any unrooted
AbstractNode subtype.
"""
function _addconnection! end
_addconnection!(node::AbstractNode{RT, NL}, branch) where {RT <: Rooted, NL} =
    error("Direction of connection must be specified for a rooted node")

"""
    _removeconnection!(node::AbstractNode, branch)

Remove a connection from an unrooted node. Must be implemented for any Unrooted
AbstractNode subtype.
"""
function _removeconnection! end
_removeconnection!(node::AbstractNode{RT, NL}, branch) where
    {RT <: Rooted, NL} =
    error("Direction of connection must be specified for a rooted node")

"""
    _hasheight(tree::AbstractTree, nodename)


"""
function _hasheight(::AbstractTree, _)
    return false
end

"""
    _getheight(tree::AbstractTree, nodename)


"""
function _getheight(::AbstractTree, _)
    throw(NullException())
    return NaN
end

"""
    _setheight!(::AbstractTree, nodename, value)


"""
function _setheight!(::AbstractTree, _, value)
    throw(NullException())
    return value
end

# AbstractBranch methods
"""
    _src(branch::AbstractBranch)

Return source node for a branch. Must be implemented for any rooted
AbstractBranch subtype.
"""
function _src end
_src(::AbstractTree, branch::AbstractBranch) = _src(branch)
_src(branch::B) where {NL, B <: AbstractBranch{Unrooted, NL}} =
    error("Unrooted branches do not have in and out connections")

"""
    _dst(branch::AbstractBranch)

Return destination node for a branch. Must be implemented for any rooted
AbstractBranch subtype.
"""
function _dst end
_dst(::AbstractTree, branch::AbstractBranch) = _dst(branch)
_dst(branch::B) where {NL, B <: AbstractBranch{Unrooted, NL}} =
    error("Unrooted branches do not have in and out connections")

"""
    _conns(branch::AbstractBranch[, exclude::AbstractNode])

Return a vector of connections for a branch, optionally excluding `exclude`
node. Must be implemented for any Unrooted AbstractBranch subtype, otherwise
can combine _src and _dst.
"""
function _conns end
_conns(branch::B) where {NL, RT <: Rooted, B <: AbstractBranch{RT, NL}} =
    [_src(branch), _dst(branch)]
_conns(branch::B, node::N) where
    {RT <: Rooted, NL, B <: AbstractBranch{RT, NL}, N <: AbstractNode{RT, NL}} =
    [_src(branch), _dst(branch)]

"""
    _getlength

Return length of a branch. May be implemented for any AbstractBranch
subtype.
"""
function _getlength end
function _getlength(::AbstractBranch)
    return NaN
end

#  - _getleafnames()
function _getleafnames(tree::AbstractTree)
    return collect(nodenamefilter(_isleaf, tree))
end

"""
    _leafinfotype(::Type{<:AbstractTree})

Returns the type of the leaf info data.
"""
function _leafinfotype end
_leafinfotype(::Type{<:AbstractTree}) = Void

"""
    _nodeinfotype(::Type{<:AbstractTree})

Returns the type of the node info data.
"""
function _nodeinfotype end
_nodeinfotype(::Type{<:AbstractTree}) = Void

"""
    _branchinfotype(::Type{<:AbstractTree})

Returns the type of the branch info data.
"""
function _branchinfotype end
_branchinfotype(::Type{<:AbstractTree}) = Void

function _addnode! end
function _removenode! end
function _addbranch! end
function _removebranch! end
function _getleafinfo end
function _setleafinfo! end
function _getnodeinfo end
function _setnodeinfo! end
function _getbranchinfo end
function _setbranchinfo! end
function _resetleaves! end
function _clearrootheight! end
