using Phylo
using Phylo: Rootedness, Rooted, TreeType
using Phylo: AbstractNode, AbstractBranch, AbstractTree
using Compat

# Need to be able to create new node and branch labels
function _newlabel(ids::Vector{Label}) where Label <: Integer
    return isempty(ids) ? 1 : maximum(ids) + 1
end
function _newlabel(ids::Vector{Vector{Label}}) where Label <: Integer
    return mapreduce(_newlabel, max, 1, ids)
end
_newlabel(names::Vector{String}, prefix::String) =
    prefix * "$(_newlabelfree(names, prefix))"
_newlabel(names::Vector{Vector{String}}, prefix::String) =
    prefix * "$(mapreduce(_newlabelfree, max, 1, names))"

function _newlabelfree(names::Vector{String}, prefix)
    num = length(names) + 1
    if "$prefix$num" ∉ names
        return num
    end
    len = length(prefix)
    goodnames = [name[(len - 1):end] for name in names
             if length(name) > len && name[1:len] == prefix]
    num = length(goodnames) + 1
    while ("$num" ∈ goodnames)
        num += 1
    end
    return num
end

"""
    _newbranchlabel(tree::AbstractTree)

Returns a new unique branch name for a tree.
"""
function _newbranchlabel end
_newbranchlabel(tree::AbstractTree) = _newlabel(collect(getbranchnames(tree)))

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
_treenametype(::Type{<: AbstractTree{OneTree}}) = Int

# AbstractTree methods
"""
    _ntrees(::AbstractTree)

Returns the number of trees in an object. Must be implemented for any ManyTrees
type.
"""
function _ntrees end
_ntrees(tree::AbstractTree{OneTree}) = 1
_ntrees(tree::AbstractTree{ManyTrees}) = length(_gettreenames(tree))

"""
    _gettrees(::AbstractTree)

Returns the trees in an object. Must be implemented for any ManyTrees type.
"""
function _gettrees end
_gettrees(tree::AbstractTree{OneTree}) = [tree]

"""
    _gettreenames(::AbstractTree)

Returns the names for the trees. Can be implemented for any ManyTrees type.
"""
function _gettreenames end
_gettreenames(tree::AbstractTree{ManyTrees}) = _gettreename.(_gettrees(tree))

"""
    _gettreename(::AbstractTree)

Returns the name for a single tree. Should be implemented for any OneTree
type where they have names.
"""
function _gettreename end
_gettreename(tree::AbstractTree{OneTree}) = 1

"""
    _gettree(::Pair{Label, AbstractTree})
    _gettree(::AbstractTree, id)

Returns a tree - either itself if it is a single tree, or the single tree
in a set with label id. Must be implemented for any ManyTrees type.
"""
function _gettree end
_gettree(treepair::Pair{Label, <: AbstractTree{OneTree}}) where Label =
    treepair[2]
function _gettree(tree::T, id) where {T <: AbstractTree{OneTree}}
    id isa _treenametype(T) || throw(TypeError(:_gettree, "index argument",
                                               _treenametype(T), id))
    id == _gettreename(tree) || throw(BoundsError(tree, id))
    return tree
end

"""
    _getnodes(tree::AbstractTree)

Returns a vector of nodes for a OneTree tree. Either _getnodes() or
_getnodenames() must be implemented for any OneTree tree type.
"""
function _getnodes end
_getnodes(tree::AbstractTree{OneTree}) =
    [_getnode(tree, name) for name in _getnodenames(tree)]

"""
    _getnodenames(tree::AbstractTree)

Returns a vector of node names for a OneTree tree. Either _getnodes() or
_getnodenames() must be implemented for any OneTree tree type.
"""
function _getnodenames end
_getnodenames(tree::AbstractTree{OneTree}) =
    [_getnodename(tree, node) for node in _getnodes(tree)]

"""
    _nnodes(::AbstractTree)

Returns the number of nodes (internal nodes and leaves) in a single tree. May
be implemented for any OneTree tree type (otherwise infers from _getnodes()).
"""
_nnodes(tree::AbstractTree{OneTree}) = length(_getnodes(tree))

"""
    _nbranches(::AbstractTree)

Returns the number of branches in a single tree. May be implemented for any OneTree tree type (otherwise infers from _getbranches()).
"""
_nbranches(tree::AbstractTree{OneTree}) = length(_getbranches(tree))

"""
    _getleafnames(::AbstractTree)

Returns the leaf names of a tree. May be implemented for any
tree type (otherwise determined from _getnodenames() and _isleaf() functions).
"""
function _getleafnames end
_getleafnames(tree::AbstractTree{OneTree}) =
    [node for node in _getnodenames(tree) if _isleaf(tree, node)]
_getleafnames(tree::AbstractTree{ManyTrees}) =
    _getleafnames(_gettree(tree, first(_gettreenames(tree))))

"""
    _getleaves(::AbstractTree)

Returns the leaves (tips) of a single tree. May be implemented for any
OneTree type (otherwise determined from _getnodes() and _isleaf() functions).
"""
function _getleaves end
_getleaves(tree::AbstractTree{OneTree}) =
    [node for node in _getnodes(tree) if _isleaf(tree, node)]

"""
    _nleaves(::AbstractTree)

Returns the number of leaves (tips) in a tree. May be implemented for any
tree type (otherwise determined from the _getleafnames() function).
"""
_nleaves(tree::AbstractTree) = length(_getleafnames(tree))

"""
    _getroots(::AbstractTree)

Returns the root(s) of a tree. May be implemented for any
OneTree type (otherwise determined from _getnodes() and _isroot() functions).
"""
_getroots(tree::AbstractTree{OneTree, <: Rooted}) =
    [node for node in _getnodes(tree) if _isroot(tree, node)]

"""
    _getroot(::AbstractTree)

Returns the unique root of a rooted tree. May be implemented for any
OneTree type (otherwise determined from _getroots()).
"""
function _getroot end
function _getroot(tree::AbstractTree{OneTree, <: Rooted})
    @assert nroots(tree) == 1 "More than one root for tree " *
    "(found $(length(roots)))"
    return first(_getroots(tree))
end

"""
    _nroots(::AbstractTree)

Returns the number of roots (subtrees) in a OneTree tree. May be implemented
for any ManyRoots type (otherwise infers from _getroots()).
"""
function _nroots end
_nroots(tree::AbstractTree{OneTree, <: Rooted}) = length(_getroots(tree))
_nroots(tree::AbstractTree{OneTree, Unrooted}) = 0

"""
    _getnode(::AbstractTree, id)
    _getnode(id)

Returns the node associated with id (which could be a name, a node or a pair)
from a tree. For some id types, it will be possible to extract the node
without reference to the tree. Must be implemented for any tree and node label
type.
"""
function _getnode end
_getnode(::AbstractTree{OneTree, RT, NL, N}, node::N) where {RT, NL, N} =
    node
_getnode(::AbstractTree{OneTree, RT, NL, N}, pair::Pair{NL, N}) where
    {RT, NL, N} = pair[2]
_getnode(pair::Pair{NL, N}) where {RT, NL, N <: AbstractNode{RT, NL}} = pair[2]

"""
    _getnodename(::AbstractTree, id)
    _getnodename(id)

Returns the name of a node associated with id (which could be a name,
a node or a pair) from a tree. For some id and node types, it will be
able to extract the node name without reference to the tree.
"""
function _getnodename end
_getnodename(::AbstractTree{OneTree, RT, NL, N, B}, pair::Pair{NL, N}) where
    {RT, NL, N, B} = pair[1]
_getnodename(::AbstractTree{OneTree, RT, NL, N, B}, nodename::NL) where
    {RT, NL, N, B} = nodename
_getnodename(pair::Pair{NL, N}) where {RT, NL, N <: AbstractNode{RT, NL}} =
    pair[1]

"""
    _getbranch(::AbstractTree, id)
    _getbranch(id)

Returns the branch associated with id (which could be a name, a branch or a
pair) from a tree. For some id types, it will be possible to extract the branch
without reference to the tree. Must be implemented for any OneTree tree type.
"""
function _getbranch end
_getbranch(::AbstractTree{OneTree, RT, NL, N, B}, branch::B) where
    {RT, NL, N, B} = branch
_getbranch(::AbstractTree{OneTree, RT, NL, N, B}, pair::Pair{Int, B}) where
    {RT, NL, N, B} = pair[2]
_getbranch(pair::Pair{Int, B}) where {RT, NL, B <: AbstractBranch{RT, NL}} =
    pair[2]

"""
    _getbranchname(::AbstractTree, id)
    _getbranchname(id)

Returns the branch associated with id (which could be a name, a branch or a
pair) from a tree. For some id and branch types, it will be able to extract
the branch name without reference to the tree.
"""
function _getbranchname end
_getbranchname(::AbstractTree{OneTree, RT, NL, N, B}, branchname::Int) where
    {RT, NL, N, B} = branchname
_getbranchname(::AbstractTree{OneTree, RT, NL, N, B}, pair::Pair{Int, B}) where
    {RT, NL, N, B} = pair[1]
_getbranchname(pair::Pair{Int, B}) where
    {RT, NL, B <: AbstractBranch{RT, NL}} = pair[1]

"""
    _hasnode(tree::AbstractTree, node[name])

May be implemented for any OneTree tree type.
"""
function _hasnode end
_hasnode(tree::AbstractTree{OneTree, RT, NL, N, B}, nodename::NL) where
    {RT, NL, N, B} = nodename ∈ _getnodenames(tree)
_hasnode(tree::AbstractTree{OneTree, RT, NL, N, B}, node::N) where
    {RT, NL, N, B} = node ∈ _getnodes(tree)

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
_getbranchnames(tree::AbstractTree{OneTree}) =
    [_getbranchname(tree, branch) for branch in _getbranches(tree)]

"""
    _hasbranch(tree::AbstractTree, branch::AbstractBranch)
    _hasbranch(tree::AbstractTree, source::AbstractNode, dest::AbstractNode)

Tests whether a branch is present in a tree. source and dest method must and
branch method may be implemented for any OneTree tree type.
"""
function _hasbranch end
_hasbranch(::AbstractTree, _) = false
_hasbranch(tree::AbstractTree{OneTree}, branch::AbstractBranch) =
    branch ∈ _getbranches(tree)
_hasbranch(tree::AbstractTree{OneTree}, branchname::Int) =
    branchname ∈ _getbranchnames(tree)

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
function _getrootheight end

"""
    _setrootheight!(::AbstractTree, value)


"""
function _setrootheight! end

"""
    _clearrootheight!(::AbstractTree)


"""
function _clearrootheight! end

"""
    _validate(::AbstractTree)


"""
function _validate(::AbstractTree)
    return true
end

"""
    _traversal(tree::AbstractTree, order::TraversalOrder, todo, sofar)

Return an iterable object containing nodes in given order - preorder, inorder,
postorder or breadthfirst
"""
function _traversal end
function _traversal(tree::AbstractTree{OneTree, <: Rooted},
                    order::TraversalOrder, todo = collect(_getroots(tree)),
                    sofar = eltype(todo)[])
    while !isempty(todo)
        if order == Phylo.breadthfirst
            append!(sofar, todo)
            children = eltype(todo)[]
            for node in todo
                append!(children, getchildren(tree, node))
            end
            todo = children
        else
            now = popfirst!(todo)
            children = getchildren(tree, now)
            if isempty(children)
                push!(sofar, now)
            else
                if order == Phylo.preorder
                    push!(sofar, now)
                    sofar = _traversal(tree, order, children, sofar)
                elseif order == Phylo.postorder
                    sofar = _traversal(tree, order, children, sofar)
                    push!(sofar, now)
                elseif order == Phylo.inorder
                    child = popfirst!(children)
                    sofar = _traversal(tree, order, [child], sofar)
                    push!(sofar, now)
                    sofar = _traversal(tree, order, children, sofar)
                end
            end
        end
    end
    return sofar
end

# AbstractNode methods
"""
    _isleaf(tree::AbstractTree, node)

Is the node a leaf? Does not need to be implemented for any node type –
inferred from _outdegree or _degree - unless tree knows which nodes are
leaves and not nodes.
"""
function _isleaf end
_isleaf(tree::AbstractTree{OneTree, <: Rooted}, node) =
    _outdegree(tree, _getnode(tree, node)) == 0
function _isleaf(tree::AbstractTree{OneTree, Unrooted}, node)
    _degree(tree, node) == 0 && return true
    _degree(tree, node) > 1 && return false
    return missing
end

"""
    _isroot(tree::AbstractTree, node)

Is the node a root node of the tree?
"""
function _isroot end
_isroot(tree::AbstractTree{OneTree, <: Rooted, NL, N}, node::N) where {NL, N} =
    !_hasinbound(tree, node)
_isroot(tree::AbstractTree{OneTree, <: Rooted, NL}, node::NL) where NL =
    _isroot(tree, _getnode(tree, node))
_isroot(tree::AbstractTree{OneTree, Unrooted}, node) = false

"""
    _isinternal(tree::AbstractTree, node::AbstractNode)

Is the node internal to the tree?
"""
function _isinternal end
_isinternal(tree::AbstractTree{OneTree, RT, NL}, node::NL) where {RT, NL} =
    _isinternal(tree, _getnode(tree, node))
_isinternal(tree::AbstractTree, node::AbstractNode{<: Rooted}) =
    _outdegree(tree, node) > 0 && _hasinbound(tree, node)
_isinternal(tree::AbstractTree, node::AbstractNode{Unrooted}) =
    _degree(tree, node) < 2 ? false : missing

"""
    _isunattached(tree::AbstractTree, node::AbstractNode)

Does the node currently form its own (sub)tree?
"""
_isunattached(tree::AbstractTree, node) =
    _degree(tree, _getnode(tree, node)) == 0

"""
    _indegree(tree::AbstractTree, node::AbstractNode)

In degree of node. Can be implemented for rooted nodes, otherwise inferred
from _hasinbound.
"""
function _indegree end
_indegree(tree::AbstractTree, node::AbstractNode) = Int(_hasinbound(tree, node))
_indegree(tree::AbstractTree, node::AbstractNode{Unrooted}) =
    _degree(tree, node) == 0 ? 0 : missing

"""
    _hasinboundspace(tree::AbstractTree, node::AbstractNode)

Is there space for a new inbound connection on a node?
"""
function _hasinboundspace end
_hasinboundspace(tree::AbstractTree, node::AbstractNode) =
    !_hasinbound(tree, node)

"""
    _outdegree(tree::AbstractTree, node::AbstractNode)

Out degree of node.
"""
function _outdegree end
_outdegree(tree::AbstractTree{OneTree, <: Rooted}, node::AbstractNode) =
    length(_getoutbounds(tree, node))
_outdegree(tree::AbstractTree{OneTree, Unrooted},
           node::AbstractNode{Unrooted}) =
    _degree(tree, node) == 0 ? 0 : missing

"""
    _hasoutboundspace(tree::AbstractTree, node::AbstractNode)

Is there space for a new outbound connection on a node? Must be implemented if
a node has a limit on the number of outbound connections (eg for a binary tree)
"""
function _hasoutboundspace end

"""
    _hasspace(tree::AbstractTree, node::AbstractNode)

Is there space for a new connection on a node? Must be implemented if
a node has a limit on the number of connections (eg for a binary tree)
"""
function _hasspace end
_hasspace(::AbstractTree{OneTree}, ::AbstractNode) =
    _hasinboundspace(tree, node) | _hasoutboundspace(tree, node)
_hasspace(::AbstractTree{OneTree, Unrooted}, ::AbstractNode) = true

"""
    _degree(tree::AbstractTree, node::AbstractNode)

Degree of node. Must be implemented for Unrooted nodes, otherwise can be
inferred from indegree and outdegree.
"""
function _degree end
_degree(tree::AbstractTree, node::AbstractNode) =
    _indegree(tree, node) + _outdegree(tree, node)

"""
    _hasinbound(tree::AbstractTree, node::AbstractNode)

Must be implemented for any AbstractNode subtype.
"""
function _hasinbound end
_hasinbound(tree::AbstractTree{OneTree, <: Rooted, NL}, name::NL) where NL =
    _hasinbound(tree, _getnode(tree, name))
_hasinbound(tree::AbstractTree{OneTree, Unrooted},
            node::AbstractNode{Unrooted}) =
    _degree(tree, node) == 0 ? false : missing

"""
    _getinbound(tree::AbstractTree, node::AbstractNode)

Get the inbound connection. Must be implemented for any rooted AbstractNode subtype.
"""
function _getinbound end
_getinbound(::AbstractTree{OneTree, Unrooted}, ::AbstractNode{Unrooted}) =
    error("Can't ask for inbound connection on Unrooted trees")
_getinbound(tree::AbstractTree{OneTree, <: Rooted, NL}, name::NL) where NL =
    _getinbound(tree, _getnode(tree, name))

"""
    _getparent(tree::AbstractTree, node)

Return the parent node for this node. Can be implemented for Rooted node types.
"""
function _getparent end
_getparent(tree::AbstractTree{OneTree, Unrooted}, _) =
    error("Can't ask for parent of Unrooted tree")
_getparent(tree::AbstractTree{OneTree, RT, NL, N, B}, node::N) where
    {RT <: Rooted, NL, N, B} = _getnode(tree,
                                        _src(tree, _getinbound(tree, node)))
_getparent(tree::AbstractTree{OneTree, <: Rooted, NL}, name::NL) where NL =
    _getnodename(tree, _getparent(tree, _getnode(tree, name)))

"""
    _addinbound!(tree::AbstractTree, node::AbstractNode, inbound)

Adds a branch to the input of a rooted node. Must be implemented for any rooted
AbstractNode subtype unless this happens when a branch is created.
"""
function _addinbound! end
_addinbound!(::AbstractTree{OneTree, <: Rooted}, ::AbstractNode{Unrooted}, _) =
    error("Unrooted trees don't have inbound connections")

"""
    _removeinbound!(tree::AbstractTree, node::AbstractNode, inbound)

Removes a branch from the input of a rooted node. Must be implemented for any
rooted AbstractNode subtype unless this happens when a branch is deleted.
"""
function _removeinbound! end
_removeinbound!(::AbstractTree{OneTree, Unrooted}, ::AbstractNode, _) =
    error("Unrooted trees don't have inbound connections")

"""
    _getoutbounds(tree::AbstractTree, node::AbstractNode)

Returns the outbound connections of a rooted node. Must be implemented for any
rooted AbstractNode subtype.
"""
function _getoutbounds end
_getoutbounds(::AbstractTree{OneTree, <: Rooted}, ::AbstractNode{Unrooted}) =
    error("Unrooted trees don't have outbound connections")

"""
    _addoutbound!(tree::AbstractTree, node::AbstractNode, branch)

Add an outbound branch to a rooted node. Must be implemented for any Rooted
AbstractNode subtype unless this happens when a branch is created.
"""
function _addoutbound! end
_addoutbound!(::AbstractTree{OneTree, Unrooted}, ::AbstractNode{Unrooted}, _) =
    error("Unrooted trees don't have inbound and outbound connections")

"""
    _removeoutbound!(tree::AbstractTree, node::AbstractNode, branch)

Remove an outbound branch from a rooted node. Must be implemented for any
AbstractNode subtype unless this happens when a branch is deleted.
"""
function _removeoutbound! end
_removeoutbound!(::AbstractTree{OneTree, Unrooted},
                 ::AbstractNode{Unrooted}, _) =
    error("Unrooted trees don't have outbound connections")

"""
    _getchildren(tree::AbstractTree, node)
    _getchildren(tree::AbstractTree, nodename)

Return the child node(s) for this node. May be implemented for any rooted
AbstractNode subtype.
"""
function _getchildren end
_getchildren(::AbstractTree{OneTree, Unrooted}, _) =
    error("Can't ask for children of Unrooted tree")
_getchildren(tree::AbstractTree{OneTree, <: Rooted, NL}, name::NL) where NL =
    [_getnodename(tree, _dst(tree, branch))
     for branch in _getoutbounds(tree, _getnode(tree, name))]

_getchildren(tree::AbstractTree{OneTree, RT, NL, N, B}, node::N) where
    {RT <: Rooted, NL, N, B} =
    [_getnode(tree, _dst(tree, branch)) for branch in _getoutbounds(tree, node)]

"""
    _getconnections(tree::AbstractTree, node::AbstractNode)

Returns all of the connections of a node. Must be implemented for any unrooted
AbstractNode subtype, can be inferred from _getinbound and _getoutbounds for
a rooted node.
"""
function _getconnections end
_getconnections(tree::AbstractTree{OneTree}, node::AbstractNode) =
    _hasinbound(tree, node) ?
        append!([_getinbound(tree, node)], _getoutbounds(tree, node)) :
        _getoutbounds(tree, node)
_getconnections(tree::AbstractTree{OneTree, <: Rooted, NL}, name::NL) where NL =
    _hasinbound(tree, name) ?
        append!([_getinbound(tree, name)], _getoutbounds(tree, name)) :
        _getoutbounds(tree, name)

"""
    _getsiblings(tree::AbstractTree, node::AbstractNode)

Returns all of the siblings of a node. May be implemented for any
AbstractNode subtype, can be inferred from _getparent and _getchildren for
a rooted node or _getconnections for an unrooted node.
"""
function _getsiblings end
_getsiblings(tree::AbstractTree{OneTree, <: Rooted}, node::AbstractNode) =
    _hasinbound(tree, node) ?
        append!([_getparent(tree, node)], _getchildren(tree, node)) :
        _getchildren(tree, node)
_getsiblings(tree::AbstractTree{OneTree, Unrooted}, node::AbstractNode) =
    [_conn(tree, branch, node) for branch in _getconnections(tree, node)]
_getsiblings(tree::AbstractTree{OneTree, <: Rooted, NL}, name::NL) where NL =
    _hasinbound(tree, name) ?
        append!([_getparent(tree, name)], _getchildren(tree, name)) :
        _getchildren(tree, name)
_getsiblings(tree::AbstractTree{OneTree, Unrooted, NL}, name::NL) where NL =
    [_getnodename(tree, _conn(tree, branch, _getnode(tree, name)))
     for branch in _getconnections(tree, name)]

"""
    _addconnection!(tree::AbstractTree, node::AbstractNode, branch)

Add a connection to an unrooted node. Must be implemented for any unrooted
AbstractNode subtype unless this happens when a branch is added.
"""
function _addconnection! end
_addconnection!(::AbstractTree, ::AbstractNode{<: Rooted}, _) =
    error("Direction of connection must be specified for a rooted node")

"""
    _removeconnection!(tree::AbstractTree, node::AbstractNode, branch)

Remove a connection from an unrooted node. Must be implemented for any Unrooted
AbstractNode subtype unless this happens when a branch is deleted.
"""
function _removeconnection! end
_removeconnection!(::AbstractTree, ::AbstractNode{<: Rooted}, _) =
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
_src(::AbstractTree, ::AbstractBranch{Unrooted}) =
    error("Unrooted branches do not have in and out connections")

"""
    _dst(branch::AbstractBranch)

Return destination node for a branch. Must be implemented for any rooted
AbstractBranch subtype.
"""
function _dst end
_dst(::AbstractTree, ::AbstractBranch{Unrooted}) =
    error("Unrooted branches do not have in and out connections")

"""
    _conns(tree::AbstractTree, branch::AbstractBranch)

Return a vector of connections for a branch. Must be implemented for any
Unrooted AbstractBranch subtype, otherwise can combine _src and _dst.
"""
function _conns end
_conns(tree::AbstractTree{OneTree, <: Rooted}, branch::AbstractBranch) =
    [_src(tree, branch), _dst(tree, branch)]

"""
    _conn(branch::AbstractBranch, exclude::AbstractNode)

Return the connection for a branch that isn't the `exclude` node. May be
implemented for any Unrooted AbstractBranch subtype, otherwise will use
_conns.
"""
function _conn end
function _conn(tree::AbstractTree{OneTree, RT, NL, N, B},
               branch::B, exclude::N) where {RT <: Rooted, NL,
                                             B <: AbstractBranch{RT, NL},
                                             N <: AbstractNode{RT, NL}}
    other = [node for node in _conns(tree, branch) if node != exclude]
    @assert length(other) == 1 _getnodename(tree, exclude) *
        " is not connected to branch $(_getbranchname(tree, branch))"
    return first(other)
end

"""
    _getlength

Return length of a branch. May be implemented for any AbstractBranch
subtype.
"""
function _getlength end
_getlength(::AbstractTree, ::AbstractBranch) = NaN
_getlength(tree::AbstractTree, bn::Int) = _getlength(tree, _getbranch(tree, bn))

"""
    _leafinfotype(::Type{<:AbstractTree})

Returns the type of the leaf info data.
"""
function _leafinfotype end
_leafinfotype(::Type{<:AbstractTree}) = Nothing

"""
    _noderecordtype(::Type{<:AbstractTree})

Returns the type of the node info data.
"""
function _noderecordtype end
_noderecordtype(::Type{<:AbstractTree}) = Nothing

"""
    _branchrecordtype(::Type{<:AbstractTree})

Returns the type of the branch info data.
"""
function _branchrecordtype end
_branchrecordtype(::Type{<:AbstractTree}) = Nothing

function _getleafinfo end
function _setleafinfo! end
function _getnoderecord end
function _setnoderecord! end
function _getbranchinfo end
function _setbranchinfo! end
function _resetleaves! end
