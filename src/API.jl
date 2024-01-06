using Phylo
using Phylo: Rootedness, Rooted, TreeType, TraversalOrder
using Phylo: AbstractElt, AbstractNode, AbstractBranch, AbstractTree
using SimpleTraits
using Unitful

#! format: off
@traitdef PreferNodeObjects{X}
@traitimpl PreferNodeObjects{X} <- _prefernodeobjects(X)

@traitdef PreferBranchObjects{X}
@traitimpl PreferBranchObjects{X} <- _preferbranchobjects(X)

@traitdef MatchNodeType{T, N}
@traitimpl MatchNodeType{T, N} <- _matchnodetype(T, N)

@traitdef MatchNodeTypes{T, N1, N2}
@traitimpl MatchNodeTypes{T, N1, N2} <- _matchnodetypes(T, N1, N2)

@traitdef HoldsNodeData{T, ND}
@traitimpl HoldsNodeData{T, ND} <- (_nodedatatype(T) ≡ ND)

@traitdef MatchBranchType{T, B}
@traitimpl MatchBranchType{T, B} <- _matchbranchtype(T, B)

@traitdef MatchBranchNodeType{T, B, N}
@traitimpl MatchBranchNodeType{T, B, N} <- _matchbranchnodetype(T, B, N)

@traitdef MatchTreeNameType{T, TN}
@traitimpl MatchTreeNameType{T, TN} <- _matchtreenametype(T, TN)
#! format: on

"""
    _prefernodeobjects(::Type{<:AbstractTree})

Does this tree or node type prefer nodes to be objects or names? Must be
implemented for every node type.
"""
function _prefernodeobjects end
function _prefernodeobjects(::Type{<:AbstractTree{TT, RT, NL, N, B}}) where
         {TT, RT, NL, N, B}
    return _prefernodeobjects(N)
end

"""
    _preferbranchobjects(::Type{<:AbstractTree})

Does this tree or branch type prefer branches to be objects or names? Must be
implemented for every branch type.
"""
function _preferbranchobjects end
function _preferbranchobjects(::Type{<:AbstractTree{TT, RT, NL, N, B}}) where
         {TT, RT, NL, N, B}
    return _preferbranchobjects(B)
end

"""
    _matchnodetype(::Type{<:AbstractTree{TT, RT, NL, N, B}}, ::Type{N})
    _matchnodetype(::Type{<:AbstractTree{TT, RT, NL, N, B}}, ::Type{NL})

Does this tree type prefer the node or node label type provided?
"""
function _matchnodetype end
function _matchnodetype(::Type{<:AbstractTree{TT, RT, NL, N, B}},
                        ::Type{N}) where {TT, RT, NL, N, B}
    return _prefernodeobjects(N)
end
function _matchnodetype(::Type{<:AbstractTree{TT, RT, NL, N, B}},
                        ::Type{NL}) where {TT, RT, NL, N, B}
    return !_prefernodeobjects(N)
end
_matchnodetypes(::Type{<:AbstractTree}, ::Type{<:Any}, ::Type{<:Any}) = false
function _matchnodetypes(::Type{T}, ::Type{N},
                         ::Type{N}) where {T <: AbstractTree, N}
    return _matchnodetype(T, N)
end

"""
    _matchbranchtype(::Type{<:AbstractTree}, ::Type{<:AbstractBranch})

Does this tree type prefer the branch or branch label type provided?
"""
function _matchbranchtype end
function _matchbranchtype(::Type{<:AbstractTree{TT, RT, NL, N, B}},
                          ::Type{B}) where {TT, RT, NL, N, B}
    return _preferbranchobjects(B)
end
function _matchbranchtype(::Type{<:AbstractTree{TT, RT, NL, N, B}},
                          ::Type{Int}) where {TT, RT, NL, N, B}
    return !_preferbranchobjects(B)
end

"""
    _matchbranchnodetype(::Type{<:AbstractTree},
                         ::Type{<:AbstractBranch},
                         ::Type{<:AbstractNode})

Does this tree type prefer the branch and node types provided?
"""
function _matchbranchnodetype end
function _matchbranchnodetype(::Type{T}, ::Type{B},
                              ::Type{N}) where {T <: AbstractTree, B, N}
    return _matchbranchtype(T, B) && _matchnodetype(T, N)
end

"""
    _matchtreenametype(::Type{<:AbstractTree}, ::Type{X})

Does this tree type prefer the node or node label type provided?
"""
_matchtreenametype(::Type{T}, ::Type{TN}) where {T <: AbstractTree, TN} = (TN ≡
                                                                           _treenametype(T))

# Need to be able to create new node and branch labels
function _newlabel(ids::AbstractVector{Label}) where {Label <: Integer}
    return isempty(ids) ? 1 : maximum(ids) + 1
end
function _newlabel(ids::AbstractVector{<:AbstractVector{Label}}) where {
                                                                        Label <:
                                                                        Integer}
    return mapreduce(_newlabel, max, 1, ids)
end
function _newlabel(names::AbstractVector{String}, prefix::String)
    return prefix * "$(_newlabelfree(names, prefix))"
end
function _newlabel(names::AbstractVector{<:AbstractVector{String}},
                   prefix::String)
    return prefix * "$(mapreduce(_newlabelfree, max, 1, names))"
end

function _newlabelfree(names::AbstractVector{String}, prefix)
    num = length(names) + 1
    if "$prefix$num" ∉ names
        return num
    end
    len = length(prefix)
    goodnames = [name[(len - 1):end]
                 for name in names
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
function _newnodelabel(tree::AbstractTree{TT, RT, String, N, B}) where
         {TT, RT, N, B}
    return _newlabel(collect(getnodenames(tree)), "Node ")
end
function _newnodelabel(tree::AbstractTree{TT, RT, I, N, B}) where
         {TT, RT, I <: Integer, N, B}
    return _newlabel(collect(getnodenames(tree)))
end

# AbstractTree type methods
"""
    _treenametype(::Type{AbstractTree})

Returns the label type for a tree type. Must be implemented for any tree type.
"""
function _treenametype end
_treenametype(::Type{<:AbstractTree{OneTree}}) = Int

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
_gettreenames(tree::AbstractTree) = _gettreename.(_gettrees(tree))

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
function _gettree(tree::T, id) where {T <: AbstractTree{OneTree}}
    id == _gettreename(tree) || throw(BoundsError(tree, id))
    return tree
end

"""
    _getnodes(tree::AbstractTree{OneTree}[, order::TraversalOrder])

Returns an interable collection of nodes for a OneTree tree. _getnodes(tree) must be
implemented for a OneTree tree type as a base mechanisms for extracting the
node list.
"""
function _getnodes end
function _getnodes(tree::AbstractTree{OneTree}, order::TraversalOrder)
    return _traversal(tree, order)
end

"""
    _getnodenames(tree::AbstractTree{OneTree})

Returns an iterable collection of node names for a OneTree tree. Can
be implemented for any OneTree tree type, especially PreferNodeObjects trees.
"""
function _getnodenames end
function _getnodenames(::T) where {T <: AbstractTree}
    return error("No _getnodes() method for tree type $T")
end
@traitfn function _getnodenames(tree::T,
                                order::TraversalOrder) where
                  {T <: AbstractTree{OneTree};
                   !PreferNodeObjects{T}}
    return _getnodes(tree, order)
end
@traitfn function _getnodenames(tree::T,
                                order::TraversalOrder) where
                  {T <: AbstractTree{OneTree}; PreferNodeObjects{T}}
    return _getnodename.(Ref(tree), _getnodes(tree, order))
end

"""
    _nnodes(::AbstractTree)

Returns the number of nodes (internal nodes and leaves) in a single tree. May
be implemented for any OneTree tree type (otherwise infers from _getnodes()).
"""
_nnodes(tree::AbstractTree{OneTree}) = length(_getnodes(tree, anyorder))

"""
    _getleafnames(::AbstractTree, ::TraversalOrder)

Returns the leaf names of a tree. May be implemented for any
tree type (otherwise determined from _getnodenames() and _isleaf() functions).
"""
function _getleafnames end
function _getleafnames(tree::AbstractTree{OneTree}, order::TraversalOrder)
    return [_getnodename(tree, node)
            for node in _traversal(tree, order) if _isleaf(tree, node)]
end
function _getleafnames(tree::AbstractTree{ManyTrees}, order::TraversalOrder)
    return _getleafnames(_gettree(tree, first(_gettreenames(tree))), order)
end

"""
    _getleaves(::AbstractTree)

Returns the leaves (tips) of a single tree. May be implemented for any
OneTree type (otherwise determined from _getnodes() and _isleaf() functions).
"""
function _getleaves end
function _getleaves(tree::AbstractTree{OneTree}, order::TraversalOrder)
    return [node for node in _traversal(tree, order) if _isleaf(tree, node)]
end

"""
    _nleaves(::AbstractTree)

Returns the number of leaves (tips) in a tree. May be implemented for any
tree type (otherwise determined from the _getleafnames() function).
"""
_nleaves(tree::AbstractTree) = length(_getleafnames(tree, anyorder))

"""
    _getroots(::AbstractTree)

Returns the root(s) of a tree. May be implemented for any
OneTree type (otherwise determined from _getnodes() and _isroot() functions).
"""
_getroots(tree::AbstractTree{OneTree, <:Rooted}) = [node
                                                    for node in _getnodes(tree)
                                                    if _isroot(tree, node)]

"""
    _getroot(::AbstractTree)

Returns the unique root of a rooted tree. May be implemented for any
OneTree type (otherwise determined from _getroots()).
"""
function _getroot end
function _getroot(tree::AbstractTree{OneTree, <:Rooted})
    @assert _nroots(tree)==1 "More than one root for tree "*
    "(found $(_nroots(tree)))"
    return first(_getroots(tree))
end

"""
    _nroots(::AbstractTree)

Returns the number of roots (subtrees) in a OneTree tree. May be implemented
for any ManyRoots type (otherwise infers from _getroots()).
"""
function _nroots end
_nroots(tree::AbstractTree{OneTree, <:Rooted}) = length(_getroots(tree))
_nroots(tree::AbstractTree{OneTree, Unrooted}) = 0

"""
    _getnode(::AbstractTree, id)

Returns the node or name associated with id (which could be a name or
a node) from a tree. Must be implemented for any PreferNodeObjects
tree and node label type.
"""
function _getnode end
@traitfn function _getnode(tree::T,
                           node::N) where {RT, NL, N,
                                           T <:
                                           AbstractTree{OneTree, RT, NL, N};
                                           PreferNodeObjects{T}}
    return node
end
@traitfn function _getnode(tree::T,
                           nodename::NL) where {RT, NL,
                                                T <:
                                                AbstractTree{OneTree, RT, NL};
                                                !PreferNodeObjects{T}}
    return nodename
end

"""
    _getnodename(::AbstractTree, id)

Returns the name of a node associated with id (which could be a name or a node)
from a tree. Must be implemented for PreferNodeObjects tree types.
"""
function _getnodename end
function _getnodename(::AbstractTree{OneTree, RT, NL, N},
                      nodename::NL) where {RT, NL, N}
    return nodename
end

"""
    _hasnode(tree::AbstractTree, node[name])

Does the tree contain this node? Must be implemented for any PreferNodeObjects
tree type with a node label.
"""
function _hasnode end
function _hasnode(tree::AbstractTree{OneTree, RT, NL},
                  node::NL) where {RT, NL}
    return node ∈ _getnodenames(tree)
end
function _hasnode(tree::AbstractTree{OneTree, RT, NL, N},
                  node::N) where {RT, NL, N}
    return node ∈ _getnodes(tree)
end

"""
    _renamenode!(tree::AbstractTree, oldnode[name], newname)

Renames a node in a tree. Optional - not implemented for most tree types.
"""
function _renamenode! end
_renamenode!(_::AbstractTree, _, _) = false

"""
    _createnode!(tree::AbstractTree, nodename[, data])

Must be implemented for any AbstractTree subtype.
"""
function _createnode! end

"""
    _deletenode!(tree::AbstractTree, nodename)

Must be implemented for any AbstractTree subtype.
"""
function _deletenode! end

"""
    _getbranch(::AbstractTree, id)

Returns the branch or name associated with id (which could be a name or
a branch) from a tree. Must be implemented for any PreferBranchObjects
tree and branch label type.
"""
function _getbranch end
@traitfn function _getbranch(::T,
                             branch::B) where {RT, NL, N, B,
                                               T <:
                                               AbstractTree{OneTree, RT, NL, N,
                                                            B};
                                               PreferBranchObjects{T}}
    return branch
end
@traitfn function _getbranch(::T,
                             branchname::Int) where {T <: AbstractTree{OneTree};
                                                     !PreferBranchObjects{T}}
    return branchname
end
@traitfn function _getbranch(tree::T, src::N1,
                             dst::N2) where
                  {T <: AbstractTree{OneTree}, N1, N2; MatchNodeTypes{T, N1,
                                                                      N2}}
    return first(b
                 for b in _getbranches(tree)
                 if _src(tree, b) == src &&
                    _dst(tree, b) == dst)
end

@traitfn function _getbranch(tree::T, src::N1,
                             dst::N2) where
                  {T <: AbstractTree{OneTree}, N1,
                   N2; !MatchNodeTypes{T, N1, N2}}
    return _getbranch(tree,
                      _getnode(tree, src),
                      _getnode(tree, dst))
end

"""
    _getbranchname(::AbstractTree, id)

Returns the name of a branch associated with id (which could be a name or a branch)
from a tree. Must be implemented for PreferBranchObjects tree types.
"""
function _getbranchname end
_getbranchname(::AbstractTree{OneTree}, id::Int) = id

"""
    _getbranches(tree::AbstractTree)

Returns a vector of branches for a OneTree tree. Either _getbranches() or
_getbranchnames() must be implemented for any OneTree tree type.
"""
function _getbranches end
_getbranches(::T) where {T} = error("No _getbranches() method for tree type $T")

"""
    _nbranches(::AbstractTree)

Returns the number of branches in a single tree. May be implemented for any
OneTree tree type (otherwise infers from _getbranches()).
"""
_nbranches(tree::AbstractTree{OneTree}) = length(_getbranches(tree))

"""
    _getbranchnames(tree::AbstractTree)

Returns a vector of branch names for a OneTree tree. Either _getbranches() or
_getbranchnames() must be implemented for any OneTree tree type.
"""
function _getbranchnames end
@traitfn function _getbranchnames(tree::T) where {T <: AbstractTree{OneTree};
                                                  !PreferBranchObjects{T}}
    return _getbranches(tree)
end
@traitfn function _getbranchnames(tree::T) where {T <: AbstractTree{OneTree};
                                                  PreferBranchObjects{T}}
    return _getbranchname.(Ref(tree), _getbranches(tree))
end

"""
    _hasbranch(tree::AbstractTree, node[name])

Does the tree contain this branch? Must be implemented for any
PreferBranchObjects tree type with a branch label.
"""
function _hasbranch end
function _hasbranch(tree::AbstractTree{OneTree},
                    branch::Int)
    return branch ∈ _getbranchnames(tree)
end
@traitfn function _hasbranch(tree::T,
                             branch::B) where {RT, NL, N, B,
                                               T <:
                                               AbstractTree{OneTree, RT, NL, N,
                                                            B};
                                               PreferBranchObjects{T}}
    return branch ∈ _getbranches(tree)
end

@traitfn function _hasbranch(tree::T, src::N1,
                             dst::N2) where
                  {T <: AbstractTree{OneTree}, N1, N2; MatchNodeTypes{T, N1,
                                                                      N2}}
    return any(_src(tree, b) == src && _dst(tree, b) == dst
               for b in _getbranches(tree))
end

"""
    _createbranch!(tree::AbstractTree, source, destination[,
                   length][, data])

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
    _validate!(::AbstractTree)

Check whether the tree is internally valid.
"""
function _validate! end
_validate!(::AbstractTree) = true

"""
    _invalidate!(::AbstractTree, state)

Confirm that the tree is no longer necessarily valid, and remove cache information.
"""
function _invalidate! end

"""
    _traversal(tree::AbstractTree, order::TraversalOrder, todo, sofar)

Return an iterable object containing nodes in given order - preorder, inorder,
postorder or breadthfirst
"""
function _traversal end
function _traversal(tree::T,
                    order::TraversalOrder) where {
                                                  T <: AbstractTree{OneTree,
                                                               <:Rooted}}
    if order == anyorder
        return _getnodes(tree)
    else
        return _traversal(tree, order, collect(_getroots(tree)))
    end
end
function _traversal(tree::AbstractTree{OneTree, <:Rooted},
                    order::TraversalOrder, todo,
                    sofar::AbstractVector = eltype(todo)[])
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
function _isleaf(tree::AbstractTree{OneTree, <:Rooted}, node)
    return _outdegree(tree, _getnode(tree, node)) == 0
end
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
_isroot(tree::AbstractTree{OneTree, <:Rooted}, node) = !_hasinbound(tree, node)
_isroot(tree::AbstractTree{OneTree, Unrooted}, node) = false

"""
    _isinternal(tree::AbstractTree, node)

Is the node internal to the tree?
"""
function _isinternal end
function _isinternal(tree::AbstractTree{OneTree, <:Rooted}, node)
    return _outdegree(tree, node) > 0 && _hasinbound(tree, node)
end
function _isinternal(tree::AbstractTree{OneTree, Unrooted}, node)
    return _degree(tree, node) < 2 ? false : missing
end

"""
    _isunattached(tree::AbstractTree, node)

Does the node currently form its own (sub)tree?
"""
_isunattached(tree::AbstractTree, node) = _degree(tree, node) == 0

"""
    _indegree(tree::AbstractTree, node)

In degree of node. Can be implemented for rooted nodes, otherwise inferred
from _hasinbound.
"""
function _indegree end
function _indegree(tree::AbstractTree{OneTree, <:Rooted}, node)
    return Int(_hasinbound(tree, node))
end
function _indegree(tree::AbstractTree{OneTree, Unrooted}, node)
    return _degree(tree, node) == 0 ? 0 : missing
end

"""
    _hasinboundspace(tree::AbstractTree, node::AbstractNode)

Is there space for a new inbound connection on a node?
"""
function _hasinboundspace end
function _hasinboundspace(tree::AbstractTree{OneTree, <:Rooted}, node)
    return !_hasinbound(tree, node)
end

"""
    _outdegree(tree::AbstractTree, node::AbstractNode)

Out degree of node.
"""
function _outdegree end
function _outdegree(tree::AbstractTree{OneTree, <:Rooted}, node)
    return length(_getoutbounds(tree, node))
end
function _outdegree(tree::AbstractTree{OneTree, Unrooted},
                    node::AbstractElt{Unrooted})
    return _degree(tree, node) == 0 ? 0 : missing
end

"""
    _hasoutboundspace(tree::AbstractTree, node::AbstractNode)

Is there space for a new outbound connection on a node? Must be implemented if
a node has a limit on the number of outbound connections (eg for a binary tree)
"""
function _hasoutboundspace end
_hasoutboundspace(::AbstractTree{OneTree, <:Rooted}, _) = true

"""
    _hasspace(tree::AbstractTree, node::AbstractNode)

Is there space for a new connection on a node? Must be implemented if
a node has a limit on the number of connections (eg for a binary tree)
"""
function _hasspace end
function _hasspace(tree::AbstractTree{OneTree, <:Rooted}, node)
    return _hasinboundspace(tree, node) | _hasoutboundspace(tree, node)
end
_hasspace(::AbstractTree{OneTree, Unrooted}, _) = true

"""
    _degree(tree::AbstractTree, node::AbstractNode)

Degree of node. Must be implemented for Unrooted nodes, otherwise can be
inferred from indegree and outdegree.
"""
function _degree end
function _degree(tree::AbstractTree{OneTree, <:Rooted}, node)
    return _indegree(tree, node) + _outdegree(tree, node)
end

"""
    _hasinbound(tree::AbstractTree, node::AbstractNode)

Must be implemented for any AbstractNode subtype.
"""
function _hasinbound end
function _hasinbound(::T, ::N) where {T, N}
    return error("No _hasinbound() function for types $T, $N")
end

"""
    _getinbound(tree::AbstractTree, node::AbstractNode)

Get the inbound connection. Must be implemented for any rooted AbstractNode subtype.
"""
function _getinbound end
function _getinbound(::T, ::B) where {T, B}
    return error("No _getinbound() function for $T, $B")
end

"""
    _getparent(tree::AbstractTree, node)

Return the parent node for this node. Can be implemented for Rooted node types.
"""
function _getparent end
function _getparent(tree::AbstractTree{OneTree}, node)
    return _getnode(tree, _src(tree, _getinbound(tree, node)))
end

"""
    _addinbound!(tree::AbstractTree, node::AbstractNode, inbound)

Adds a branch to the input of a rooted node. Must be implemented for any rooted
AbstractNode subtype unless this happens when a branch is created.
"""
function _addinbound! end

"""
    _removeinbound!(tree::AbstractTree, node::AbstractNode, inbound)

Removes a branch from the input of a rooted node. Must be implemented for any
rooted AbstractNode subtype unless this happens when a branch is deleted.
"""
function _removeinbound! end

"""
    _getoutbounds(tree::AbstractTree, node::AbstractNode)

Returns the outbound connections of a rooted node. Must be implemented for any
rooted AbstractNode subtype.
"""
function _getoutbounds end
function _getoutbounds(::T, ::N) where {T, N}
    return error("No _getoutbounds() function for $T, $N")
end

"""
    _addoutbound!(tree::AbstractTree, node::AbstractNode, branch)

Add an outbound branch to a rooted node. Must be implemented for any Rooted
AbstractNode subtype unless this happens when a branch is created.
"""
function _addoutbound! end

"""
    _removeoutbound!(tree::AbstractTree, node::AbstractNode, branch)

Remove an outbound branch from a rooted node. Must be implemented for any
AbstractNode subtype unless this happens when a branch is deleted.
"""
function _removeoutbound! end

"""
    _getchildren(tree::AbstractTree, node)
    _getchildren(tree::AbstractTree, nodename)

Return the child node(s) for this node. May be implemented for any rooted
AbstractNode subtype.
"""
function _getchildren end
function _getchildren(tree::AbstractTree{OneTree, RT},
                      node::N) where {RT <: Rooted, N <: AbstractElt{RT}}
    return N[_dst(tree, branch) for branch in _getoutbounds(tree, node)]
end
function _getchildren(tree::AbstractTree{OneTree, RT, NL},
                      node::NL) where {RT <: Rooted, NL}
    return NL[_dst(tree, branch) for branch in _getoutbounds(tree, node)]
end

"""
    _getconnections(tree::AbstractTree, node::AbstractNode)

Returns all of the connections of a node. Must be implemented for any unrooted
AbstractNode subtype, can be inferred from _getinbound and _getoutbounds for
a rooted node.
"""
function _getconnections end
function _getconnections(tree::AbstractTree{OneTree, <:Rooted}, node, exclude)
    return filter(∉(exclude),
                  _hasinbound(tree, node) ?
                  append!([_getinbound(tree, node)],
                          _getoutbounds(tree, node)) :
                  _getoutbounds(tree, node))
end

"""
    _getsiblings(tree::AbstractTree, node::AbstractNode)

Returns all of the siblings (actually immediate connections) of a node.
May be implemented for any AbstractNode subtype, can be inferred from
_getparent and _getchildren for a rooted node or _getconnections for an
unrooted node.
"""
function _getsiblings end
function _getsiblings(tree::AbstractTree{OneTree, <:Rooted}, node)
    return _hasinbound(tree, node) ?
           append!([_getparent(tree, node)], _getchildren(tree, node)) :
           _getchildren(tree, node)
end
function _getsiblings(tree::AbstractTree{OneTree, Unrooted}, node)
    return [_conn(tree, branch, node)
            for branch in _getconnections(tree, node, [])]
end

"""
    _addconnection!(tree::AbstractTree, node::AbstractNode, branch)

Add a connection to an unrooted node. Must be implemented for any unrooted
AbstractNode subtype unless this happens when a branch is added.
"""
function _addconnection! end
function _addconnection!(::AbstractTree{OneTree, <:Rooted}, _, _)
    return error("Direction of connection must be specified for a rooted node")
end

"""
    _removeconnection!(tree::AbstractTree, node::AbstractNode, branch)

Remove a connection from an unrooted node. Must be implemented for any Unrooted
AbstractNode subtype unless this happens when a branch is deleted.
"""
function _removeconnection! end
function _removeconnection!(::AbstractTree{OneTree, <:Rooted}, _, _)
    return error("Direction of connection must be specified for a rooted node")
end

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
    return missing
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
    _src(tree, branch)

Return source node for a branch. Must be implemented for any rooted
branch type.
"""
function _src end
_src(::T, ::B) where {T, B} = error("No _src() function for $T, $B")

"""
    _dst(branch::AbstractBranch)

Return destination node for a branch. Must be implemented for any rooted
AbstractBranch subtype.
"""
function _dst end
_dst(::T, ::B) where {T, B} = error("No _dst() function for $T, $B")

"""
    _conns(tree::AbstractTree, branch::AbstractBranch)

Return a vector of connections for a branch. Must be implemented for any
Unrooted AbstractBranch subtype, otherwise can combine _src and _dst.
"""
function _conns end
function _conns(tree::AbstractTree{OneTree, <:Rooted}, branch)
    return [_src(tree, branch), _dst(tree, branch)]
end

"""
    _conn(branch::AbstractBranch, exclude::AbstractNode)

Return the connection for a branch that isn't the `exclude` node. May be
implemented for any Unrooted AbstractBranch subtype, otherwise will use
_conns.
"""
function _conn end
function _conn(tree::AbstractTree{OneTree, <:Rooted}, branch, exclude)
    other = [node for node in _conns(tree, branch) if node !== exclude]
    @assert length(other)==1 _getnodename(tree,
                                          exclude)*
    " is not connected to branch $(_getbranchname(tree, branch))"
    return first(other)
end

"""
    _getlength

Return length of a branch. May be implemented for any AbstractBranch
subtype.
"""
function _getlength end
_getlength(::AbstractTree, _) = missing

"""
    _haslength

Return length of a branch. May be implemented for any AbstractBranch
subtype.
"""
function _haslength end
_haslength(t::AbstractTree, b) = !ismissing(_getlength(t, b))

"""
    _leafinfotype(::Type{<:AbstractTree})

Returns the type of the leaf info data.
"""
function _leafinfotype end
_leafinfotype(::Type{<:AbstractTree}) = Nothing

"""
    _nodedatatype(::Type{<:AbstractTree})

Returns the type of the node info data.
"""
function _nodedatatype end
_nodedatatype(::Type{<:AbstractTree}) = Nothing

"""
    _branchdims(::Type{<:AbstractTree})

Returns the dimensions of the branch lengths for the tree.
"""
function _branchdims end
_branchdims(::Type{<:AbstractTree}) = NoDims

"""
    _branchdatatype(::Type{<:AbstractTree})

Returns the type of the branch info data.
"""
function _branchdatatype end
_branchdatatype(::Type{<:AbstractTree}) = Nothing

function _getnodedata end
function _getnodedata(::T, ::N) where {T, N}
    return error("No _getnodedata() function for $T, $N")
end
function _getnodedata(tree::AbstractTree{OneTree}, node, label)
    return _getnodedata(tree, node)[label]
end
function _setnodedata! end
function _setnodedata!(tree::AbstractTree{OneTree}, node, label, value)
    return (_getnodedata(tree, node)[label] = value)
end

function _getbranchdata end
function _getbranchdata(::T, ::B) where {T, B}
    return error("No _getbranchdata() function for $T, $B")
end
function _getbranchdata(tree::AbstractTree, branch, label)
    return _getbranchdata(tree, branch)[label]
end
function _setbranchdata! end
function _setbranchdata!(tree::AbstractTree, branch, label, value)
    return (_getbranchdata(tree, branch)[label] = value)
end

function _getleafinfo end
function _setleafinfo! end

"""
    _gettreeinfo(tree::AbstractTree)
    _gettreeinfo(tree::AbstractTree, treename)

Returns the info data associated with the tree(s).
"""
function _gettreeinfo end
function _gettreeinfo(tree::AbstractTree{OneTree})
    return Dict(_gettreename(tree) => Dict{String, Any}())
end
function _gettreeinfo(tree::AbstractTree{OneTree}, treename)
    @assert _gettreename(tree)==treename "No tree called $treename"
    return Dict{String, Any}()
end

"""
    _resetleaves!(::AbstractTree)

Fixes leaf naming after creation or deletion of nodes or branches. Must be
implemented by tree types where this is handled separately.
"""
function _resetleaves! end
_resetleaves!(::AbstractTree) = true
