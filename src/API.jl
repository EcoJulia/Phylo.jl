using Phylo
using Phylo: Rootedness, Rooted, TreeType, TraversalOrder
using Phylo: AbstractElt, AbstractNode, AbstractBranch, AbstractTree
using SimpleTraits
using Unitful

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

"""
    _prefernodeobjects(::Type{<:AbstractTree})

Does this tree or node type prefer nodes to be objects or names? Must be
implemented for every node type.
"""
function _prefernodeobjects end
_prefernodeobjects(::Type{<:AbstractTree{TT, RT, NL, N, B}}) where
{TT, RT, NL, N, B} = _prefernodeobjects(N)

"""
    _preferbranchobjects(::Type{<:AbstractTree})

Does this tree or branch type prefer branches to be objects or names? Must be
implemented for every branch type.
"""
function _preferbranchobjects end
_preferbranchobjects(::Type{<:AbstractTree{TT, RT, NL, N, B}}) where
{TT, RT, NL, N, B} = _preferbranchobjects(B)

"""
    _matchnodetype(::Type{<:AbstractTree{TT, RT, NL, N, B}}, ::Type{N})
    _matchnodetype(::Type{<:AbstractTree{TT, RT, NL, N, B}}, ::Type{NL})

Does this tree type prefer the node or node label type provided?
"""
function _matchnodetype end
_matchnodetype(::Type{<:AbstractTree{TT, RT, NL, N, B}},
               ::Type{N}) where {TT, RT, NL, N, B} = _prefernodeobjects(N)
_matchnodetype(::Type{<:AbstractTree{TT, RT, NL, N, B}},
               ::Type{NL}) where {TT, RT, NL, N, B} = !_prefernodeobjects(N)
_matchnodetypes(::Type{<:AbstractTree}, ::Type{<:Any}, ::Type{<:Any}) = false
_matchnodetypes(::Type{T}, ::Type{N}, ::Type{N}) where {T <: AbstractTree, N} =
    _matchnodetype(T, N)

    """
    _matchbranchtype(::Type{<:AbstractTree}, ::Type{<:AbstractBranch})

Does this tree type prefer the branch or branch label type provided?
"""
function _matchbranchtype end
_matchbranchtype(::Type{<:AbstractTree{TT, RT, NL, N, B}},
                 ::Type{B}) where {TT, RT, NL, N, B} =
                     _preferbranchobjects(B)
_matchbranchtype(::Type{<:AbstractTree{TT, RT, NL, N, B}},
                 ::Type{Int}) where {TT, RT, NL, N, B} =
                     !_preferbranchobjects(B)

"""
    _matchbranchnodetype(::Type{<:AbstractTree},
                         ::Type{<:AbstractBranch},
                         ::Type{<:AbstractNode})

Does this tree type prefer the branch and node types provided?
"""
function _matchbranchnodetype end
_matchbranchnodetype(::Type{T}, ::Type{B}, ::Type{N}) where {T <: AbstractTree, B, N} =
    _matchbranchtype(T, B) && _matchnodetype(T, N)

"""
    _matchtreenametype(::Type{<:AbstractTree}, ::Type{X})

Does this tree type prefer the node or node label type provided?
"""
_matchtreenametype(::Type{T}, ::Type{TN}) where {T <: AbstractTree, TN} =
    (TN ≡ _treenametype(T))
                                                   
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
_gettree(treepair::Pair{Label, <: AbstractTree{OneTree}}) where Label =
    treepair[2]
@traitfn function _gettree(tree::T,
                           id::TN) where {T <: AbstractTree{OneTree}, TN;
                                          MatchTreeNameType{T, TN}}
    id == _gettreename(tree) || throw(BoundsError(tree, id))
    return tree
end

"""
    _getnodes(tree::AbstractTree{OneTree})

Returns a vector of nodes for a OneTree tree. Either _getnodes() must be
implemented for any OneTree tree type.
"""
function _getnodes end
_getnodes(::T) where T <: AbstractTree = error("No _getnodes() method for tree type $T")
_getnodes(tree::AbstractTree{OneTree}, order::TraversalOrder) =
    _traversal(tree, order)

"""
    _getnodenames(tree::AbstractTree{OneTree})

Returns a vector of node names for a OneTree tree. Can
be implemented for any OneTree tree type, especially PreferNodeObjects trees.
"""
function _getnodenames end
_getnodenames(::T) where T <: AbstractTree = error("No _getnodes() method for tree type $T")
@traitfn _getnodenames(tree::T, order::TraversalOrder) where
{T <: AbstractTree{OneTree}; !PreferNodeObjects{T}} = _getnodes(tree, order)
@traitfn _getnodenames(tree::T, order::TraversalOrder) where
{T <: AbstractTree{OneTree}; PreferNodeObjects{T}} =
    _getnodename.(tree, _getnodes(tree, order))

"""
    _nnodes(::AbstractTree)

Returns the number of nodes (internal nodes and leaves) in a single tree. May
be implemented for any OneTree tree type (otherwise infers from _getnodes()).
"""
_nnodes(tree::AbstractTree{OneTree}) = length(_getnodes(tree, anyorder))

"""
    _getleafnames(::AbstractTree)

Returns the leaf names of a tree. May be implemented for any
tree type (otherwise determined from _getnodenames() and _isleaf() functions).
"""
function _getleafnames end
_getleafnames(tree::AbstractTree{OneTree}, order::TraversalOrder) =
    [_getnodename(tree, node)
     for node in _traversal(tree, order) if _isleaf(tree, node)]
_getleafnames(tree::AbstractTree{ManyTrees}, order::TraversalOrder) =
    _getleafnames(_gettree(tree, first(_gettreenames(tree))), order)

"""
    _getleaves(::AbstractTree)

Returns the leaves (tips) of a single tree. May be implemented for any
OneTree type (otherwise determined from _getnodes() and _isleaf() functions).
"""
function _getleaves end
_getleaves(tree::AbstractTree{OneTree}, order::TraversalOrder) =
    [node for node in _traversal(tree, order) if _isleaf(tree, node)]

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
_getroots(tree::AbstractTree{OneTree, <: Rooted}) =
    [node for node in _getnodes(tree) if _isroot(tree, node)]

"""
    _getroot(::AbstractTree)

Returns the unique root of a rooted tree. May be implemented for any
OneTree type (otherwise determined from _getroots()).
"""
function _getroot end
function _getroot(tree::AbstractTree{OneTree, <: Rooted})
    @assert _nroots(tree) == 1 "More than one root for tree " *
        "(found $(_nroots(tree)))"
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

Returns the node or name associated with id (which could be a name,
a node or a pair) from a tree. Must be implemented for any PreferNodeObjects
tree and node label type.
"""
function _getnode end
@traitfn _getnode(tree::T,
                  node::N) where {RT, NL, N,
                                  T <: AbstractTree{OneTree, RT, NL, N};
                                  PreferNodeObjects{T}} = node
@traitfn _getnode(tree::T,
                  pair::Pair{NL, N}) where {RT, NL, N,
                                            T <: AbstractTree{OneTree,
                                                              RT, NL, N};
                                            PreferNodeObjects{T}} = pair[2]
@traitfn _getnode(tree::T,
                  nodename::NL) where {RT, NL,
                                       T <: AbstractTree{OneTree, RT, NL};
                                       !PreferNodeObjects{T}} = nodename
@traitfn _getnode(tree::T,
                  pair::Pair{NL, N}) where {RT, NL, N,
                                            T <: AbstractTree{OneTree,
                                                              RT, NL, N};
                                            !PreferNodeObjects{T}} = pair[1]

"""
    _getnodename(::AbstractTree, id)

Returns the name of a node associated with id (which could be a name, a node
or a pair) from a tree. Must be implemented for PreferNodeObjects tree types.
"""
function _getnodename end
_getnodename(::AbstractTree{OneTree, RT, NL, N},
             pair::Pair{NL, N}) where {RT, NL, N} = pair[1]
_getnodename(::AbstractTree{OneTree, RT, NL, N},
             nodename::NL) where {RT, NL, N} = nodename

"""
    _hasnode(tree::AbstractTree, node[name])

Does the tree contain this node? Must be implemented for any PreferNodeObjects
tree type with a node label.
"""
function _hasnode end
_hasnode(tree::AbstractTree{OneTree, RT, NL},
         node::NL) where {RT, NL} = node ∈ _getnodenames(tree)
_hasnode(tree::AbstractTree{OneTree, RT, NL, N},
         node::N) where {RT, NL, N} = node ∈ _getnodes(tree)

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

Returns the branch or name associated with id (which could be a name,
a branch or a pair) from a tree. Must be implemented for any PreferBranchObjects
tree and branch label type.
"""
function _getbranch end
@traitfn _getbranch(::T,
                    branch::B) where {RT, NL, N, B,
                                      T <: AbstractTree{OneTree, RT, NL, N, B};
                                      PreferBranchObjects{T}} = branch
@traitfn _getbranch(::T,
                    pair::Pair{Int, B}) where {RT, NL, N, B,
                                               T <: AbstractTree{OneTree, RT,
                                                                 NL, N, B};
                                               PreferBranchObjects{T}} = pair[2]
@traitfn _getbranch(::T,
                    branchname::Int) where {T <: AbstractTree{OneTree};
                                            !PreferBranchObjects{T}} =
                                                branchname
@traitfn _getbranch(::T,
                    pair::Pair{Int, B}) where {RT, NL, N, B,
                                               T <: AbstractTree{OneTree,
                                                                 RT, NL, N, B};
                                               !PreferBranchObjects{T}} =
                                                   pair[1]
@traitfn _getbranch(tree::T, src::N1, dst::N2) where
    {T <: AbstractTree{OneTree}, N1, N2; MatchNodeTypes{T, N1, N2}} =
    first(b for b in _getbranches(tree) if _src(tree, b) == src && _dst(tree, b) == dst)

@traitfn _getbranch(tree::T, src::N1, dst::N2) where
    {T <: AbstractTree{OneTree}, N1, N2; !MatchNodeTypes{T, N1, N2}} =
    _getbranch(tree, _getnode(tree, src), _getnode(tree, dst))

"""
    _getbranchname(::AbstractTree, id)

Returns the name of a branch associated with id (which could be a name, a branch
or a pair) from a tree. Must be implemented for PreferBranchObjects tree types.
"""
function _getbranchname end
_getbranchname(::AbstractTree{OneTree, RT, NL, N},
               pair::Pair{NL, N}) where {RT, NL, N} = pair[1]
_getbranchname(::AbstractTree{OneTree}, id::Int) = id

"""
    _getbranches(tree::AbstractTree)

Returns a vector of branches for a OneTree tree. Either _getbranches() or
_getbranchnames() must be implemented for any OneTree tree type.
"""
function _getbranches end
_getbranches(::T) where T = error("No _getbranches() method for tree type $T")

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
@traitfn _getbranchnames(tree::T) where {T <: AbstractTree{OneTree};
                                         !PreferBranchObjects{T}} =
                                             _getbranches(tree)
@traitfn _getbranchnames(tree::T) where {T <: AbstractTree{OneTree};
                                         PreferBranchObjects{T}} =
                                             [_getbranchname(tree, branch)
                                              for branch in _getbranches(tree)]

"""
    _hasbranch(tree::AbstractTree, node[name])

Does the tree contain this branch? Must be implemented for any
PreferBranchObjects tree type with a branch label.
"""
function _hasbranch end
_hasbranch(tree::AbstractTree{OneTree},
           branch::Int) = branch ∈ _getbranchnames(tree)
@traitfn _hasbranch(tree::T,
                    branch::B) where {RT, NL, N, B,
                                      T <: AbstractTree{OneTree, RT, NL, N, B};
                                      PreferBranchObjects{T}} =
                                          branch ∈ _getbranches(tree)

@traitfn function _hasbranch(tree::T, src::N1, dst::N2) where
    {T <: AbstractTree{OneTree}, N1, N2; MatchNodeTypes{T, N1, N2}}
    return any(_src(tree, b) == src && _dst(tree, b) == dst for b in _getbranches(tree))
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


"""
function _validate! end
_validate!(::AbstractTree) = true

"""
    _traversal(tree::AbstractTree, order::TraversalOrder, todo, sofar)

Return an iterable object containing nodes in given order - preorder, inorder,
postorder or breadthfirst
"""
function _traversal end
function _traversal(tree::T, order::TraversalOrder) where {T <: AbstractTree{OneTree, <: Rooted}}
    if order == anyorder
        return _getnodes(tree)
    else
        return _traversal(tree, order, collect(_getroots(tree)))
    end
end
function _traversal(tree::AbstractTree{OneTree, <: Rooted},
                    order::TraversalOrder, todo::Vector,
                    sofar::Vector = eltype(todo)[])
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
_isroot(tree::AbstractTree{OneTree, <: Rooted}, node) =
    !_hasinbound(tree, node)
_isroot(tree::AbstractTree{OneTree, Unrooted}, node) = false

"""
    _isinternal(tree::AbstractTree, node)

Is the node internal to the tree?
"""
function _isinternal end
_isinternal(tree::AbstractTree{OneTree, <: Rooted}, node) =
    _outdegree(tree, node) > 0 && _hasinbound(tree, node)
_isinternal(tree::AbstractTree{OneTree, Unrooted}, node) =
    _degree(tree, node) < 2 ? false : missing

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
_indegree(tree::AbstractTree{OneTree, <: Rooted}, node) =
    Int(_hasinbound(tree, node))
_indegree(tree::AbstractTree{OneTree, Unrooted}, node) =
    _degree(tree, node) == 0 ? 0 : missing

"""
    _hasinboundspace(tree::AbstractTree, node::AbstractNode)

Is there space for a new inbound connection on a node?
"""
function _hasinboundspace end
_hasinboundspace(tree::AbstractTree{OneTree}, node) =
    !_hasinbound(tree, node)

"""
    _outdegree(tree::AbstractTree, node::AbstractNode)

Out degree of node.
"""
function _outdegree end
_outdegree(tree::AbstractTree{OneTree, <: Rooted}, node) =
    length(_getoutbounds(tree, node))
_outdegree(tree::AbstractTree{OneTree, Unrooted},
           node::AbstractElt{Unrooted}) =
               _degree(tree, node) == 0 ? 0 : missing

"""
    _hasoutboundspace(tree::AbstractTree, node::AbstractNode)

Is there space for a new outbound connection on a node? Must be implemented if
a node has a limit on the number of outbound connections (eg for a binary tree)
"""
function _hasoutboundspace end
_hasoutboundspace(::AbstractTree{OneTree, <: Rooted}, _) = true

"""
    _hasspace(tree::AbstractTree, node::AbstractNode)

Is there space for a new connection on a node? Must be implemented if
a node has a limit on the number of connections (eg for a binary tree)
"""
function _hasspace end
_hasspace(tree::AbstractTree{OneTree, <: Rooted}, node) =
    _hasinboundspace(tree, node) | _hasoutboundspace(tree, node)
_hasspace(::AbstractTree{OneTree, Unrooted}, _) = true

"""
    _degree(tree::AbstractTree, node::AbstractNode)

Degree of node. Must be implemented for Unrooted nodes, otherwise can be
inferred from indegree and outdegree.
"""
function _degree end
_degree(tree::AbstractTree{OneTree, <: Rooted}, node) =
    _indegree(tree, node) + _outdegree(tree, node)

"""
    _hasinbound(tree::AbstractTree, node::AbstractNode)

Must be implemented for any AbstractNode subtype.
"""
function _hasinbound end
_hasinbound(::T, ::N) where {T, N} = error("No _hasinbound() function for types $T, $N")

"""
    _getinbound(tree::AbstractTree, node::AbstractNode)

Get the inbound connection. Must be implemented for any rooted AbstractNode subtype.
"""
function _getinbound end
_getinbound(::T, ::B) where {T, B} = error("No _getinbound() function for $T, $B")

"""
    _getparent(tree::AbstractTree, node)

Return the parent node for this node. Can be implemented for Rooted node types.
"""
function _getparent end
_getparent(tree::AbstractTree{OneTree}, node) =
    _getnode(tree, _src(tree, _getinbound(tree, node)))

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
_getoutbounds(::T, ::N) where {T, N} = error("No _getoutbounds() function for $T, $N")

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
_getchildren(tree::AbstractTree{OneTree, RT}, node::N) where {RT <: Rooted, N <: AbstractElt{RT}} =
    N[_dst(tree, branch) for branch in _getoutbounds(tree, node)]
_getchildren(tree::AbstractTree{OneTree, RT, NL}, node::NL) where {RT <: Rooted, NL} =
    NL[_dst(tree, branch) for branch in _getoutbounds(tree, node)]

"""
    _getconnections(tree::AbstractTree, node::AbstractNode)

Returns all of the connections of a node. Must be implemented for any unrooted
AbstractNode subtype, can be inferred from _getinbound and _getoutbounds for
a rooted node.
"""
function _getconnections end
_getconnections(tree::AbstractTree{OneTree, <: Rooted}, node) =
    _hasinbound(tree, node) ?
        append!([_getinbound(tree, node)], _getoutbounds(tree, node)) :
        _getoutbounds(tree, node)

"""
    _getsiblings(tree::AbstractTree, node::AbstractNode)

Returns all of the siblings of a node. May be implemented for any
AbstractNode subtype, can be inferred from _getparent and _getchildren for
a rooted node or _getconnections for an unrooted node.
"""
function _getsiblings end
_getsiblings(tree::AbstractTree{OneTree, <: Rooted}, node) =
    _hasinbound(tree, node) ?
    append!([_getparent(tree, node)], _getchildren(tree, node)) :
    _getchildren(tree, node)
_getsiblings(tree::AbstractTree{OneTree, Unrooted}, node) =
    [_conn(tree, branch, node) for branch in _getconnections(tree, node)]

"""
    _addconnection!(tree::AbstractTree, node::AbstractNode, branch)

Add a connection to an unrooted node. Must be implemented for any unrooted
AbstractNode subtype unless this happens when a branch is added.
"""
function _addconnection! end
_addconnection!(::AbstractTree{OneTree, <: Rooted}, _, _) =
    error("Direction of connection must be specified for a rooted node")

"""
    _removeconnection!(tree::AbstractTree, node::AbstractNode, branch)

Remove a connection from an unrooted node. Must be implemented for any Unrooted
AbstractNode subtype unless this happens when a branch is deleted.
"""
function _removeconnection! end
_removeconnection!(::AbstractTree{OneTree, <: Rooted}, _, _) =
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
_conns(tree::AbstractTree{OneTree, <: Rooted}, branch) =
    [_src(tree, branch), _dst(tree, branch)]

"""
    _conn(branch::AbstractBranch, exclude::AbstractNode)

Return the connection for a branch that isn't the `exclude` node. May be
implemented for any Unrooted AbstractBranch subtype, otherwise will use
_conns.
"""
function _conn end
function _conn(tree::AbstractTree{OneTree, <: Rooted}, branch, exclude)
    other = [node for node in _conns(tree, branch) if node !== exclude]
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
_getnodedata(::T, ::N) where {T, N} = error("No _getnodedata() function for $T, $N")
_getnodedata(tree::AbstractTree{OneTree}, node, label) =
    _getnodedata(tree, node)[label]
function _setnodedata! end
_setnodedata!(tree::AbstractTree{OneTree}, node, label, value) =
    (_getnodedata(tree, node)[label] = value)

function _getbranchdata end
_getbranchdata(::T, ::B) where {T, B} = error("No _getbranchdata() function for $T, $B")
_getbranchdata(tree::AbstractTree, branch, label) =
    _getbranchdata(tree, branch)[label]
function _setbranchdata! end
_setbranchdata!(tree::AbstractTree, branch, label, value) =
    (_getbranchdata(tree, branch)[label] = value)

function _getleafinfo end
function _setleafinfo! end

"""
    _gettreeinfo(tree::AbstractTree)
    _gettreeinfo(tree::AbstractTree, treename)

Returns the info data associated with the tree(s).
"""
function _gettreeinfo end
_gettreeinfo(tree::AbstractTree{OneTree}) =
    Dict(_gettreename(tree) => Dict{String, Any}())
function _gettreeinfo(tree::AbstractTree{OneTree}, treename)
    @assert _gettreename(tree) == treename "No tree called $treename"
    return Dict{String, Any}()
end

"""
    _resetleaves!(::AbstractTree)

Fixes leaf naming after creation or deletion of nodes or branches. Must be
implemented by tree types where this is handled separately.
"""
function _resetleaves! end
_resetleaves!(::AbstractTree) = true
