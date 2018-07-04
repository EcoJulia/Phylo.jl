using Phylo

function _newlabel(ids::Vector{Label}) where Label <: Integer
    return isempty(ids) ? 1 : maximum(ids) + 1
end

function _newlabel(names::Vector{String}, prefix)
    names = collect(Iterators.filter(n -> length(n) > length(prefix) &&
                                     n[1:length(prefix)]==prefix, names))
    start = length(names) + 1
    name = prefix * "$start"
    while (name ∈ names)
        start += 1
        name = prefix * "$start"
    end
    return name
end

function _ntrees(::AbstractTree)
    return 1
end

function _extractnode(::T, node::N) where {T <: AbstractTree, N <: AbstractNode}
    return node
end

function _extractnode(tree::T, nodename::NL) where {NL, BL, T <: AbstractTree{NL, BL}}
    return getnode(tree, nodename)
end

function _extractnode(::T, pair::Pair{NL, N}) where {NL, BL, T <: AbstractTree{NL, BL},
                                                     N <: AbstractNode}
    return pair[2]
end

function _extractnodename(::T, pair::Pair{NL, N}) where {NL, BL, T <: AbstractTree{NL, BL},
                                                         N <: AbstractNode}
    return pair[1]
end

function _extractbranch(::T, branch::B) where {T <: AbstractTree, B <: AbstractBranch}
    return branch
end

function _extractbranch(tree::T, branchname::BL) where {NL, BL, T <: AbstractTree{NL, BL}}
    return getbranch(tree, branchname)
end

function _extractbranch(::T, pair::Pair{BL, B}) where {NL, BL, T <: AbstractTree{NL, BL},
                                                       B <: AbstractBranch}
    return pair[2]
end

function _extractbranchname(::T, pair::Pair{BL, B}) where {NL, BL, T <: AbstractTree{NL, BL},
                                                           B <: AbstractBranch}
    return pair[1]
end

# AbstractTree methods
"""
    _nodetype(::AbstractTree)

Returns type of nodes in a tree.
"""
function _nodetype end

"""
    _branchtype(::AbstractTree)

Returns type of branches in a tree.
"""
function _branchtype end


"""
    _newbranchlabel(tree::AbstractTree)


"""
function _newbranchlabel end

function _newbranchlabel(tree::AbstractTree{NL, String}) where NL
    return _newlabel(_getbranchnames(tree), "Branch ")
end

function _newbranchlabel(tree::AbstractTree{NL, I}) where {NL, I <: Integer}
    return _newlabel(_getbranchnames(tree))
end

"""
    _addbranch!(tree::AbstractTree, source, destination;
                length::Float64 = NaN,
                branchname = _newbranchlabel(tree))

Must be implemented for any AbstractTree subtype.
"""
function _addbranch! end

"""
    _deletebranch!(tree::AbstractTree, branchname)

Must be implemented for any AbstractTree subtype.
"""
function _deletebranch! end

"""
    _branch!(tree::AbstractTree, source, length::Float64, destination, branchname)


"""
function _branch!(tree::AbstractTree, source, length::Float64, destination, branchname)
    _addbranch!(tree, source, _addnode!(tree, destination), length, branchname)
    return destination
end

"""
    _newnodelabel(tree::AbstractTree)


"""
function _newnodelabel end

function _newnodelabel(tree::T) where {BL, T <: AbstractTree{String, BL}}
    return _newlabel(_getnodenames(tree), "Node ")
end

function _newnodelabel(tree::T) where {I <: Integer, BL, T <: AbstractTree{I, BL}}
    return _newlabel(_getnodenames(tree))
end

"""
    _getnodenames(tree::AbstractTree)

Must be implemented for any AbstractTree subtype.
"""
function _getnodenames end

"""
    _getnodes(tree::AbstractTree)

Must be implemented for any AbstractTree subtype.
"""
function _getnodes end

"""
    _addnode!(tree::AbstractTree, nodename)

Must be implemented for any AbstractTree subtype.
"""
function _addnode! end

"""
    _addnodes!(tree::AbstractTree, nodenames::AbstractVector)
    _addnodes!(tree:AbstractTree, count::Integer)


"""
function _addnodes! end

function _addnodes!(tree::AbstractTree, nodenames::AbstractVector)
    return map(name -> _addnode!(tree, name), nodenames)
end

function _addnodes!(tree::AbstractTree, count::Integer)
    return map(name -> addnode!(tree), 1:count)
end

"""
    _deletenode!(tree::AbstractTree, nodename)

Must be implemented for any AbstractTree subtype.
"""
function _deletenode! end

"""
    _hasnode(tree::AbstractTree, nodename)

Must be implemented for any AbstractTree subtype.
"""
function _hasnode end

"""
    _getnode(tree::AbstractTree, nodename)

Must be implemented for any AbstractTree subtype.
"""
function _getnode end

"""
    _getbranchnames(tree::AbstractTree)

Must be implemented for any AbstractTree subtype.
"""
function _getbranchnames end

"""
    _getbranches(tree::AbstractTree)

Must be implemented for any AbstractTree subtype.
"""
function _getbranches end

"""
    _hasbranch(tree::AbstractTree, branchname)

Must be implemented for any AbstractTree subtype.
"""
function _hasbranch end

"""
    _getbranch(tree::AbstractTree, branchname)

Must be implemented for any AbstractTree subtype.
"""
function _getbranch end

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
    _isleaf(node::AbstractNode)


"""
_isleaf(node::AbstractNode) = _outdegree(node) == 0

"""
    _isroot(node::AbstractNode)


"""
_isroot(node::AbstractNode) = !_hasinbound(node)

"""
    _isinternal(node::AbstractNode)


"""
_isinternal(node::AbstractNode) = _outdegree(node) > 0 && _hasinbound(node)

"""
    _isunattached(node::AbstractNode)


"""
_isunattached(node::AbstractNode) = _outdegree(node) == 0 && !_hasinbound(node)

"""
    _indegree(node::AbstractNode)


"""
_indegree(node::AbstractNode) = _hasinbound(node) ? 1 : 0

"""
    _hasinboundspace(node::AbstractNode)


"""
_hasinboundspace(node::AbstractNode) = !_hasinbound(node)

"""
    _outdegree(node::AbstractNode)

Must be implemented for any AbstractNode subtype.
"""
function _outdegree end

"""
    _hasoutboundspace(node::AbstractNode)


"""
function _hasoutboundspace end

"""
    _hasinbound(node::AbstractNode)

Must be implemented for any AbstractNode subtype.
"""
function _hasinbound end

"""
    _getinbound(node::AbstractNode)

Must be implemented for any AbstractNode subtype.
"""
function _getinbound end

"""
    _setinbound!(node::AbstractNode, inbound)

Must be implemented for any AbstractNode subtype.
"""
function _setinbound! end

"""
    _deleteinbound!(node::AbstractNode, inbound)

Must be implemented for any AbstractNode subtype.
"""
function _deleteinbound! end

"""
    _getoutbounds(node::AbstractNode)

Must be implemented for any AbstractNode subtype.
"""
function _getoutbounds end

"""
    _addinbound!(node::AbstractNode, outbound)

Must be implemented for any AbstractNode subtype.
"""
function _addoutbound! end

"""
    _deleteoutbound!(node::AbstractNode, outbound)

Must be implemented for any AbstractNode subtype.
"""
function _deleteoutbound! end

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
    _src

Return source node for a branch. Must be implemented for any
AbstractBranch subtype.
"""
function _src end

"""
    _dst

Return destination node for a branch. Must be implemented for any
AbstractBranch subtype.
"""
function _dst end

"""
    _getlength

Return length of a branch. May be implemented for any AbstractBranch
subtype.
"""
function _getlength end
function _getlength(::AbstractBranch)
    return NaN
end

"""
    _setsrc!

Set source node for a branch. Must be implemented for any
AbstractBranch subtype.
"""
function _setsrc! end

"""
    _setdst!

Set destination node for a graph. Must be implemented for any
AbstractBranch subtype.
"""
function _setdst! end

#  - _getleafnames()
function _getleafnames(tree::AbstractTree)
    return collect(nodenamefilter(_isleaf, tree))
end

function _getleafinfo end
function _setleafinfo! end
function _getnoderecord end
function _setnoderecord! end
function _resetleaves! end
function _clearrootheight! end
function _setnode! end
function _setbranch! end
function _leafinfotype end
function _nleaves end
