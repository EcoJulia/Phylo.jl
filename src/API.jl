using Compat

@compat abstract type AbstractNode end
@compat abstract type AbstractBranch end
@compat abstract type AbstractTree{NodeLabel, BranchLabel} end

function _newlabel{Label <: Integer}(ids::Vector{Label}, _)
    return isempty(ids) ? 1 : maximum(ids) + 1
end

function _newlabel(names::Vector{String}, prefix)
    names = collect(Compat.Iterators.filter(n -> length(n) > length(prefix) &&
                                            n[1:length(prefix)]==prefix,
                                            names))
    start = length(names) + 1
    name = prefix * "$start"
    while (name ∈ names)
        start += 1
        name = prefix * "$start"
    end
    return name
end

function _extractnode{N <: AbstractNode}(node::N)
    return node
end

function _extractnode{L, N <: AbstractNode}(pair::Pair{L, N})
    return pair[2]
end

function _nodetype end

function _extractbranch{B <: AbstractBranch}(branch::B)
    return branch
end

function _extractbranch{L, B <: AbstractBranch}(pair::Pair{L, B})
    return pair[2]
end

function _branchtype end

# AbstractTree methods
"""
    _newbranchlabel(tree::AbstractTree)


"""
function _newbranchlabel end
function _newbranchlabel(tree::AbstractTree)
    return _newlabel(_getbranchnames(tree), "Branch ")
end

"""
    _addbranch!(tree::AbstractTree, source, target;
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
    _branch!(tree::AbstractTree, source, length::Float64;
             nodename = _newnodelabel(tree),
             branchname = _newbranchlabel(tree))


"""
function _branch!(tree::AbstractTree, source,
                  length::Float64 = NaN;
                  nodename = _newnodelabel(tree),
                  branchname = _newbranchlabel(tree))
    target = _addnode!(tree, nodename)
    _addbranch!(tree, source, target, length, branchname)
    return target
end

"""
    _newnodelabel(tree::AbstractTree)


"""
function _newnodelabel end
function _newnodelabel(tree::AbstractTree)
    return _newlabel(_getnodenames(tree), "Node ")
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
    _addnodes!(tree::AbstractTree, count::Integer)


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
_isleaf(node::AbstractNode) = _outdegree(node) == 0 && _hasinbound(node)

"""
    _isroot(node::AbstractNode)


"""
_isroot(node::AbstractNode) = !_hasinbound(node) && _outdegree(node) > 0

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
function _hasoutboundspace(node::AbstractNode)
    return _outdegree < 2
end

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
    _getoutbounds(node::AbstractNode)

Must be implemented for any AbstractNode subtype.
"""
function _getoutbounds end

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
    _getsource

Must be implemented for any AbstractBranch subtype.
"""
function _getsource end

"""
    _gettarget

Must be implemented for any AbstractBranch subtype.
"""
function _gettarget end

"""
    _getlength

Must be implemented for any AbstractBranch subtype.
"""
function _getlength end

"""
    _setsource!

Must be implemented for any AbstractBranch subtype.
"""
function _setsource! end

"""
    _settarget!

Must be implemented for any AbstractBranch subtype.
"""
function _settarget! end
