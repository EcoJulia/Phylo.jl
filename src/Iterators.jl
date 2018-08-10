using Compat
using Compat: Random, mapreduce
import Compat.IteratorSize, Base.length, Compat.IteratorEltype, Base.eltype

using Phylo.API

abstract type AbstractTreeIterator{T <: AbstractTree} end

function IteratorSize(::Type{AbstractTreeIterator})
    return HasLength()
end

function IteratorEltype(::Type{AbstractTreeIterator})
    return HasEltype()
end

abstract type AbstractNodeIterator{T <: AbstractTree} <: AbstractTreeIterator{T} end

function length(ni::It) where It <: AbstractNodeIterator
    return ni.filterfn === nothing ? length(_getnodes(ni.tree)) :
        mapreduce(val -> ni.filterfn(_extractnode(ni.tree, val)) ? 1 : 0,
                  +, ni; init = 0)
end

abstract type AbstractBranchIterator{T <: AbstractTree} <: AbstractTreeIterator{T} end

function length(bi::It) where It <: AbstractBranchIterator
    return bi.filterfn === nothing ? length(_getbranches(bi.tree)) :
        mapreduce(val -> bi.filterfn(_extractbranch(bi.tree, val)) ? 1 : 0,
                  +, bi; init = 0)
end

struct NodeIterator{T <: AbstractTree} <: AbstractNodeIterator{T}
    tree::T
    filterfn::Union{Function, Nothing}
end
"""
    nodeiter(tree::AbstractTree)

Returns an iterator over the nodes of any tree.
"""
nodeiter(tree::T) where T <: AbstractTree =
    NodeIterator{T}(tree, nothing)

"""
    nodefilter(filterfn::Function, tree::AbstractTree)

Returns an iterator over the nodes of any tree, where the
`AbstractNode` is filtered by the function `filterfn`.
"""
nodefilter(filterfn::Function, tree::T) where T <: AbstractTree =
    NodeIterator{T}(tree, filterfn)

eltype(ni::NodeIterator{T}) where T <: AbstractTree = nodetype(ni.tree)

struct NodeNameIterator{T <: AbstractTree} <: AbstractNodeIterator{T}
    tree::T
    filterfn::Union{Function, Nothing}
end

"""
    nodenameiter(tree::AbstractTree)

Returns an iterator over the names of the nodes of any tree.
"""
nodenameiter(tree::T) where T <: AbstractTree =
    NodeNameIterator{T}(tree, nothing)

"""
    nodenamefilter(filterfn::Function, tree::AbstractTree)

Returns an iterator over the nodenames of any tree, where the
`AbstractNode` itself is filtered by the function `filterfn`.
"""
nodenamefilter(filterfn::Function, tree::T) where T <: AbstractTree =
    NodeNameIterator{T}(tree, filterfn)

eltype(ni::NodeNameIterator{T}) where T <: AbstractTree = nodenametype(ni.tree)

struct BranchIterator{T <: AbstractTree} <: AbstractBranchIterator{T}
    tree::T
    filterfn::Union{Function, Nothing}
end

"""
    branchiter(tree::AbstractTree)

Returns an iterator over the branches of any tree.
"""
branchiter(tree::T) where T <: AbstractTree =
    BranchIterator{T}(tree, nothing)

"""
    branchfilter(filterfn::Function, tree::AbstractTree)

Returns an iterator over the branches of any tree, where the
`AbstractBranch` is filtered by the function `filterfn`.
"""
branchfilter(filterfn::Function, tree::T) where T <: AbstractTree =
    BranchIterator{T}(tree, filterfn)

eltype(bi::BranchIterator{T}) where T <: AbstractTree = branchtype(bi.tree)

struct BranchNameIterator{T <: AbstractTree} <: AbstractBranchIterator{T}
    tree::T
    filterfn::Union{Function, Nothing}
end

"""
    branchnameiter(tree::AbstractTree)

Returns an iterator over the names of branches of any tree.
"""
branchnameiter(tree::T) where T <: AbstractTree =
    BranchNameIterator{T}(tree, nothing)

"""
    branchnamefilter(filterfn::Function, tree::AbstractTree)

Returns an iterator over the names of the branches of any tree, where
the `AbstractBranch` is filtered by the function `filterfn`.
"""
branchnamefilter(filterfn::Function, tree::T) where T <: AbstractTree =
    BranchNameIterator{T}(tree, filterfn)

eltype(bi::BranchNameIterator{T}) where T <: AbstractTree = branchnametype(bi.tree)

@static if VERSION >= v"0.7.0"
import Base.iterate
function iterate(ni::NodeIterator, state = nothing)
    nodes = _getnodes(ni.tree)
    if state === nothing
        result = iterate(nodes)
    else
        result = iterate(nodes, state)
    end

    result === nothing && return nothing

    if ni.filterfn === nothing
        return _extractnode(ni.tree, result[1]), result[2]
    end

    val, state = result
    node = _extractnode(ni.tree, val)
    while !ni.filterfn(node)
        result = iterate(nodes, state)
        result === nothing && return nothing
        val, state = result
        node = _extractnode(ni.tree, val)
    end

    return node, state
end

function iterate(ni::NodeNameIterator, state = nothing)
    nodes = _getnodes(ni.tree)
    if state === nothing
        result = iterate(nodes)
    else
        result = iterate(nodes, state)
    end

    result === nothing && return nothing

    if ni.filterfn === nothing
        return _extractnodename(ni.tree, result[1]), result[2]
    end

    val, state = result
    node = _extractnode(ni.tree, val)
    while !ni.filterfn(node)
        result = iterate(nodes, state)
        result === nothing && return nothing
        val, state = result
        node = _extractnode(ni.tree, val)
    end

    name = _extractnodename(ni.tree, val)
    return name, state
end

function iterate(bi::BranchIterator, state = nothing)
    branches = _getbranches(bi.tree)
    if state === nothing
        result = iterate(branches)
    else
        result = iterate(branches, state)
    end

    result === nothing && return nothing

    if bi.filterfn === nothing
        return _extractbranch(bi.tree, result[1]), result[2]
    end

    val, state = result
    branch = _extractbranch(bi.tree, val)
    while !bi.filterfn(branch)
        result = iterate(branches, state)
        result === nothing && return nothing
        val, state = result
        branch = _extractbranch(bi.tree, val)
    end

    return branch, state
end

function iterate(bi::BranchNameIterator, state = nothing)
    branches = _getbranches(bi.tree)
    if state === nothing
        result = iterate(branches)
    else
        result = iterate(branches, state)
    end

    result === nothing && return nothing

    if bi.filterfn === nothing
        return _extractbranchname(bi.tree, result[1]), result[2]
    end

    val, state = result
    branch = _extractbranch(bi.tree, val)
    while !bi.filterfn(branch)
        result = iterate(branches, state)
        result === nothing && return nothing
        val, state = result
        branch = _extractbranch(bi.tree, val)
    end
    name = _extractbranchname(bi.tree, val)

    return name, state
end

else
import Base.start, Base.next, Base.done
function start(ni::It) where It <: AbstractNodeIterator
    nodes = _getnodes(ni.tree)
    state = start(nodes)

    if ni.filterfn === nothing || done(nodes, state)
        return state
    end

    val, after = next(nodes, state)
    while !ni.filterfn(_extractnode(ni.tree, val))
        state = after
        if done(nodes, state)
            return state
        end
        val, after = next(nodes, state)
    end

    return state
end

function done(ni::It, state) where It <: AbstractNodeIterator
    return done(_getnodes(ni.tree), state)
end

function next(ni::NodeIterator, state)
    nodes = _getnodes(ni.tree)
    val, state = next(nodes, state)
    node = _extractnode(ni.tree, val)

    if ni.filterfn === nothing || done(ni, state)
        return node, state
    end

    val, after = next(nodes, state)
    while !ni.filterfn(_extractnode(ni.tree, val))
        state = after
        if done(nodes, state)
            return node, state
        end
        val, after = next(nodes, state)
    end

    return node, state
end

function next(ni::NodeNameIterator, state)
    nodes = getnodes(ni.tree)
    val, state = next(nodes, state)
    node = _extractnode(ni.tree, val)
    name = _extractnodename(ni.tree, val)

    if ni.filterfn === nothing || done(ni, state)
        return name, state
    end

    val, after = next(nodes, state)
    while !ni.filterfn(_extractnode(ni.tree, val))
        state = after
        if done(nodes, state)
            return name, state
        end
        val, after = next(nodes, state)
    end

    return name, state
end

function start(bi::It) where It <: AbstractBranchIterator
    branches = _getbranches(bi.tree)
    state = start(branches)

    if bi.filterfn === nothing || done(branches, state)
        return state
    end

    val, after = next(branches, state)
    while !bi.filterfn(_extractbranch(bi.tree, val))
        state = after
        if done(branches, state)
            return state
        end
        val, after = next(branches, state)
    end

    return state
end

function done(bi::It, state) where It <: AbstractBranchIterator
    return done(_getbranches(bi.tree), state)
end

function next(bi::BranchIterator, state)
    branches = _getbranches(bi.tree)
    val, state = next(branches, state)
    branch = _extractbranch(bi.tree, val)

    if bi.filterfn === nothing || done(bi, state)
        return branch, state
    end

    val, after = next(branches, state)
    while !bi.filterfn(_extractbranch(bi.tree, val))
        state = after
        if done(branches, state)
            return branch, state
        end
        val, after = next(branches, state)
    end

    return branch, state
end

function next(bi::BranchNameIterator, state)
    branches = _getbranches(bi.tree)
    val, state = next(branches, state)
    branch = _extractbranch(bi.tree, val)
    name = _extractbranchname(bi.tree, val)

    if bi.filterfn === nothing || done(bi, state)
        return name, state
    end

    val, after = next(branches, state)
    while !bi.filterfn(_extractbranch(bi.tree, val))
        state = after
        if done(branches, state)
            return name, state
        end
        val, after = next(branches, state)
    end

    return name, state
end

end
