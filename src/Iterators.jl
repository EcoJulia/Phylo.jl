# SPDX-License-Identifier: BSD-2-Clause

import Base: length, eltype
using Base: HasLength, HasEltype
using Phylo.API

# For a single tree
import Base.IteratorSize
function IteratorSize(::Type{<:AbstractTree{OneTree}})
    return HasLength()
end

import Base.IteratorEltype
function IteratorEltype(::Type{<:AbstractTree{OneTree}})
    return HasEltype()
end
length(::AbstractTree{OneTree}) = 1
eltype(ni::T) where {T <: AbstractTree} = T

abstract type AbstractTreeIterator{T <: AbstractTree} end

function IteratorSize(::Type{AbstractTreeIterator})
    return HasLength()
end

function IteratorEltype(::Type{AbstractTreeIterator})
    return HasEltype()
end

abstract type AbstractNodeIterator{T <: AbstractTree} <: AbstractTreeIterator{T} end

function length(ni::It) where {It <: AbstractNodeIterator}
    return isnothing(ni.filterfn) ? _nnodes(ni.tree) :
           count(val -> ni.filterfn(ni.tree, _getnode(ni.tree, val)), ni)
end

abstract type AbstractBranchIterator{T <: AbstractTree} <:
              AbstractTreeIterator{T} end

function length(bi::It) where {It <: AbstractBranchIterator}
    return isnothing(bi.filterfn) ? _nbranches(bi.tree) :
           count(val -> bi.filterfn(bi.tree, _getbranch(bi.tree, val)), bi)
end

"""
    NodeIterator

The struct representing an iterator for nodes of a phylogenetic tree
"""
struct NodeIterator{T <: AbstractTree} <: AbstractNodeIterator{T}
    tree::T
    filterfn::Union{Function, Nothing}
end
"""
    nodeiter(tree::AbstractTree)

Returns an iterator over the nodes of any tree.
"""
nodeiter(tree::T) where {T <: AbstractTree} = NodeIterator{T}(tree, nothing)

"""
    nodefilter(filterfn::Function, tree::AbstractTree)

Returns an iterator over the nodes of any tree, where the
`AbstractNode` is filtered by the function `filterfn`.
"""
nodefilter(filterfn::Function, tree::T) where {T <: AbstractTree} = NodeIterator{T}(tree,
                                                                                    filterfn)

eltype(ni::NodeIterator{T}) where {T <: AbstractTree} = nodetype(T)

"""
    NodeNameIterator

The struct representing an iterator for nodenames of a phylogenetic tree
"""
struct NodeNameIterator{T <: AbstractTree} <: AbstractNodeIterator{T}
    tree::T
    filterfn::Union{Function, Nothing}
end

"""
    nodenameiter(tree::AbstractTree)

Returns an iterator over the names of the nodes of any tree.
"""
nodenameiter(tree::T) where {T <: AbstractTree} = NodeNameIterator{T}(tree,
                                                                      nothing)

"""
    nodenamefilter(filterfn::Function, tree::AbstractTree)

Returns an iterator over the nodenames of any tree, where the
`AbstractNode` itself is filtered by the function `filterfn`.
"""
nodenamefilter(filterfn::Function, tree::T) where {T <: AbstractTree} = NodeNameIterator{T}(tree,
                                                                                            filterfn)

eltype(ni::NodeNameIterator{T}) where {T <: AbstractTree} = nodenametype(T)

"""
    BranchIterator

The struct representing an iterator for branches of a phylogenetic tree
"""
struct BranchIterator{T <: AbstractTree} <: AbstractBranchIterator{T}
    tree::T
    filterfn::Union{Function, Nothing}
end

"""
    branchiter(tree::AbstractTree)

Returns an iterator over the branches of any tree.
"""
branchiter(tree::T) where {T <: AbstractTree} = BranchIterator{T}(tree, nothing)

"""
    branchfilter(filterfn::Function, tree::AbstractTree)

Returns an iterator over the branches of any tree, where the
`AbstractBranch` is filtered by the function `filterfn`.
"""
branchfilter(filterfn::Function, tree::T) where {T <: AbstractTree} = BranchIterator{T}(tree,
                                                                                        filterfn)

eltype(bi::BranchIterator{T}) where {T <: AbstractTree} = branchtype(T)

"""
    BranchNameIterator

The struct representing an iterator for branchnames of a phylogenetic tree
"""
struct BranchNameIterator{T <: AbstractTree} <: AbstractBranchIterator{T}
    tree::T
    filterfn::Union{Function, Nothing}
end

"""
    branchnameiter(tree::AbstractTree)

Returns an iterator over the names of branches of any tree.
"""
branchnameiter(tree::T) where {T <: AbstractTree} = BranchNameIterator{T}(tree,
                                                                          nothing)

"""
    branchnamefilter(filterfn::Function, tree::AbstractTree)

Returns an iterator over the names of the branches of any tree, where
the `AbstractBranch` is filtered by the function `filterfn`.
"""
branchnamefilter(filterfn::Function, tree::T) where {T <: AbstractTree} = BranchNameIterator{T}(tree,
                                                                                                filterfn)

eltype(bi::BranchNameIterator{T}) where {T <: AbstractTree} = branchnametype(T)

import Base: iterate
function iterate(tree::AbstractTree, state = nothing)
    if isnothing(state)
        return first(gettrees(tree)), 1
    elseif ntrees(tree) > state
        return collect(gettrees(tree))[state], state + 1
    else
        return nothing
    end
end

function iterate(ni::NodeIterator, state = nothing)
    nodes = _getnodes(ni.tree)
    if isnothing(state)
        result = iterate(nodes)
    else
        result = iterate(nodes, state)
    end

    isnothing(result) && return nothing

    if isnothing(ni.filterfn)
        return _getnode(ni.tree, result[1]), result[2]
    end

    val, state = result
    node = _getnode(ni.tree, val)
    while !ni.filterfn(ni.tree, node)
        result = iterate(nodes, state)
        isnothing(result) && return nothing
        val, state = result
        node = _getnode(ni.tree, val)
    end

    return node, state
end

function iterate(ni::NodeNameIterator, state = nothing)
    nodes = _getnodes(ni.tree)
    if isnothing(state)
        result = iterate(nodes)
    else
        result = iterate(nodes, state)
    end

    isnothing(result) && return nothing

    if isnothing(ni.filterfn)
        return _getnodename(ni.tree, result[1]), result[2]
    end

    val, state = result
    node = _getnode(ni.tree, val)
    while !ni.filterfn(ni.tree, node)
        result = iterate(nodes, state)
        isnothing(result) && return nothing
        val, state = result
        node = _getnode(ni.tree, val)
    end

    name = _getnodename(ni.tree, val)
    return name, state
end

function iterate(bi::BranchIterator, state = nothing)
    branches = _getbranches(bi.tree)
    if isnothing(state)
        result = iterate(branches)
    else
        result = iterate(branches, state)
    end

    isnothing(result) && return nothing

    if isnothing(bi.filterfn)
        return _getbranch(bi.tree, result[1]), result[2]
    end

    val, state = result
    branch = _getbranch(bi.tree, val)
    while !bi.filterfn(bi.tree, branch)
        result = iterate(branches, state)
        isnothing(result) && return nothing
        val, state = result
        branch = _getbranch(bi.tree, val)
    end

    return branch, state
end

function iterate(bi::BranchNameIterator, state = nothing)
    branches = _getbranches(bi.tree)
    if isnothing(state)
        result = iterate(branches)
    else
        result = iterate(branches, state)
    end

    isnothing(result) && return nothing

    if isnothing(bi.filterfn)
        return _getbranchname(bi.tree, result[1]), result[2]
    end

    val, state = result
    branch = _getbranch(bi.tree, val)
    while !bi.filterfn(bi.tree, branch)
        result = iterate(branches, state)
        isnothing(result) && return nothing
        val, state = result
        branch = _getbranch(bi.tree, val)
    end
    name = _getbranchname(bi.tree, val)

    return name, state
end
