using Compat
using Phylo.API

abstract type AbstractTreeIterator{T <: AbstractTree} end

function iteratorsize(::Type{AbstractTreeIterator})
    return Base.HasLength()
end

function iteratoreltype(::Type{AbstractTreeIterator})
    return Base.HasEltype()
end

abstract type AbstractNodeIterator{T <: AbstractTree} <: AbstractTreeIterator{T} end

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

function length(ni::It) where It <: AbstractNodeIterator
    return ni.filterfn === nothing ? length(_getnodes(ni.tree)) :
        mapreduce(val -> ni.filterfn(_extractnode(ni.tree, val)) ? 1 : 0, +, 0, ni)
end


abstract type AbstractBranchIterator{T <: AbstractTree} <: AbstractTreeIterator{T} end

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

function length(bi::It) where It <: AbstractBranchIterator
    return bi.filterfn === nothing ? length(_getbranches(bi.tree)) :
        mapreduce(val -> bi.filterfn(_extractbranch(bi.tree, val)) ? 1 : 0, +, 0, bi)
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
    NodeIterator{T}(tree, Union{Function, Nothing}(filterfn))

eltype(ni::NodeIterator{T}) where T <: AbstractTree = nodetype(ni.tree)

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
    NodeNameIterator{T}(tree, Union{Function, Nothing}(filterfn))

eltype(ni::NodeNameIterator{T}) where T <: AbstractTree = nodenametype(ni.tree)

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
    BranchIterator{T}(tree, Union{Function, Nothing}(filterfn))

eltype(bi::BranchIterator{T}) where T <: AbstractTree = branchtype(bi.tree)

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
    BranchNameIterator{T}(tree, Union{Function, Nothing}(filterfn))

eltype(bi::BranchNameIterator{T}) where T <: AbstractTree = branchnametype(bi.tree)

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
