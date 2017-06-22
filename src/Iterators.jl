using Compat
using Phylo.API

@compat abstract type AbstractTreeIterator{T <: AbstractTree} end

function iteratorsize(::Type{AbstractTreeIterator})
    return Base.HasLength()
end

function iteratoreltype(::Type{AbstractTreeIterator})
    return Base.HasEltype()
end

@compat abstract type AbstractNodeIterator{T <: AbstractTree} <: AbstractTreeIterator{T} end

function start{It <: AbstractNodeIterator}(ni::It)
    nodes = _getnodes(ni.tree)
    state = start(nodes)
    
    if isnull(ni.filterfn) || done(nodes, state)
        return state
    end

    fn = get(ni.filterfn)
    val, after = next(nodes, state)
    while !fn(_extractnode(ni.tree, val))
        state = after
        if done(nodes, state)
            return state
        end
        val, after = next(nodes, state)
    end
    
    return state
end

function done{It <: AbstractNodeIterator}(ni::It, state)
    return done(_getnodes(ni.tree), state)
end

function length{It <: AbstractNodeIterator}(ni::It)
    return isnull(ni.filterfn) ? length(_getnodes(ni.tree)) :
        mapreduce(val -> get(ni.filterfn)(_extractnode(ni.tree, val)) ? 1 : 0, +, 0, ni)
end


@compat abstract type AbstractBranchIterator{T <: AbstractTree} <: AbstractTreeIterator{T} end

function start{It <: AbstractBranchIterator}(bi::It)
    branches = _getbranches(bi.tree)
    state = start(branches)
    
    if isnull(bi.filterfn) || done(branches, state)
        return state
    end

    fn = get(bi.filterfn)
    val, after = next(branches, state)
    while !fn(_extractbranch(bi.tree, val))
        state = after
        if done(branches, state)
            return state
        end
        val, after = next(branches, state)
    end
    
    return state
end

function done{It <: AbstractBranchIterator}(bi::It, state)
    return done(_getbranches(bi.tree), state)
end

function length{It <: AbstractBranchIterator}(bi::It)
    return isnull(bi.filterfn) ? length(_getbranches(bi.tree)) :
        mapreduce(val -> get(bi.filterfn)(_extractbranch(bi.tree, val)) ? 1 : 0, +, 0, bi)
end

immutable NodeIterator{T <: AbstractTree} <: AbstractNodeIterator{T}
    tree::T
    filterfn::Nullable{Function}
end

NodeIterator{T <: AbstractTree}(tree::T) =
    NodeIterator{T}(tree, Nullable{Function}())

NodeIterator{T <: AbstractTree}(tree::T, filterfn::Function) =
    NodeIterator{T}(tree, Nullable{Function}(filterfn))

eltype{T <: AbstractTree}(ni::NodeIterator{T}) = nodetype(ni.tree)

function next(ni::NodeIterator, state)
    nodes = _getnodes(ni.tree)
    val, state = next(nodes, state)
    node = _extractnode(ni.tree, val)
    
    if isnull(ni.filterfn) || done(ni, state)
        return node, state
    end

    fn = get(ni.filterfn)
    val, after = next(nodes, state)
    while !fn(_extractnode(ni.tree, val))
        state = after
        if done(nodes, state)
            return node, state
        end
        val, after = next(nodes, state)
    end
    
    return node, state
end

immutable NodeNameIterator{T <: AbstractTree} <: AbstractNodeIterator{T}
    tree::T
    filterfn::Nullable{Function}
end

NodeNameIterator{T <: AbstractTree}(tree::T) =
    NodeNameIterator{T}(tree, Nullable{Function}())

NodeNameIterator{T <: AbstractTree}(tree::T, filterfn::Function) =
    NodeNameIterator{T}(tree, Nullable{Function}(filterfn))

eltype{T <: AbstractTree}(ni::NodeNameIterator{T}) = nodenametype(ni.tree)

function next(ni::NodeNameIterator, state)
    nodes = getnodes(ni.tree)
    val, state = next(nodes, state)
    node = _extractnode(ni.tree, val)
    name = _extractnodename(ni.tree, val)
    
    if isnull(ni.filterfn) || done(ni, state)
        return name, state
    end

    fn = get(ni.filterfn)
    val, after = next(nodes, state)
    while !fn(_extractnode(ni.tree, val))
        state = after
        if done(nodes, state)
            return name, state
        end
        val, after = next(nodes, state)
    end
    
    return name, state
end

immutable BranchIterator{T <: AbstractTree} <: AbstractBranchIterator{T}
    tree::T
    filterfn::Nullable{Function}
end

BranchIterator{T <: AbstractTree}(tree::T) =
    BranchIterator{T}(tree, Nullable{Function}())

BranchIterator{T <: AbstractTree}(tree::T, filterfn::Function) =
    BranchIterator{T}(tree, Nullable{Function}(filterfn))

eltype{T <: AbstractTree}(bi::BranchIterator{T}) = branchtype(bi.tree)

function next(bi::BranchIterator, state)
    branches = _getbranches(bi.tree)
    val, state = next(branches, state)
    branch = _extractbranch(bi.tree, val)
    
    if isnull(bi.filterfn) || done(bi, state)
        return branch, state
    end

    fn = get(bi.filterfn)
    val, after = next(branches, state)
    while !fn(_extractbranch(bi.tree, val))
        state = after
        if done(branches, state)
            return branch, state
        end
        val, after = next(branches, state)
    end
    
    return branch, state
end

immutable BranchNameIterator{T <: AbstractTree} <: AbstractBranchIterator{T}
    tree::T
    filterfn::Nullable{Function}
end

BranchNameIterator{T <: AbstractTree}(tree::T) =
    BranchNameIterator{T}(tree, Nullable{Function}())

BranchNameIterator{T <: AbstractTree}(tree::T, filterfn::Function) =
    BranchNameIterator{T}(tree, Nullable{Function}(filterfn))

eltype{T <: AbstractTree}(bi::BranchNameIterator{T}) = branchnametype(bi.tree)

function next(bi::BranchNameIterator, state)
    branches = _getbranches(bi.tree)
    val, state = next(branches, state)
    branch = _extractbranch(bi.tree, val)
    name = _extractbranchname(bi.tree, val)

    if isnull(bi.filterfn) || done(bi, state)
        return name, state
    end

    fn = get(bi.filterfn)
    val, after = next(branches, state)
    while !fn(_extractbranch(bi.tree, val))
        state = after
        if done(branches, state)
            return name, state
        end
        val, after = next(branches, state)
    end
    
    return name, state
end
