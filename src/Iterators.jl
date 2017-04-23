using AbstractPhylo.API

import Base.start, Base.next, Base.done
import Base.iteratorsize, Base.iteratoreltype, Base.eltype, Base.length

immutable NodeIterator{T <: AbstractTree}
    tree::T
    filterfn::Nullable{Function}
end

NodeIterator{T <: AbstractTree}(tree::T) =
    NodeIterator{T}(tree, Nullable{Function}())

NodeIterator{T <: AbstractTree}(tree::T, filterfn::Function) =
    NodeIterator{T}(tree, Nullable{Function}(filterfn))

function start(ni::NodeIterator)
    nodes = _getnodes(ni.tree)
    state = start(nodes)
    
    if isnull(ni.filterfn) || done(nodes, state)
        return state
    end

    fn = get(ni.filterfn)
    val, after = next(nodes, state)
    while !fn(_extractnode(val))
        state = after
        if done(nodes, state)
            return state
        end
        val, after = next(nodes, state)
    end
    
    return state
end

function next(ni::NodeIterator, state)
    nodes = _getnodes(ni.tree)
    val, state = next(nodes, state)
    node = _extractnode(val)
    
    if isnull(ni.filterfn) || done(ni, state)
        return node, state
    end

    fn = get(ni.filterfn)
    val, after = next(nodes, state)
    while !fn(_extractnode(val))
        state = after
        if done(nodes, state)
            return node, state
        end
        val, after = next(nodes, state)
    end
    
    return node, state
end

function done(ni::NodeIterator, state)
    return done(_getnodes(ni.tree), state)
end

function iteratorsize(ni::Type{NodeIterator})
    return HasLength()
end

function iteratoreltype(ni::Type{NodeIterator})
    return HasEltype()
end

function eltype(ni::NodeIterator)
    return _nodetype(ni.tree)
end

function length(ni::NodeIterator)
    if isnull(ni.filterfn)
        return length(getnodes(ni.tree))
    end
    
    len = 0
    state = start(ni)
    fn = get(ni.filterfn)
    while !done(ni, state)
        val, state = next(ni, state)
        if fn(_extractnode(val))
            len += 1
        end
    end
    
    return len
end

immutable BranchIterator{T <: AbstractTree}
    tree::T
    filterfn::Nullable{Function}
end

BranchIterator{T <: AbstractTree}(tree::T) =
    BranchIterator{T}(tree, Nullable{Function}())

BranchIterator{T <: AbstractTree}(tree::T, filterfn::Function) =
    BranchIterator{T}(tree, Nullable{Function}(filterfn))

function start(bi::BranchIterator)
    branches = _getbranches(bi.tree)
    state = start(branches)
    
    if isnull(bi.filterfn) || done(branches, state)
        return state
    end

    fn = get(bi.filterfn)
    val, after = next(branches, state)
    while !fn(_extractbranch(val))
        state = after
        if done(branches, state)
            return state
        end
        val, after = next(branches, state)
    end
    
    return state
end

function next(bi::BranchIterator, state)
    branches = _getbranches(bi.tree)
    val, state = next(branches, state)
    branch = _extractbranch(val)
    
    if isnull(bi.filterfn) || done(bi, state)
        return branch, state
    end

    fn = get(bi.filterfn)
    val, after = next(branches, state)
    while !fn(_extractbranch(val))
        state = after
        if done(branches, state)
            return branch, state
        end
        val, after = next(branches, state)
    end
    
    return branch, state
end

function done(bi::BranchIterator, state)
    return done(_getbranches(bi.tree), state)
end

function iteratorsize(bi::Type{BranchIterator})
    return HasLength()
end

function iteratoreltype(bi::Type{BranchIterator})
    return HasEltype()
end

function eltype(bi::BranchIterator)
    return _branchtype(bi.tree)
end

function length(bi::BranchIterator)
    if isnull(bi.filterfn)
        return length(getbranches(bi.tree))
    end
    
    len = 0
    state = start(bi)
    fn = get(bi.filterfn)
    while !done(bi, state)
        val, state = next(bi, state)
        if fn(_extractbranch(val))
            len += 1
        end
    end
    
    return len
end
