using Phylo.API

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

function eltype{T <: AbstractTree}(ni::NodeIterator{T})
    return _nodetype(T)
end

function length(ni::NodeIterator)
    if isnull(ni.filterfn)
        return length(_getnodes(ni.tree))
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

immutable NodeNameIterator{T <: AbstractTree}
    tree::T
    filterfn::Nullable{Function}
end

NodeNameIterator{T <: AbstractTree}(tree::T) =
    NodeNameIterator{T}(tree, Nullable{Function}())

NodeNameIterator{T <: AbstractTree}(tree::T, filterfn::Function) =
    NodeNameIterator{T}(tree, Nullable{Function}(filterfn))

function start(ni::NodeNameIterator)
    nodes = getnodes(ni.tree)
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

function next(ni::NodeNameIterator, state)
    nodes = getnodes(ni.tree)
    val, state = next(nodes, state)
    node = _extractnode(val)
    name = _extractnodename(val)
    
    if isnull(ni.filterfn) || done(ni, state)
        return name, state
    end

    fn = get(ni.filterfn)
    val, after = next(nodes, state)
    while !fn(_extractnode(val))
        state = after
        if done(nodes, state)
            return name, state
        end
        val, after = next(nodes, state)
    end
    
    return name, state
end

function done(ni::NodeNameIterator, state)
    return done(getnodes(ni.tree), state)
end

function iteratorsize(ni::Type{NodeNameIterator})
    return HasLength()
end

function iteratoreltype(ni::Type{NodeNameIterator})
    return HasEltype()
end

function eltype{T <: AbstractTree}(ni::NodeNameIterator{T})
    return keytype(_getnodes(ni.tree))
end

function length(ni::NodeNameIterator)
    nodes = _getnodes(ni.tree)
    if isnull(ni.filterfn)
        return length(nodes)
    end
    
    len = 0
    state = start(nodes)
    fn = get(ni.filterfn)
    while !done(nodes, state)
        val, state = next(nodes, state)
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

function eltype{T <: AbstractTree}(bi::BranchIterator{T})
    return _branchtype(T)
end

function length(bi::BranchIterator)
    return isnull(bi.filterfn) ? length(_getbranches(bi.tree)) :
        mapreduce(val -> fn(_extractbranch(val)) ? 1 : 0, +, 0, bi)
end

immutable BranchNameIterator{T <: AbstractTree}
    tree::T
    filterfn::Nullable{Function}
end

BranchNameIterator{T <: AbstractTree}(tree::T) =
    BranchNameIterator{T}(tree, Nullable{Function}())

BranchNameIterator{T <: AbstractTree}(tree::T, filterfn::Function) =
    BranchNameIterator{T}(tree, Nullable{Function}(filterfn))

function start(bi::BranchNameIterator)
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

function next(bi::BranchNameIterator, state)
    branches = _getbranches(bi.tree)
    val, state = next(branches, state)
    branch = _extractbranch(val)
    name = _extractbranchname(val)

    if isnull(bi.filterfn) || done(bi, state)
        return name, state
    end

    fn = get(bi.filterfn)
    val, after = next(branches, state)
    while !fn(_extractbranch(val))
        state = after
        if done(branches, state)
            return name, state
        end
        val, after = next(branches, state)
    end
    
    return name, state
end

function done(bi::BranchNameIterator, state)
    return done(_getbranches(bi.tree), state)
end

function iteratorsize(bi::Type{BranchNameIterator})
    return HasLength()
end

function iteratoreltype(bi::Type{BranchNameIterator})
    return HasEltype()
end

function eltype{T <: AbstractTree}(bi::BranchNameIterator{T})
    return keytype(_getbranches(bi.tree))
end

function length(bi::BranchNameIterator)
    return isnull(bi.filterfn) ? length(_getbranches(bi.tree)) :
        mapreduce(val -> fn(_extractbranch(val)) ? 1 : 0, +, 0, bi)
end
