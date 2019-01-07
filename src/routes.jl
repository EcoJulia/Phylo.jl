using Compat: mapreduce

function _treehistory(tree::T, nodename::NL) where
    {RT <: Rooted, NL, N, B, T <: AbstractTree{OneTree, RT, NL, N, B}}
    branches, nodes = _treehistory(tree, getnode(tree, nodename))
    return [getbranchname(tree, b) for b in branches],
           [getnodename(tree, n) for n in nodes]
end

function _treehistory(tree::T, node::N) where
    {RT <: Rooted, NL, N, B, T <: AbstractTree{OneTree, RT, NL, N, B}}
    branches = B[]
    nodestoprocess = N[node]
    nodesprocessed = N[]
    while !isempty(nodestoprocess)
        nextnode = pop!(nodestoprocess)
        push!(nodesprocessed, nextnode)
        if hasinbound(tree, nextnode)
            push!(branches, getinbound(tree, nextnode))
            push!(nodestoprocess, getparent(tree, nextnode))
        end
    end
    return branches, nodesprocessed
end

function _treefuture(tree::T, nodename::NL) where
    {RT <: Rooted, NL, N, B, T <: AbstractTree{OneTree, RT, NL, N, B}}
    branches, nodes = _treefuture(tree, getnode(tree, nodename))
    return [getbranchname(tree, b) for b in branches],
           [getnodename(tree, n) for n in nodes]
end

function _treefuture(tree::T, node::N) where
    {RT <: Rooted, NL, N, B, T <: AbstractTree{OneTree, RT, NL, N, B}}
    branches = B[]
    nodestoprocess = N[node]
    nodesprocessed = N[]
    while !isempty(nodestoprocess)
        nextnode = pop!(nodestoprocess)
        push!(nodesprocessed, nextnode)
        append!(branches, getoutbounds(tree, nextnode))
        append!(nodestoprocess, getchildren(tree, nextnode))
    end
    return branches, nodesprocessed
end

"""
    branchhistory(tree::AbstractTree, node)

Find the branch route between a node on a tree and its root
"""
function branchhistory(tree::T, node::Union{NL, N}) where
    {RT <: Rooted, NL, N, B, T <: AbstractTree{OneTree, RT, NL, N, B}}
    return _treehistory(tree, node)[1]
end

"""
    nodehistory(tree::AbstractTree, node)

Find the node route between a node on a tree and its root
"""
function nodehistory(tree::T, node::Union{NL, N}) where
    {RT <: Rooted, NL, N, B, T <: AbstractTree{OneTree, RT, NL, N, B}}
    return _treehistory(tree, node)[2]
end

"""
    branchfuture(tree::AbstractTree, node)

Find the branches between a node on a tree and its leaves
"""
function branchfuture(tree::T, node::Union{NL, N}) where
    {RT <: Rooted, NL, N, B, T <: AbstractTree{OneTree, RT, NL, N, B}}
    return _treefuture(tree, node)[1]
end

"""
    nodefuture(tree::AbstractTree, node)

Find the nodes between a node on a tree and its leaves
"""
function nodefuture(tree::T, node::Union{NL, N}) where
    {RT <: Rooted, NL, N, B, T <: AbstractTree{OneTree, RT, NL, N, B}}
    return _treefuture(tree, node)[2]
end

"""
    branchfuture(tree::AbstractTree, node)

Find the branches between a node on a tree and its leaves
"""
function branchfuture(tree::T, node::NL) where {NL, BL, T <: AbstractTree{NL, BL}}
    return _treefuture(tree, node)[1]
end

"""
    nodefuture(tree::AbstractTree, node)

Find the nodes between a node on a tree and its leaves
"""
function nodefuture(tree::T, node::NL) where {NL, BL, T <: AbstractTree{NL, BL}}
    return _treefuture(tree, node)[2]
end

"""
    branchroute(tree::AbstractTree, node1, node2)

Find the branch route between two nodes on a tree
"""
function branchroute(tree::T, node1::NL, node2::NL) where
    {RT <: Rooted, NL, N, B, T <: AbstractTree{OneTree, RT, NL, N, B}}
    route = branchroute(tree, getnode(tree, node1), getnode(tree, node2))
    return [getbranchname(tree, b) for b in route]
end
function branchroute(tree::T, node1::N, node2::N) where
    {RT <: Rooted, NL, N, B, T <: AbstractTree{OneTree, RT, NL, N, B}}
    branches1, nodes1 = _treehistory(tree, node1)
    branches2, nodes2 = _treehistory(tree, node2)
    nodes1[end] == nodes2[end] ||
        return error("No route between nodes")
    common = branches1 ∩ branches2
    nodes1[end] == nodes2[end] || error("No route between nodes in tree")
    return append!([b for b in branches1 if b ∉ common],
                   [b for b in reverse(branches2) if b ∉ common])
end

"""
    noderoute(tree::AbstractTree, node1, node2)

Find the node route between two nodes on a tree
"""
function noderoute(tree::T, node1::NL, node2::NL) where
    {RT <: Rooted, NL, N, B, T <: AbstractTree{OneTree, RT, NL, N, B}}
    route = noderoute(tree, getnode(tree, node1), getnode(tree, node2))
    return [getnodename(tree, n) for n in route]
end
function noderoute(tree::T, node1::N, node2::N) where
    {RT <: Rooted, NL, N, B, T <: AbstractTree{OneTree, RT, NL, N, B}}
    branches1, nodes1 = _treehistory(tree, node1)
    branches2, nodes2 = _treehistory(tree, node2)
    nodes1[end] == nodes2[end] || error("No route between nodes in tree")
    common = nodes1[end]
    while min(length(nodes1), length(nodes2)) > 0 && nodes1[end] == nodes2[end]
        common = nodes1[end]
        pop!(nodes1)
        pop!(nodes2)
    end
    push!(nodes1, common)
    return append!(nodes1, reverse(nodes2))
end

"""
    distance(tree::AbstractTree, node1, node2)

Distance between two nodes on a tree
"""
function distance(tree::AbstractTree, node1, node2)
    branches = branchroute(tree, node1, node2)
    return mapreduce(branch -> getlength(tree, branch), +, branches;
    init = 0.0)
end

"""
    distances(tree::AbstractTree)

Pairwise distances between all leaf nodes on a tree
"""
function distances(tree::AbstractTree)
    leaves = getleaves(tree)
    return [distance(tree, li, lj) for li in leaves, lj in leaves]
end

"""
    height(tree::AbstractTree, node)

Height of a node of the tree above the root
"""
function heighttoroot(tree::AbstractTree, node)
    return mapreduce(branch -> getlength(tree, branch), +,
                     branchhistory(tree, node); init = 0.0)
end

"""
    heights(tree::AbstractTree)

Height of all of the leaves of the tree above the root
"""
function heightstoroot(tree::AbstractTree)
    return [heighttoroot(tree, leaf) for leaf in getleaves(tree)]
end
