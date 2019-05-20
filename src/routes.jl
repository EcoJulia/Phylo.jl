using SimpleTraits
using Phylo.API
using Compat: mapreduce
using AxisArrays

@traitfn function _treehistory(tree::T, node::N1) where
    {N1, NL, N, B, T <: AbstractTree{OneTree, <:Rooted, NL, N, B};
     PreferBranchObjects{T}}
    branches = B[]
    nodestoprocess = N1[node]
    nodesprocessed = N1[]
    while !isempty(nodestoprocess)
        nextnode = pop!(nodestoprocess)
        push!(nodesprocessed, nextnode)
        if hasinbound(tree, nextnode)
            push!(branches, _getinbound(tree, nextnode))
            push!(nodestoprocess, _getparent(tree, nextnode))
        end
    end
    return branches, nodesprocessed
end

@traitfn function _treehistory(tree::T, node::N1) where
    {N1, NL, N, B, T <: AbstractTree{OneTree, <:Rooted, NL, N, B};
     !PreferBranchObjects{T}}
    branches = Int[]
    nodestoprocess = N1[node]
    nodesprocessed = N1[]
    while !isempty(nodestoprocess)
        nextnode = pop!(nodestoprocess)
        push!(nodesprocessed, nextnode)
        if hasinbound(tree, nextnode)
            push!(branches, _getinbound(tree, nextnode))
            push!(nodestoprocess, _getparent(tree, nextnode))
        end
    end
    return branches, nodesprocessed
end

@traitfn function _treefuture(tree::T, node::N1) where
    {N1, NL, N, B, T <: AbstractTree{OneTree, <:Rooted, NL, N, B};
     PreferBranchObjects{T}}
    branches = B[]
    nodestoprocess = N1[node]
    nodesprocessed = N1[]
    while !isempty(nodestoprocess)
        nextnode = pop!(nodestoprocess)
        push!(nodesprocessed, nextnode)
        append!(branches, _getoutbounds(tree, nextnode))
        append!(nodestoprocess, _getchildren(tree, nextnode))
    end
    return branches, nodesprocessed
end

@traitfn function _treefuture(tree::T, node::N1) where
    {N1, NL, N, B, T <: AbstractTree{OneTree, <:Rooted, NL, N, B};
     !PreferBranchObjects{T}}
    branches = Int[]
    nodestoprocess = N1[node]
    nodesprocessed = N1[]
    while !isempty(nodestoprocess)
        nextnode = pop!(nodestoprocess)
        push!(nodesprocessed, nextnode)
        append!(branches, _getoutbounds(tree, nextnode))
        append!(nodestoprocess, _getchildren(tree, nextnode))
    end
    return branches, nodesprocessed
end

"""
    branchhistory(tree::AbstractTree, node)

Find the branch route between a node on a tree and its root
"""
function branchhistory end
@traitfn function branchhistory(tree::T, node::N) where
    {N, T <: AbstractTree{OneTree, <:Rooted}; MatchNodeType{T, N}}
    return _treehistory(tree, node)[1]
end
@traitfn function branchhistory(tree::T, node::N) where
    {N, T <: AbstractTree{OneTree, <:Rooted, N}; !MatchNodeType{T, N}}
    return _treehistory(tree, getnode(tree, node))[1]
end

"""
    nodehistory(tree::AbstractTree, node)

Find the node route between a node on a tree and its root
"""
function nodehistory end
@traitfn function nodehistory(tree::T, node::N) where
    {N, T <: AbstractTree{OneTree, <:Rooted}; MatchNodeType{T, N}}
    return _treehistory(tree, node)[2]
end
@traitfn function nodehistory(tree::T, node::N) where
    {N, T <: AbstractTree{OneTree, <:Rooted, N}; !MatchNodeType{T, N}}
    return getnodename.(tree, _treehistory(tree, getnode(tree, node))[2])
end

"""
    branchfuture(tree::AbstractTree, node)

Find the branches between a node on a tree and its leaves
"""
function branchfuture end
@traitfn function branchfuture(tree::T, node::N) where
    {N, T <: AbstractTree{OneTree, <:Rooted}; MatchNodeType{T, N}}
    return _treefuture(tree, node)[1]
end
@traitfn function branchfuture(tree::T, node::N) where
    {N, T <: AbstractTree{OneTree, <:Rooted, N}; !MatchNodeType{T, N}}
    return _treefuture(tree, getnode(tree, node))[1]
end

"""
    nodefuture(tree::AbstractTree, node)

Find the nodes between a node on a tree and its leaves
"""
function nodefuture end
@traitfn function nodefuture(tree::T, node::N) where
    {N, T <: AbstractTree{OneTree, <:Rooted}; MatchNodeType{T, N}}
    return _treefuture(tree, node)[2]
end
@traitfn function nodefuture(tree::T, node::N) where
    {N, T <: AbstractTree{OneTree, <:Rooted, N}; !MatchNodeType{T, N}}
    return getnodename.(tree, _treefuture(tree, getnode(tree, node))[2])
end

"""
    branchroute(tree::AbstractTree, node1, node2)

Find the branch route between two nodes on a tree
"""
function branchroute end
@traitfn function branchroute(tree::T, node1::N, node2::N) where
    {RT, N, T <: AbstractTree{OneTree, RT, N}; !MatchNodeType{T, N}}
    return branchroute(tree, getnode(tree, node1), getnode(tree, node2))
end
@traitfn function branchroute(tree::T, node1::N, node2::N) where
    {RT, N, T <: AbstractTree{OneTree, RT}; MatchNodeType{T, N}}
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
function noderoute end
@traitfn function noderoute(tree::T, node1::N, node2::N) where
    {RT, N, T <: AbstractTree{OneTree, RT, N}; !MatchNodeType{T, N}}
    route = noderoute(tree, getnode(tree, node1), getnode(tree, node2))
    return [getnodename(tree, n) for n in route]
end
@traitfn function noderoute(tree::T, node1::N, node2::N) where
    {RT, N, T <: AbstractTree{OneTree}; MatchNodeType{T, N}}
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
    branches = branchroute(tree, getnode(tree, node1), getnode(tree, node2))
    return mapreduce(branch -> getlength(tree, branch), +, branches; init = 0.0)
end

"""
    distances(tree::AbstractTree)

Pairwise distances between all leaf nodes on a tree
"""
function distances(tree::AbstractTree)
    leaves = [node for node in traversal(tree, preorder) if isleaf(tree, node)]
    names = getnodename.(tree, leaves)
    return AxisArray([distance(tree, li, lj) for li in leaves, lj in leaves],
                     Axis{:x}(names), Axis{:y}(names))
end

"""
    height(tree::AbstractTree, node)

Height of a node of the tree above the root
"""
function heighttoroot(tree::AbstractTree{OneTree, <:Rooted}, node)
    return mapreduce(branch -> getlength(tree, branch), +,
                     branchhistory(tree, node); init = 0.0)
end

"""
    heights(tree::AbstractTree)

Height of all of the leaves of the tree above the root
"""
function heightstoroot(tree::AbstractTree)
    leaves = [node for node in traversal(tree, preorder) if isleaf(tree, node)]
    return AxisArray([heighttoroot(tree, leaf) for leaf in leaves],
                     Axis{:x}(getnodename.(tree, leaves)))
end
