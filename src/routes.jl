using SimpleTraits
using Phylo.API
using AxisArrays

@traitfn function _treehistory(tree::T,
                               node::N1) where
                  {N1, NL, N, B,
                   T <: AbstractTree{OneTree, <:Rooted, NL, N, B};
                   PreferBranchObjects{T}}
    branches = B[]
    nodestoprocess = N1[node]
    nodesprocessed = N1[]
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

@traitfn function _treehistory(tree::T,
                               node::N1) where
                  {N1, NL, N, B,
                   T <: AbstractTree{OneTree, <:Rooted, NL, N, B};
                   !PreferBranchObjects{T}}
    branches = Int[]
    nodestoprocess = N1[node]
    nodesprocessed = N1[]
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

@traitfn function _treefuture(tree::T,
                              node::N1) where
                  {N1, NL, N, B,
                   T <: AbstractTree{OneTree, <:Rooted, NL, N, B};
                   PreferBranchObjects{T}}
    branches = B[]
    nodestoprocess = N1[node]
    nodesprocessed = N1[]
    while !isempty(nodestoprocess)
        nextnode = pop!(nodestoprocess)
        push!(nodesprocessed, nextnode)
        append!(branches, getoutbounds(tree, nextnode))
        append!(nodestoprocess, getchildren(tree, nextnode))
    end
    return branches, nodesprocessed
end

@traitfn function _treefuture(tree::T,
                              node::N1) where
                  {N1, NL, N, B,
                   T <: AbstractTree{OneTree, <:Rooted, NL, N, B};
                   !PreferBranchObjects{T}}
    branches = Int[]
    nodestoprocess = N1[node]
    nodesprocessed = N1[]
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
function branchhistory end
@traitfn function branchhistory(tree::T,
                                node::N) where
                  {N, T <: AbstractTree{OneTree, <:Rooted}; MatchNodeType{T, N}}
    return _treehistory(tree, node)[1]
end
@traitfn function branchhistory(tree::T,
                                node::N) where
                  {N,
                   T <: AbstractTree{OneTree, <:Rooted, N}; !MatchNodeType{T,
                                                                           N}}
    return _treehistory(tree, getnode(tree, node))[1]
end

"""
    nodehistory(tree::AbstractTree, node)

Find the node route between a node on a tree and its root
"""
function nodehistory end
@traitfn function nodehistory(tree::T,
                              node::N) where
                  {N, T <: AbstractTree{OneTree, <:Rooted}; MatchNodeType{T, N}}
    return _treehistory(tree, node)[2]
end
@traitfn function nodehistory(tree::T,
                              node::N) where
                  {N,
                   T <: AbstractTree{OneTree, <:Rooted, N}; !MatchNodeType{T,
                                                                           N}}
    return getnodename.(tree, _treehistory(tree, getnode(tree, node))[2])
end

"""
    branchfuture(tree::AbstractTree, node)

Find the branches between a node on a tree and its leaves
"""
function branchfuture end
@traitfn function branchfuture(tree::T,
                               node::N) where
                  {N, T <: AbstractTree{OneTree, <:Rooted}; MatchNodeType{T, N}}
    return _treefuture(tree, node)[1]
end
@traitfn function branchfuture(tree::T,
                               node::N) where
                  {N,
                   T <: AbstractTree{OneTree, <:Rooted, N}; !MatchNodeType{T,
                                                                           N}}
    return _treefuture(tree, getnode(tree, node))[1]
end

"""
    nodefuture(tree::AbstractTree, node)

Find the nodes between a node on a tree and its leaves
"""
function nodefuture end
@traitfn function nodefuture(tree::T,
                             node::N) where
                  {N, T <: AbstractTree{OneTree, <:Rooted}; MatchNodeType{T, N}}
    return _treefuture(tree, node)[2]
end
@traitfn function nodefuture(tree::T,
                             node::N) where
                  {N,
                   T <: AbstractTree{OneTree, <:Rooted, N}; !MatchNodeType{T,
                                                                           N}}
    return getnodename.(tree, _treefuture(tree, getnode(tree, node))[2])
end

"""
    branchroute(tree::AbstractTree, node1, node2)

Find the branch route between two nodes on a tree
"""
function branchroute end
@traitfn function branchroute(tree::T, node1::N,
                              node2::N) where
                  {RT, N, T <: AbstractTree{OneTree, RT}; !MatchNodeType{T, N}}
    return branchroute(tree, getnode(tree, node1), getnode(tree, node2))
end
@traitfn function branchroute(tree::T, node1::N,
                              node2::N) where
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
@traitfn function noderoute(tree::T, node1::N,
                            node2::N) where
                  {RT, N,
                   T <: AbstractTree{OneTree, RT, N}; !MatchNodeType{T, N}}
    route = noderoute(tree, getnode(tree, node1), getnode(tree, node2))
    return [getnodename(tree, n) for n in route]
end
@traitfn function noderoute(tree::T, node1::N,
                            node2::N) where
                  {N, T <: AbstractTree{OneTree}; MatchNodeType{T, N}}
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
