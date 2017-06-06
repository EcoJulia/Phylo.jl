function _treehistory{NL, BL}(tree::AbstractTree{NL, BL}, node::NL)
    branches = BL[]
    nodes = NL[node]
    while hasinbound(tree, node)
        inbound = getinbound(tree, node)
        branches = push!(branches, inbound)
        node = getsource(tree, inbound)
        nodes = push!(nodes, node)
    end
    return branches, nodes
end

"""
    branchhistory(tree::AbstractTree, node)

Find the branch route between a node on a tree and its root
"""
function branchhistory{NL, BL}(tree::AbstractTree{NL, BL}, node::NL)
    return _treehistory(tree, node)[1]
end

"""
    nodehistory(tree::AbstractTree, node)

Find the node route between a node on a tree and its root
"""
function nodehistory{NL, BL}(tree::AbstractTree{NL, BL}, node::NL)
    return _treehistory(tree, node)[2]
end

"""
    branchroute(tree::AbstractTree, node1, node2)

Find the branch route between two nodes on a tree
"""
function branchroute{NL, BL}(tree::AbstractTree{NL, BL}, node1::NL, node2::NL)
    branches1, nodes1 = _treehistory(tree, node1)
    branches2, nodes2 = _treehistory(tree, node2)
    nodes1[end] == nodes2[end] ||
        return Nullable{Vector{BL}}()
    common = branches1 ∩ branches2
    return Nullable(append!(filter(b -> b ∉ common, branches1),
                            filter(b -> b ∉ common, reverse(branches2))))
end

"""
    noderoute(tree::AbstractTree, node1, node2)

Find the node route between two nodes on a tree
"""
function noderoute{NL, BL}(tree::AbstractTree{NL, BL}, node1::NL, node2::NL)
    branches1, nodes1 = _treehistory(tree, node1)
    branches2, nodes2 = _treehistory(tree, node2)
    nodes1[end] == nodes2[end] ||
        return Nullable{Vector{NL}}()
    common = nodes1[end]
    while min(length(nodes1), length(nodes2)) > 0 && nodes1[end] == nodes2[end]
        common = nodes1[end]
        pop!(nodes1)
        pop!(nodes2)
    end
    push!(nodes1, common)
    return Nullable(append!(nodes1, reverse(nodes2)))
end

"""
    distance(tree::AbstractTree, node1, node2)

Distance between two nodes on a tree
"""
function distance(tree::AbstractTree, node1, node2)
    branches = branchroute(tree, node1, node2)
    return isnull(branches) ? Inf :
        mapreduce(branch -> getlength(tree, branch), +, 0.0, get(branches))
end

"""
    distances(tree::AbstractTree)

Pairwise distances between all leaf nodes on a tree
"""
function distances(tree::AbstractTree)
    leaves = NodeNameIterator(tree, isleaf)
    return [distance(tree, i, j) for i in leaves, j in leaves]
end

"""
    height(tree::AbstractTree, node)

Height of a node of the tree above the root 
"""
function heighttoroot(tree::AbstractTree, node)
    return mapreduce(branch -> getlength(tree, branch), +, 0.0,
                     branchhistory(tree, node))
end

"""
    heights(tree::AbstractTree)

Height of all of the leaves of the tree above the root 
"""
function heightstoroot(tree::AbstractTree)
    return [heighttoroot(tree, i) for i in NodeNameIterator(tree, isleaf)]
end
