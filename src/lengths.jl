function _historyandroot{NL, BL}(tree::AbstractTree{NL, BL}, node::NL)
    branches = BL[]
    while hasinbound(tree, node)
        inbound = getinbound(tree, node)
        branches = push!(branches, inbound)
        node = getsource(tree, inbound)
    end
    return branches, node
end

"""
    treehistory(tree::AbstractTree, node)

Find the branch route between a node on a tree and its root
"""
function treehistory{NL, BL}(tree::AbstractTree{NL, BL}, node::NL)
    return _historyandroot(tree, node)[1]
end

"""
    treepath(tree::AbstractTree, node1, node2)

Find the branch route between two nodes on a tree
"""
function treepath{NL, BL}(tree::AbstractTree{NL, BL}, node1::NL, node2::NL)
    branches1, root1 = _historyandroot(tree, node1)
    branches2, root2 = _historyandroot(tree, node2)
    root1 == root2 ||
        return Nullable{Vector{BL}}()
    common = branches1 ∩ branches2
    return Nullable(append!(filter(b -> b ∉ common, branches1),
                            filter(b -> b ∉ common, reverse(branches2))))
end

"""
    distance(tree::AbstractTree, node1, node2)

Distance between two nodes on a tree
"""
function distance(tree::AbstractTree, node1, node2)
    branches = treepath(tree, node1, node2)
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
                     treehistory(tree, node))
end

"""
    heights(tree::AbstractTree)

Height of all of the leaves of the tree above the root 
"""
function heightstoroot(tree::AbstractTree)
    return [heighttoroot(tree, i) for i in NodeNameIterator(tree, isleaf)]
end
