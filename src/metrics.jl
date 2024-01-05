"""
    mrca(tree::AbstractTree, target)

Returns the node within `tree` that is the Most Recent Common Ancestor of the 
leaves (or internal nodes) defined by `target`. `target` can be an iterator
over nodes or an AbstractArray of nodes. The return value has the same node
type as the elements of `target`.

### Examples
≡≡≡≡≡≡≡≡≡≡≡

julia> tree = open(parsenewick, Phylo.path("H1N1.newick"))
RootedTree with 507 tips, 1013 nodes and 1012 branches.
Leaf names are 227, 294, 295, 110, 390, ... [501 omitted] ... and 418


julia> tips = rand(collect(nodefilter(isleaf, tree)), 3)
3-element Vector{LinkNode{OneRoot, String, Dict{String, Any}, LinkBranch{OneRoot, String, Dict{String, Any}, Float64}}}:
 LinkNode 153, a tip of the tree with an incoming connection (branch 888).

 LinkNode 120, a tip of the tree with an incoming connection (branch 195).

 LinkNode 504, a tip of the tree with an incoming connection (branch 44).


julia> mrca(tree, tips)
LinkNode Node 1003, an internal node with 1 inbound and 2 outbound connections (branches 1001 and 999, 1000)
"""
function mrca(tree::AbstractTree, target)
    ancestors = getancestors(tree, first(target))
    ranks = Dict(j => i for (i, j) in enumerate(ancestors))
    checked = Set(ancestors)
    oldest = 1
    for species in target
        while !(species ∈ checked)
            push!(checked, species)
            species = getparent(tree, species)
        end
        oldest = max(oldest, get(ranks, species, 0))
    end
    return ancestors[oldest]
end

"""
    nodeheights(tree::Phylo.AbstractTree; onlyleaves = false, noleaves = false)

Returns an AxisArray of the height of all nodes in `tree` over the root node.
`onlyleaves` and `noleaves` filter the nodes for leaves and internal nodes,
respectively
"""
function nodeheights(tree::Phylo.AbstractTree; onlyleaves = false,
                     noleaves = false)
    function findheights!(clade::String, parentheight::Float64 = 0.0)
        myheight = parentheight
        if hasinbound(tree, clade)
            myheight += getlength(tree, getinbound(tree, clade))
        end
        height[clade] = myheight
        for ch in getchildren(tree, clade)
            findheights!(ch, myheight)
        end
    end
    names = getnodename.((tree,), traversal(tree, preorder))
    height = AxisArray(Vector{Float64}(undef, nnodes(tree)), x = names)
    root = getnodename(tree, getroot(tree))
    findheights!(root)
    onlyleaves && return height[filter(t -> isleaf(tree, t), names)]
    noleaves && return height[filter(t -> !isleaf(tree, t), names)]
    return height
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
    leaves = getleaves(tree)
    names = getnodename.(tree, leaves)
    return AxisArray([distance(tree, li, lj) for li in leaves, lj in leaves],
                     Axis{:x}(names), Axis{:y}(names))
end

"""
    heighttoroot(tree::AbstractTree, node)

Height of a node of the tree above the root
"""
function heighttoroot(tree::AbstractTree{OneTree, <:Rooted}, node)
    return mapreduce(branch -> getlength(tree, branch), +,
                     branchhistory(tree, node); init = 0.0)
end

"""
    heightstoroot(tree::AbstractTree)

Height of all of the leaves of the tree above the root
"""
function heightstoroot(tree::AbstractTree)
    return nodeheights(tree, onlyleaves = true)
end
