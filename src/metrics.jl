# SPDX-License-Identifier: BSD-2-Clause

using LinearAlgebra: diag

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
3-element Vector{RecursiveNode{OneRoot, String, Dict{String, Any}, Dict{String, Any}, PolytomousBranching, Float64}}:
 leaf node '414'
 leaf node '232'
 leaf node '216'

julia> mrca(tree, tips)
internal node 'Node 980'
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

Pairwise distances between all leaf nodes on a tree, as an `AxisArray` over the
leaf names in the order of `getleaves(tree)`.

For a rooted tree whose branches all have lengths, each distance is computed as
``d_{ij} = h_i + h_j - 2 h_{\\mathrm{mrca}(i, j)}`` from the leaf heights and the
shared paths of [`vcv`](@ref), in a single pass over the tree. Other trees fall
back to summing the branches on the route between each pair of leaves.
"""
function distances(tree::AbstractTree)
    leaves = getleaves(tree)
    names = getnodename.(tree, leaves)
    if _hasalllengths(tree)
        shared = _sharedpaths(tree, names)
        h = diag(shared)
        return AxisArray(h .+ h' .- 2 .* shared, Axis{:x}(names),
                         Axis{:y}(names))
    end
    return AxisArray([distance(tree, li, lj) for li in leaves, lj in leaves],
                     Axis{:x}(names), Axis{:y}(names))
end

"""
    vcv(tree::AbstractTree{OneTree, OneRoot})

The phylogenetic variance-covariance matrix of the leaves of a rooted tree: for
each pair of leaves, the length of the path from the root that they share (the
height of their most recent common ancestor), with each leaf's height on the
diagonal. This is the covariance of trait values at the leaves under Brownian
motion with unit rate. Returns an `AxisArray` over the leaf names in the order
of `getleaves(tree)`. Every branch must have a length.
"""
function vcv(tree::AbstractTree{OneTree, OneRoot})
    _hasalllengths(tree) ||
        throw(ArgumentError("vcv needs every branch of the tree to have a length"))
    names = getnodename.(tree, getleaves(tree))
    return AxisArray(_sharedpaths(tree, names), Axis{:x}(names),
                     Axis{:y}(names))
end

# The fast path needs a single rooted tree with a length on every branch
_hasalllengths(::AbstractTree) = false
function _hasalllengths(tree::AbstractTree{OneTree, OneRoot})
    return all(branch -> haslength(tree, branch), getbranches(tree))
end

# Shared root-path lengths between the leaves `names` (the height of each pair's
# most recent common ancestor, and each leaf's own height on the diagonal), filled
# in one postorder pass: at every internal node, the leaves below different
# children meet at that node. Each pair is written once, so the cost is
# proportional to the number of pairs.
function _sharedpaths(tree::AbstractTree{OneTree, OneRoot}, names)
    index = Dict(name => i for (i, name) in enumerate(names))
    height = nodeheights(tree)
    shared = zeros(Float64, length(names), length(names))
    below = Dict{eltype(names), Vector{Int}}()
    for node in traversal(tree, postorder)
        name = getnodename(tree, node)
        if isleaf(tree, node)
            i = index[name]
            shared[i, i] = height[name]
            below[name] = [i]
        else
            h = height[name]
            leaves = Int[]
            for child in getchildren(tree, name)
                clade = pop!(below, child)
                for i in leaves, j in clade
                    shared[i, j] = shared[j, i] = h
                end
                append!(leaves, clade)
            end
            below[name] = leaves
        end
    end
    return shared
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
