function mrca(tree, target)
    ancestors = getancestors(tree, first(target))
    ranks = Dict(j => i for (i,j) in enumerate(ancestors))
    checked = Set(ancestors)
    oldest = 1
    for species in target
        while !(species ∈ checked)
            push!(checked, species)
            species = getparent(tree, species)
        end
        oldest = max(oldest, get(ranks, species, 0))
    end
    ancestors[oldest]
end

function nodeheights(tree::Phylo.AbstractTree; onlyleaves = false, noleaves = false)
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
    names = getnodename.((tree, ), traversal(tree, preorder))
    height = AxisArray(Vector{Float64}(undef, nnodes(tree)), x = names)
    root = getnodename(tree, getroot(tree))
    findheights!(root)
    onlyleaves && return height[filter(t->isleaf(tree, t), names)]
    noleaves && return height[filter(t->!isleaf(tree, t), names)]
    height
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
    return nodeheights(tree, onlyleaves = true)
end
