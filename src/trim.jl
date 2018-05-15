"""
    getinternalnodes(t::AbstractTree)
Function to retrieve only the internal nodes from a tree, `t`, which does not
include tips or root.

"""
function getinternalnodes(t::AbstractTree)
    return collect(nodenamefilter(x->!isleaf(x)& !isroot(x), t))
end
"""
    drop_tip!(t::AbstractTree, tip::Vector{String})
Function to drop tips from a phylogenetic tree `t`, which are not found in
the vector of tip names, `tip`.

"""
function drop_tip!(t::AbstractTree, tip::Vector{String})
    tree_names = getleafnames(t)
    cut_names = setdiff(tree_names, tip)
    # Remove nodes that are not in tip names
    for i in cut_names[1:end]
        deletenode!(t, i)
    end
    # Remove nodes that are related to the cut tips
    while length(setdiff(collect(nodenamefilter(isleaf, t)), tip)) > 0
        nodes = setdiff(collect(nodenamefilter(isleaf, t)), tip)
        map(x -> deletenode!(t, x), nodes)
    end
    # Merge internal nodes that no longer have multiple children
    while sum(map(x-> length(getchildren(t, x)).< 2, getinternalnodes(t))) > 0
        inner_nodes = getinternalnodes(t)
        remove_nodes = find(map(x-> length(getchildren(t, x)).< 2, inner_nodes))
        for i in remove_nodes
            parent = getparent(t, inner_nodes[i])
            parentbranch = getinbound(getnode(t, inner_nodes[i]))

            child = getchildren(t, inner_nodes[i])[1]
            childbranch = getoutbounds(getnode(t, inner_nodes[i]))[1]

            len = distance(t, parent, child)

            deletebranch!(t, parentbranch)
            deletebranch!(t, childbranch)
            delete!(getnodes(t), inner_nodes[i])
            delete!(t.noderecords, inner_nodes[i])

            addbranch!(t, parent, child, len)
        end
    end
    # Remove root if it no longer has two branches
    root = collect(nodenamefilter(isroot, t))[1]
    if length(getchildren(t, root)) < 2
        deletenode!(t, root)
    end
    map(x -> delete!(t.leafinfos, x), cut_names)
end
