using IterableTables
using IterableTables: getiterator
using Compat: findall

"""
    getinternalnodes(t::AbstractTree)
Function to retrieve only the internal nodes from a tree, `t`, which does not
include tips or root.

"""
function getinternalnodes(t::AbstractTree)
    return collect(nodenamefilter(x->!isleaf(x) & !isroot(x), t))
end
"""
    droptips!(t::T, tips::Vector{NL}) where {NL, BL, T <: AbstractTree{NL, BL}}
Function to drop tips from a phylogenetic tree `t`, which are found in
the vector of tip names, `tips`.

"""
function droptips!(t::T, tips::Vector{NL}) where {NL, BL, T <: AbstractTree{NL, BL}}
    tree_names = getleafnames(t)
    keep_tips = setdiff(tree_names, tips)
    # Remove nodes that are not in tip names
    for i in tips
        deletenode!(t, i)
    end
    # Remove nodes that are related to the cut tips
    while length(setdiff(collect(nodenamefilter(isleaf, t)), keep_tips)) > 0
        nodes = setdiff(collect(nodenamefilter(isleaf, t)), keep_tips)
        map(x -> deletenode!(t, x), nodes)
    end
    # Merge internal nodes that no longer have multiple children
    while any(map(x -> length(getchildren(t, x)) .< 2,
                  getinternalnodes(t)))
        inner_nodes = getinternalnodes(t)
        remove_nodes = findall(map(x->length(getchildren(t, x)) .< 2,
                                   inner_nodes))
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

    if !isempty(getleafinfo(t))
        li = leafinfotype(t)(Iterators.filter(line -> line[1] ∉ tips,
                                              getiterator(getleafinfo(t))))
        setleafinfo!(t, li)
    end
    return tips
end

"""
    keeptips!(t::T, tips::Vector{NL}) where {NL, BL, T <: AbstractTree{NL, BL}}
Function to keep only the tips in a phylogenetic tree, `t`, that are found in
the vector of tip names, `tip`.

"""
function keeptips!(t::T, tips::Vector{NL}) where {NL, BL, T <: AbstractTree{NL, BL}}
    tree_names = getleafnames(t)
    cut_names = setdiff(tree_names, tips)
    droptips!(t, cut_names)
end
