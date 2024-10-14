# SPDX-License-Identifier: BSD-2-Clause

using IterableTables
using IterableTables: getiterator
using Phylo.API

"""
    getinternalnodes(t::AbstractTree)

Function to retrieve only the internal nodes from a tree, `t`, which does not
include tips or root.
"""
getinternalnodes(t::AbstractTree) = collect(nodenamefilter(isinternal, t))

"""
    droptips!(tree::AbstractTree{OneTree}, tips; keep = false)

Function to drop tips from a phylogenetic tree `tree`, which are found in
the vector of tips or tip names, `tips`. `keep` determines whether to keep
internal and root nodes that now only have one child (default is `false`).
Internal nodes with no children will always be removed.
"""
function droptips! end
@traitfn function droptips!(tree::T,
                            tips::AbstractVector{N};
                            keep = false) where
                  {N, RT, NL,
                   T <: AbstractTree{OneTree, RT, NL}; !MatchNodeType{T, N}}
    return isempty(tips) ? NL[] :
           droptips!(tree, [getnode(tree, tip) for tip in tips], keep = keep)
end

@traitfn function droptips!(tree::T,
                            tips::AbstractVector{N};
                            keep = false) where
                  {N, RT, NL,
                   T <: AbstractTree{OneTree, RT, NL}; MatchNodeType{T, N}}
    keep_tips = N[tip for tip in getleaves(tree) if tip ∉ tips]
    tipnames = NL[getnodename(tree, tip) for tip in tips]

    # Remove nodes that are not in tip names
    for tip in tips
        _deletenode!(tree, tip)
    end

    # Remove nodes that are related to the cut tips
    while length(setdiff(collect(nodefilter(isleaf, tree)), keep_tips)) > 0
        nodes = setdiff(collect(nodefilter(isleaf, tree)), keep_tips)
        map(x -> deletenode!(tree, x), nodes)
    end

    if !keep
        # Merge internal nodes that no longer have multiple children
        while any(map(x -> length(getchildren(tree, x)) .< 2,
                      getinternalnodes(tree)))
            inner_nodes = getinternalnodes(tree)
            remove_nodes = findall(map(x -> length(getchildren(tree, x)) .< 2,
                                       inner_nodes))
            for i in remove_nodes
                parent = getparent(tree, inner_nodes[i])
                parentbranch = getinbound(tree, inner_nodes[i])

                child = getchildren(tree, inner_nodes[i])[1]
                childbranch = getoutbounds(tree, getnode(tree, inner_nodes[i]))[1]

                len = distance(tree, parent, child)

                deletebranch!(tree, parentbranch)
                deletebranch!(tree, childbranch)
                deletenode!(tree, inner_nodes[i])

                createbranch!(tree, parent, child, len)
            end
        end
    end

    # Remove root if it no longer has two branches
    root = first(nodefilter(isroot, tree))
    if length(getchildren(tree, root)) < 2
        deletenode!(tree, root)
    end

    if !isempty(getleafinfo(tree))
        li = leafinfotype(typeof(tree))(Iterators.
                                        filter(line -> line[1] ∉ tipnames,
                                               getiterator(getleafinfo(tree))))
        setleafinfo!(tree, li)
    end

    return tipnames
end

"""
    keeptips!(tree::AbstractTree{OneTree}, tips)

Function to keep only the tips in a phylogenetic tree, `tree`, that are found in
the vector of tips or tip names, `tips`.
"""
function keeptips! end
@traitfn function keeptips!(tree::T,
                            tips::AbstractVector{N}) where
                  {N, RT,
                   T <: AbstractTree{OneTree, RT, N}; !MatchNodeType{T, N}}
    return keeptips!(tree, [getnode(tree, tip) for tip in tips])
end

@traitfn function keeptips!(tree::T,
                            tips::AbstractVector{N}) where
                  {N, RT, T <: AbstractTree{OneTree, RT}; MatchNodeType{T, N}}
    drop_tips = [tip for tip in getleaves(tree) if tip ∉ tips]
    return droptips!(tree, drop_tips)
end
