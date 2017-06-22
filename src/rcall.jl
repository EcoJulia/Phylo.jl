using Phylo
using RCall
using RCall: protect, unprotect, rcall_p, RClass, isObject, isS4

import RCall.rcopy

function rcopy{T <: AbstractTree}(::Type{T}, rt::Ptr{VecSxp})
    if !isObject(rt) || isS4(rt) || rcopy(rcall_p(:class, rt)) != "phylo"
        error("Object is not of S3 phylo class, aborting")
    end

    if !rcopy(rcall_p(Symbol("is.rooted"), rt))
        error("Cannot currently translate unrooted trees")
    end

    dict = rcopy(Dict{Symbol, Any}, rt)
    nodes = dict[Symbol("tip.label")]
    tree = NamedTree(nodes)
    edges = dict[:edge]
    nnode = dict[:Nnode]
    lengths = dict[Symbol("edge.length")]
    nontips = nnode
    append!(nodes, addnodes!(tree, nontips))
    
    for edge in 1:size(edges, 1)
        addbranch!(tree,
                   nodes[edges[edge, 1]], nodes[edges[edge, 2]],
                   lengths[edge])
    end

    validate(tree) || warn("Tree does not internally validate")
    return tree
end

import RCall.rcopytype

rcopytype(::Type{RClass{:phylo}}, s::Ptr{VecSxp}) = NamedTree

import RCall.sexp

function sexp(tree::AbstractTree)
    validate(tree) || warn("Tree does not internally validate")

    tipnames = collect(NodeNameIterator(tree, isleaf))
    root = collect(NodeNameIterator(tree, isroot))
    if (length(root) != 1)
        error("Can't currently translate tree with > 1 roots")
    end
    nontips = collect(NodeNameIterator(tree, isinternal))
    tor = Dict{Symbol, Any}()
    tor[:Nnode] = length(nontips) + length(root)
    tor[Symbol("tip.label")] = tipnames
    nodes = copy(tipnames)
    push!(nodes, root[1])
    append!(nodes, nontips)
    bi = BranchIterator(tree)
    lengths = Vector{Float64}(length(bi))
    edges = Matrix{Int32}(length(lengths), 2)
    index = 1
    for branch in bi
        lengths[index] = getlength(branch)
        edges[index, :] = indexin([src(branch), dst(branch)], nodes)
        index += 1
    end
    tor[:edge] = edges
    tor[Symbol("edge.length")] = lengths
    sobj = protect(sexp(tor))
    #setattrib!(sobj, :order, sexp("cladewise"))
    setclass!(sobj, sexp("phylo"))
    unprotect(1)
    return sobj
end
