# SPDX-License-Identifier: BSD-2-Clause

module PhyloRCallExt

if isdefined(Base, :get_extension)
    using RCall
    using RCall: protect, unprotect, rcall_p, RClass, isObject, isS4
    import RCall: rcopy
    import RCall: rcopytype
    import RCall: sexp
else
    using ..RCall
    using ..RCall: protect, unprotect, rcall_p, RClass, isObject, isS4
    import ..RCall: rcopy
    import ..RCall: rcopytype
    import ..RCall: sexp
end

using Phylo

function rcopy(::Type{T}, rt::Ptr{VecSxp}) where {T <: AbstractTree}
    if !isObject(rt) || isS4(rt) ||
       rcopy(String, rcall_p(:class, rt)) != "phylo"
        error("Object is not of S3 phylo class, aborting")
    end

    if !rcopy(Bool, rcall_p(Symbol("is.rooted"), rt))
        error("Cannot currently translate unrooted trees")
    end

    dict = rcopy(Dict{Symbol, Any}, rt)
    nodes = dict[:tip_label]
    tree = T(nodes)
    edges = dict[:edge]
    nnode = dict[:Nnode]
    lengths = dict[:edge_length]
    nontips = nnode
    append!(nodes, getnodename.(tree, createnodes!(tree, nontips)))

    for edge in Base.axes(edges, 1)
        createbranch!(tree,
                      nodes[edges[edge, 1]], nodes[edges[edge, 2]],
                      lengths[edge])
    end

    validate!(tree) || @warn "Tree does not internally validate"
    return tree
end

rcopytype(::Type{RClass{:phylo}}, s::Ptr{VecSxp}) = RootedTree

function sexp(tree::T) where {T <: AbstractTree}
    validate!(tree) || @warn "Tree does not internally validate"
    leafnames = getleafnames(tree)
    root = getnodename.(tree, getroots(tree))
    if (length(root) != 1)
        error("Can't currently translate tree with > 1 roots")
    end
    nontips = [getnodename(tree, node)
               for node in traversal(tree, preorder)
               if isinternal(tree, node)]
    tor = Dict{Symbol, Any}()
    tor[:Nnode] = length(nontips) + length(root)
    tor[Symbol("tip.label")] = leafnames
    nodes = copy(leafnames)
    push!(nodes, root[1])
    append!(nodes, nontips)
    bi = branchiter(tree)
    lengths = Vector{Float64}(undef, length(bi))
    edges = Matrix{Int32}(undef, length(lengths), 2)
    index = 1
    for branch in bi
        lengths[index] = getlength(tree, branch)
        edges[index, :] = indexin([getnodename(tree, src(tree, branch)),
                                      getnodename(tree, dst(tree, branch))], nodes)
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

end
