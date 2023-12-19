using Phylo
using Phylo.API
using Printf

abstract type OutputType end
struct CompactOutput <: OutputType end
struct StandardOutput <: OutputType end

function outputnode!(io::IO, tree::TREE, node, ::CompactOutput,
                     _ = nodedatatype(TREE)) where {TT, TREE <: AbstractTree{TT, Unrooted, String}}
    d = degree(tree, node)
    opening = d == 0 ? "unattached" : (d > 1 ? "node" : "leaf")
    print(io, "$opening \"$(getnodename(tree, node))\"")
    return nothing
end

function outputnode!(io::IO, tree::TREE, node, ::CompactOutput,
                     _ = nodedatatype(TREE)) where {TT, TREE <: AbstractTree{TT, Unrooted, <: Number}}
    d = degree(tree, node)
    opening = d == 0 ? "unattached" : (d > 1 ? "node" : "leaf")
    print(io, "$opening $(getnodename(tree, node))")
    return nothing
end

function outputnode!(io::IO, tree::TREE, node, ::CompactOutput,
                     _ = nodedatatype(TREE)) where {TT, TREE <: AbstractTree{TT, <: Rooted, String}}
    od = outdegree(tree, node)
    opening = hasinbound(tree, node) ?
        (od > 0 ? "root node" : "unattached node") :
        (od > 0 ? "internal node" : "leaf node")
    print(io, "$opening \"$(getnodename(tree, node))\"")
    return nothing
end

function outputnode!(io::IO, tree::TREE, node, ::CompactOutput,
                     _ = nodedatatype(TREE)) where {TT, TREE <: AbstractTree{TT, <: Rooted, <: Number}}
    od = outdegree(tree, node)
    opening = hasinbound(tree, node) ?
        (od > 0 ? "root node" : "unattached node") :
        (od > 0 ? "internal node" : "leaf node")
    print(io, "$opening $(getnodename(tree, node))")
    return nothing
end

function outputnode!(io::IO, tree::AbstractTree, node, ::StandardOutput, ot)
    ios = IOBuffer()
    od = outdegree(tree, node)
    if hasinbound(tree, node)
        print(ios, "(")
        outputbranch!(ios, tree, getinbound(tree, node), CompactOutput(), Nothing)
        print(ios, ") --> ")
    end
    print(ios, "[")
    outputnode!(ios, tree, node, CompactOutput(), ot)
    print(ios, "]")
    opening = String(take!(ios))

    blank = repeat(" ", length(opening))
    for (i, b) in enumerate(getoutbounds(tree, node))
        print(io, i == 1 ? "$opening --> (" :  "$blank --> (")
        outputbranch!(io, tree, b, CompactOutput(), Nothing)
        println(io, ")")
    end
    return nothing
end

function outputnode(tree::TREE, node, ot::OutputType) where TREE <: AbstractTree
    ios = IOBuffer()
    outputnode!(ios, tree, node, ot, nodedatatype(TREE))
    return String(take!(ios))
end

function outputbranch!(io::IO, tree::TREE, branch, ::CompactOutput,
                       _ = branchdatatype(TREE)) where TREE <: AbstractTree
    print(io, "branch $(getbranchname(tree, branch))")
    if haslength(tree, branch)
        print(io, ": length $(getlength(tree, branch))")
    end
    return nothing
end

function outputbranch!(io::IO, tree::TREE, branch, ::StandardOutput,
                       _ = branchdatatype(TREE)) where {TT, TREE <: AbstractTree{TT, <: Rooted}}
    print(io, "[")
    outputnode!(io, tree, src(tree, branch), CompactOutput(), Nothing)
    print(io, "] --> (")
    outputbranch!(io, tree, branch, CompactOutput(), Nothing)
    print(io, ") --> [")
    outputnode!(io, tree, dst(tree, branch), CompactOutput(), Nothing)
    println(io, "]")
    return nothing
end

function outputbranch(tree::AbstractTree, branch, ot::OutputType)
    ios = IOBuffer()
    outputbranch!(ios, tree, branch, ot)
    return String(take!(ios))
end

function outputtree!(io::IO, tree::TREE, ::CompactOutput) where TREE <: AbstractTree{OneTree}
    print(io, "$TREE with $(nleaves(tree)) tips and $(nroots(tree)) $(nroots(tree) == 1 ? "root" : "roots"). ")
    ln = getleafnames(tree)
    if length(ln) < 10
        print(io, "Leaf names are " * join(ln, ", ", " and "))
    else
        print(io, "Leaf names are " * join(ln[1:5], ", ") *
              ", ... [$(length(ln) - 6) omitted] ... and $(ln[end])")
    end
end

function outputtree!(io::IO, tree::TREE, ::StandardOutput) where TREE <: AbstractTree{OneTree}
    outputtree!(io, tree, CompactOutput())
    n = nnodes(tree)
    ns = collect(getnodes(tree))
    print(io, "\n\n$n nodes: [")
    if n < 10
        println(io, join(ns, ", ", " and ") * "]")
    else
        println(io, join(ns[1:5], ", ") *
                " ... $(n - 6) missing ... $(ns[end])]")
    end

    b = nbranches(tree)
    bs = collect(getbranches(tree))
    print(io, "\n$b branches: [")
    if n < 10
        println(io, join(bs, ", ", " and ") * "]")
    else
        println(io, join(bs[1:5], ", ") *
                " ... $(b - 6) missing ... $(collect(bs)[end])]")
    end

    if nodedatatype(TREE) ≢ Nothing
        print(io, "\nNode records: ")
        nn = getnodenames(tree)
        if nnodes(tree) == 1
            println(io, Dict(nn[1] => getnodedata(tree, nn[1])))
        elseif nnodes(tree) == 2
            println(io, nn[1] => getnodedata(tree, nn[1]), nn[2] => getnodedata(tree, nn[2]))
        elseif nnodes(tree) > 2
            println(io, nn[1] => getnodedata(tree, nn[1]), " ... ", nn[end] => getnodedata(tree, nn[end]))
        end
    end
end

function outputtree!(io::IO, treeset::TREE, ::CompactOutput) where TREE <: AbstractTree{ManyTrees}
    println(io, "$TREE with $(ntrees(treeset)) tree(s), each with $(nleaves(treeset)) tips.")
    tn = collect(gettreenames(treeset))
    if length(tn) == 1
        print(io, "Tree name is $(tn[1])")
    elseif length(tn) < 10
        print(io, "Tree names are " * join(tn, ", ", " and "))
    else
        print(io, "Tree names are " * join(tn[1:5], ", ") *
              " ... $(length(tn) - 6) missing ... $(tn[end])")
    end
    println(io, ". $(nnodes(treeset)) nodes and $(nbranches(treeset)) branches.")
    return nothing
end

function outputtree!(io::IO, treeset::TREE, ::StandardOutput) where TREE <: AbstractTree{ManyTrees}
    outputtree!(io, treeset, CompactOutput())
    nt = ntrees(treeset)
    if nt < 10
        for treename in gettreenames(treeset)
            print(io, "\n$treename: ")
            outputtree!(io, treeset[treename], CompactOutput())
        end
    else
        for (index, treename) in enumerate(gettreenames(treeset))
            if index ≤ 5 || index == nt
                print(io, "\n$treename: ")
                outputtree!(io, treeset[treename], CompactOutput())
            elseif index == 6
                println(io, "\n[$(nt - 6) trees omitted]\n")
            end
        end
    end
    return nothing
end

function outputtree(tree::AbstractTree, ot::OutputType)
    ios = IOBuffer()
    outputtree!(ios, tree, ot)
    return String(take!(ios))
end

import Base: show
function show(io::IO, tree::TREE) where TREE <: AbstractTree{OneTree}
    if get(io, :compact, false)
        outputtree!(io, tree, CompactOutput())
    else
        outputtree!(io, tree, StandardOutput())
    end
end

function show(io::IO, treeset::TREE) where TREE <: AbstractTree{ManyTrees}
    if get(io, :compact, false)
        outputtree!(io, treeset, CompactOutput())
    else
        outputtree!(io, treeset, StandardOutput())
    end
end
