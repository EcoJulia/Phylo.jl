using Phylo
using Phylo.API
using Printf

abstract type OutputType end
abstract type NewickLike <: OutputType end
struct Newick <: NewickLike end
struct Nexus <: NewickLike end
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

function outputnode!(io::IO, tree::AbstractTree{TT, RT, String}, node, ::Newick,
                     ::Type{Nothing}) where {TT, RT}
    print(io, "\"", getnodename(tree, node), "\"")
    return nothing
end

function outputnode!(io::IO, tree::AbstractTree{TT, RT, <: Number}, node, ::Newick,
                     ::Type{Nothing}) where {TT, RT}
    print(io, getnodename(tree, node))
    return nothing
end

function outputnode!(io::IO, tree::AbstractTree, node, ::Newick, ::Type{<: Dict})
    print(io, "\"", getnodename(tree, node), "\"")
    nd = getnodedata(tree, node)
    if !isempty(nd)
        print(io, "[&")
        for (i, key) in enumerate(keys(nd))
            if i > 1
                print(io, ",")
            end
            value = nd[key]
            if value isa String
                print(io, key, "=\"", value, "\"")
            elseif value isa Number
                print(io, key, "=", value)
            elseif value isa Vector
                print(io, key, "={")
                for (j, elt) in enumerate(value)
                    if j > 1
                        print(io, ",")
                    end
                    if elt isa String
                        print(io, "\"", elt, "\"")
                    elseif elt isa Number
                        print(io, elt)
                    end
                end
                print(io, "}")
            end
        end
        print(io, "]")
    end
    return nothing
end

function outputnode(tree::AbstractTree, node, ot::OutputType)
    ios = IOBuffer()
    outputnode!(ios, tree, node, ot)
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
    outputnode!(io, tree, branch.in, CompactOutput(), Nothing)
    print(io, "] --> (")
    outputbranch!(io, tree, branch, CompactOutput(), Nothing)
    print(io, ") --> [")
    outputnode!(io, tree, branch.conns[1], CompactOutput(), Nothing)
    prinln(io, "]")
    return nothing
end

function outputbranch!(io::IO, tree::AbstractTree, branch, ::NewickLike, ::Type{Nothing})
    if haslength(tree, branch)
        print(io, ":", getlength(tree, branch))
    end
    return nothing
end

function outputbranch!(io::IO, tree::AbstractTree, branch, ::NewickLike, ::Type{<: Dict})
    if haslength(tree, branch)
        print(io, ":")
        bd = getbranchdata(tree, branch)
        if !isempty(bd)
            print(io, "[&")
            for (i, key) in enumerate(keys(bd))
                if i > 1
                    print(io, ",")
                end
                value = nd[key]
                if value isa String
                    print(io, key, "=\"", value, "\"")
                elseif value isa Number
                    print(io, key, "=", value)
                elseif value isa Vector
                    print(io, key, "={")
                    for (j, elt) in enumerate(value)
                        if j > 1
                            print(io, ",")
                        end
                        if elt isa String
                            print(io, "\"", elt, "\"")
                        elseif elt isa Number
                            print(io, elt)
                        end
                    end
                    print(io, key, "}")
                end
            end
            print(io, "]")
        end
        print(io, getlength(tree, branch))
    end
    return nothing
end

function outputbranch(tree::AbstractTree, branch, ot::OutputType)
    ios = IOBuffer()
    outputbranch!(ios, tree, branch, ot)
    return String(take!(ios))
end

function outputsubtree!(io::IO, tree::T, node, ot::Newick) where T <: AbstractTree{OneTree, OneRoot}
    if !isleaf(tree, node)
        print(io, "(")
        for (i, branch) in enumerate(getoutbounds(tree, node))
            if i > 1
                print(io, ",")
            end
            outputsubtree!(io, tree, dst(tree, branch), ot)
            outputbranch!(io, tree, branch, ot, branchdatatype(T))
        end
        print(io, ")")
    end
    outputnode!(io, tree, node, ot, nodedatatype(T))
    return nothing
end

function outputsubtree!(io::IO, tree::T, node, ot::Newick, exclude = []) where T <: AbstractTree{OneTree, Unrooted}
    cs = getconnections(tree, node, exclude)
    if !isempty(cs)
        print(io, "(")
        for (i, branch) in cs
            if i > 1
                print(io, ",")
            end
            outputsubtree!(io, tree, dst(tree, branch), ot, [branch])
            outputbranch!(io, tree, branch, ot, branchdatatype(T))
        end
        print(io, ")")
    end
    outputnode!(io, tree, node, ot, nodedatatype(T))
    return nothing
end

function outputtree!(io::IO, tree::AbstractTree{OneTree}, node, ot::Newick)
    outputsubtree!(io::IO, tree::AbstractTree{OneTree}, node, ot::Newick)
    print(io, ";")
    return nothing
end

outputtree!(io::IO, tree::AbstractTree{OneTree, OneRoot}, ot::Newick) =
    outputtree!(io, tree, getroot(tree), ot)

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
