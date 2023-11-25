using Phylo
using Phylo.API
using Printf

function show(io::IO, object::NamedTuple{(:tree, :node), Tuple{T, N}}) where 
              {RT <: Rooted, NL, T <: AbstractTree{OneTree, RT, NL}, N}
    node = getnodename(object.tree, object.node)
    od = outdegree(object.tree, object.node)
    if !hasinbound(object.tree, object.node)
        if od > 0
            blank = repeat(" ", length("[root $node]"))
            for (i, bn) in enumerate(getbranchname.(Ref(object.tree), getoutbounds(object.tree, object.node)))
                b = typeof(bn) <: Number ? "branch $bn" : "\"$bn\""
                if od == 1
                    print(io, "[root $node] --> (branch $b)")
                elseif get(io, :compact, false)
                    if i == 1
                        print(io, "[root $node] --> (branches $b,")
                    elseif i < od
                        print(io, "$b,")
                    else
                        print(io, "and $b)")
                    end
                else # multiline view
                    if i == 1
                        println(io, "[root $node] --> (branch $b)")
                    elseif i < od
                        println(io, "$blank --> (branch $b)")
                    else
                        print(io, "$blank --> (branch $b)")
                    end
                end
            end
        else
            print(io, "[unattached $node]")
        end
    else # hasinbound
        bn = getbranchname(object.tree, getinbound(object.tree, object.node))
        inb = typeof(bn) <: Number ? "$bn" : "\"bn\""
        if od == 0
            print(io, "(branch $inb) --> [leaf $node]")
        else
            blank = repeat(" ", length("(branch $inb) --> [internal $node]"))
            for (i, bn) in enumerate(getbranchname.(Ref(object.tree), getoutbounds(object.tree, object.node)))
                b = typeof(bn) <: Number ? "$bn" : "\"$bn\""
                if od == 1
                    print(io, "(branch $inb) --> [internal $node] --> (branch $b)")
                elseif get(io, :compact, false)
                    if i == 1
                        print(io, "(branch $inb) --> [internal $node] --> (branches $b,")
                    elseif i < od
                        print(io, "$b,")
                    else
                        print(io, "and $b)")
                    end
                else # multiline view
                    if i == 1
                        println(io, "(branch $inb) --> [internal $node] --> (branch $b)")
                    elseif i < od
                        println(io, "$blank --> (branch $b)")
                    else
                        print(io, "$blank --> (branch $b)")
                    end
                end
            end
        end
    end
end

function show(io::IO, object::NamedTuple{(:tree, :branch), Tuple{T, B}}) where 
    {RT <: Rooted, NL, T <: AbstractTree{OneTree, RT, NL}, B}
    source = NL <: Number ? "node $(getnodename(object.tree, src(object.tree, object.branch)))" :
        "\"$(getnodename(object.tree, src(object.tree, object.branch)))\""
    destination = NL <: Number ? "node $(getnodename(object.tree, dst(object.tree, object.branch)))" :
        "\"$(getnodename(object.tree, dst(object.tree, object.branch)))\""
    if haslength(object.tree, object.branch)
        print(io, "[$source] --> [$(getlength(object.tree, object.branch)) " *
              "length branch $(getbranchname(object.tree, object.branch))] --> [$destination]")
    else
        print(io, "[$source] --> (branch $(getbranchname(object.tree, object.branch))) --> [$destination]")
    end
end

function showsimple(io::IO, object::TreeSet)
    @printf(io, "TreeSet with %d trees, each with %d tips.\n",
            ntrees(object), nleaves(object))
    tn = collect(gettreenames(object))
    if length(tn) < 10
        println(io, "Tree names are " * join(tn, ", ", " and "))
    else
        println(io, "Tree names are " * join(tn[1:5], ", ") *
                " ... $(length(tn) - 6) missing ... $(tn[end])")
    end
end

function show(io::IO, object::TreeSet)
    showsimple(io, object)
    tn =  sort(collect(gettreenames(object)))
    index = 0
    for name in tn
        index += 1
        if index ≤ 5 || index == length(tn)
            @printf(io, "\n%s: ", name)
            showsimple(io, object[name])
        elseif index == 6
            println("\n[$(length(tn)-6) trees omitted]\n")
        end
    end
end

function showsimple(io::IO, object::TREE) where TREE <: AbstractTree
    print(io, "$TREE with $(nleaves(object)) tips and $(nroots(object)) roots. ")
    ln = getleafnames(object)
    if length(ln) < 10
        print(io, "Leaf names are " * join(ln, ", ", " and "))
    else
        print(io, "Leaf names are " * join(ln[1:5], ", ") *
              ", ... [$(length(ln) - 6) omitted] ... and $(ln[end])")
    end
end

function show(io::IO, object::AbstractTree)
    showsimple(io, object)
    if get(io, :compact, false)
        println("$(nnodes(object)) nodes and $(nbranches(object)) branches.")
    else
        println(io)
        println(io, "$(nnodes(object)) nodes:")
        for node in getnodes(object)
            println(io, (tree=object, node=node))
        end
        println(io, "$(nbranches(object)) branches:")
        for branch in getbranches(object)
            println(io, (tree=object, branch=branch))
        end
        if nodedatatype(typeof(object)) ≢ Nothing
            println(io, "Node records:")
            println(io, Dict(name => getnodedata(object, name)
                             for name in getnodenames(object)))
        end
    end
end
