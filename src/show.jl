using Phylo
using Phylo.API
using Compat
using Compat.Printf

function show(io::IO, object::AbstractNode, n::String = "")
    node = "node"
    if !isempty(n)
        node *= " $n"
    end
    if !_hasinbound(object)
        if _outdegree(object) > 0
            blank = repeat(" ", length(" [root $node]") + (isempty(n) ? 0 : 1))
            for (i, bn) in zip(1:_outdegree(object), _getoutbounds(object))
                b = typeof(bn) <: Number ? "$bn" : "\"$bn\""
                if _outdegree(object) == 1
                    print(io, "[root $node]-->[branch $b]")
                elseif get(io, :compact, false)
                    if i == 1
                        print(io, "[root $node]-->[branches $b")
                    elseif i < _outdegree(object)
                        print(io, ", $b")
                    else
                        print(io, " and $b]")
                    end
                else # multiline view
                    if i == 1
                        println(io, "[root $node]-->[branch $b]")
                    elseif i < _outdegree(object)
                        println(io, "$blank-->[branch $b]")
                    else
                        print(io, "$blank-->[branch $b]")
                    end
                end
            end
        else
            print(io, "[unattached $node]")
        end
    else # hasinbound
        inb = typeof(_getinbound(object)) <: Number ?
            "$(_getinbound(object))" : "\"$(_getinbound(object))\""
        if _outdegree(object) == 0
            print(io, "[branch $inb]-->[leaf $node]")
        elseif _hasinbound(object)
            blank = repeat(" ",
                           length(" [branch $inb]-->[internal $node]") +
                           (isempty(n) ? 0 : 1))
            for (i, bn) in zip(1:_outdegree(object), _getoutbounds(object))
                b = typeof(bn) <: Number ? "$bn" : "\"$bn\""
                if _outdegree(object) == 1
                    print(io, "[branch $inb]-->[internal $node]-->[branch $b]")
                elseif get(io, :compact, false)
                    if i == 1
                        print(io, "[branch $inb]-->[internal $node]-->" *
                              "[branches $b")
                    elseif i < _outdegree(object)
                        print(io, ", $b")
                    else
                        print(io, " and $b]")
                    end
                else # multiline view
                    if i == 1
                        println(io, "[branch $inb]-->[internal $node]-->" *
                                "[branch $b]")
                    elseif i < _outdegree(object)
                        println(io, "$blank-->[branch $b]")
                    else
                        print(io, "$blank-->[branch $b]")
                    end
                end
            end
        end
    end
end

function show(io::IO, p::Pair{NT, N}) where {N <: AbstractNode, NT}
    n = NT <: Number ? "$(p[1])" : "\"$(p[1])\""
    show(io, p[2], "$n")
end

function show(io::IO, object::B) where {B <: AbstractBranch}
    NT = typeof(_src(object))
    source = NT <: Number ? "$(_src(object))" : "\"$(_src(object))\""
    destination = NT <: Number ? "$(_dst(object))" : "\"$(_dst(object))\""
    println(io, "[node $source]-->[$(_getlength(object)) length branch]-->" *
            "[node $destination]")
end

function show(io::IO, p::Pair{BT, B}) where {BT, B <: AbstractBranch}
    NT = typeof(_src(p[2]))
    source = NT <: Number ? "$(_src(p[2]))" : "\"$(_src(p[2]))\""
    destination = NT <: Number ? "$(_dst(p[2]))" : "\"$(_dst(p[2]))\""
    branch = BT <: Number ? "$(p[1])" : "\"$(p[1])\""
    println(io, "[node $source]-->[$(_getlength(p[2])) length branch $branch]-->[node $destination]")
end

function showsimple(io::IO, object::TreeSet)
    @printf(io, "TreeSet with %d trees, each with %d tips.\n",
            ntrees(object), nleaves(object))
    tn = collect(treenameiter(object))
    if length(tn) < 10
        println(io, "Tree names are " *
                join(tn, ", ", " and "))
    else
        println(io, "Tree names are " *
                join(tn[1:5], ", ") *
                " ... $(length(tn) - 6) missing ... " *
                "$(tn[end])")
    end
end

function show(io::IO, object::TreeSet)
    showsimple(io, object)
    for name in treenameiter(object)
        @printf(io, "\n%s: ", name)
        showsimple(io, object[name])
    end
end

function showsimple(io::IO, object::TREE) where TREE <: AbstractBranchTree
    println(io, "$TREE with $(nleaves(object)) tips, " *
            "$(length(_getnodes(object))) nodes and " *
            "$(length(_getbranches(object))) branches.")
    ln = getleafnames(object)
    if length(ln) < 10
        println(io, "Leaf names are " *
                join(ln, ", ", " and "))
    else
        println(io, "Leaf names are " *
                join(ln[1:5], ", ") *
                ", ... [$(length(ln) - 6) omitted] ... and " *
                "$(ln[end])")
    end
end

function show(io::IO, object::TREE) where TREE <: AbstractBranchTree
    showsimple(io, object)
    if !get(io, :compact, true)
        println(io, "Nodes:")
        println(io, _getnodes(object))
        println(io, "Branches:")
        println(io, _getbranches(object))
    end
end

function show(io::IO, object::BinaryTree{LI, ND}) where {LI, ND}
    showsimple(io, object)
    if !get(io, :compact, true)
        println(io, "Nodes:")
        println(io, _getnodes(object))
        println(io, "Branches:")
        println(io, _getbranches(object))
        if ND !== Nothing
            println(io, "NodeData:")
            println(io, Dict(map(nodename ->
                                nodename => getnoderecord(object, nodename),
                                nodenameiter(object))))
        end
    end
end

function show(io::IO, object::PolytomousTree{LI, ND}) where {LI, ND}
    showsimple(io, object)
    if !get(io, :compact, true)
        println(io, "Nodes:")
        println(io, _getnodes(object))
        println(io, "Branches:")
        println(io, _getbranches(object))
        if ND !== Nothing
            println(io, "NodeData:")
            println(io, Dict(map(nodename ->
                                nodename => getnoderecord(object, nodename),
                                nodenameiter(object))))
        end
    end
end
