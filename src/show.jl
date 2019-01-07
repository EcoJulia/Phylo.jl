using Phylo
using Phylo.API
using Compat
using Compat.Printf

function show(io::IO, object::Tuple{<: AbstractTree, <: AbstractNode},
              n::String = "")
    node = "node"
    if !isempty(n)
        node *= " $n"
    end
    if !_hasinbound(object[1], object[2])
        if _outdegree(object[1], object[2]) > 0
            blank = repeat(" ", length(" [root $node]") + (isempty(n) ? 0 : 1))
            for (i, bn) in zip(Base.OneTo(_outdegree(object[1], object[2])),
                               _getoutbounds(object[1], object[2]))
                b = typeof(bn) <: Number ? "$bn" : "\"$bn\""
                if _outdegree(object[1], object[2]) == 1
                    print(io, "[root $node]-->[branch $b]")
                elseif get(io, :compact, false)
                    if i == 1
                        print(io, "[root $node]-->[branches $b")
                    elseif i < _outdegree(object[1], object[2])
                        print(io, ", $b")
                    else
                        print(io, " and $b]")
                    end
                else # multiline view
                    if i == 1
                        println(io, "[root $node]-->[branch $b]")
                    elseif i < _outdegree(object[1], object[2])
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
        inb = typeof(_getinbound(object[1], object[2])) <: Number ?
            "$(_getinbound(object[1], object[2]))" :
            "\"$(_getinbound(object[1], object[2]))\""
        if _outdegree(object[1], object[2]) == 0
            print(io, "[branch $inb]-->[leaf $node]")
        elseif _hasinbound(object[1], object[2])
            blank = repeat(" ",
                           length(" [branch $inb]-->[internal $node]") +
                           (isempty(n) ? 0 : 1))
            for (i, bn) in zip(Base.OneTo(_outdegree(object[1], object[2])),
                               _getoutbounds(object[1], object[2]))
                b = typeof(bn) <: Number ? "$bn" : "\"$bn\""
                if _outdegree(object[1], object[2]) == 1
                    print(io, "[branch $inb]-->[internal $node]-->[branch $b]")
                elseif get(io, :compact, false)
                    if i == 1
                        print(io, "[branch $inb]-->[internal $node]-->" *
                              "[branches $b")
                    elseif i < _outdegree(object[1], object[2])
                        print(io, ", $b")
                    else
                        print(io, " and $b]")
                    end
                else # multiline view
                    if i == 1
                        println(io, "[branch $inb]-->[internal $node]-->" *
                                "[branch $b]")
                    elseif i < _outdegree(object[1], object[2])
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

function show(io::IO, object::Tuple{<: AbstractTree{OneTree, RT, NL},
                                    <: AbstractBranch{RT, NL}}) where {RT, NL}
    source = NL <: Number ? "$(_src(object[1], object[2]))" :
        "\"$(_src(object[1], object[2]))\""
    destination = NL <: Number ? "$(_dst(object[1], object[2]))" :
        "\"$(_dst(object[1], object[2]))\""
    println(io, "[node $source]-->[$(_getlength(object[1], object[2])) " *
        "length branch]--> [node $destination]")
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
    for name in _gettreenames(object)
        @printf(io, "\n%s: ", name)
        showsimple(io, object[name])
    end
end

function showsimple(io::IO, object::TREE) where {TREE <: AbstractTree}
    println(io, "$TREE with $(_nleaves(object)) tips, " *
            "$(_nnodes(object)) nodes and " *
            "$(_nbranches(object)) branches.")
    ln = getleafnames(object)
    if length(ln) < 10
        println(io, "Leaf names are " * join(ln, ", ", " and "))
    else
        println(io, "Leaf names are " * join(ln[1:5], ", ") *
                ", ... [$(length(ln) - 6) omitted] ... and $(ln[end])")
    end
end

function show(io::IO, object::AbstractTree)
    showsimple(io, object)
    if !get(io, :compact, true)
        println(io, "Nodes:")
        println(io, [(object, node) for node in _getnodes(object)])
        println(io, "Branches:")
        println(io, [(object, branch) for branch in _getbranches(object)])
        if _noderecordtype(TREE) !== Nothing
            println(io, "Node records:")
            println(io, Dict(name => _getnoderecord(object, name)
                             for name in _getnodenames(object)))
        end
    end
end
