using Phylo
using Phylo.API
using Compat

function show(io::IO, object::AbstractNode, n::String = "")
    node = "node"
    if !isempty(n)
        node *= " $n"
    end
    if !_hasinbound(object)
        if _outdegree(object) > 0
            blank = repeat(" ", length("[root $node]") + (isempty(n) ? 0 : 1))
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
                        print(io, "[root $node]-->[branch $b]\n")
                    elseif i < _outdegree(object)
                        print(io, "$blank-->[branch $b]\n")
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
                           length("[branch $inb]-->[internal $node]") +
                           (isempty(n) ? 0 : 1))
            for (i, bn) in zip(1:_outdegree(object), _getoutbounds(object))
                b = typeof(bn) <: Number ? "$bn" : "\"$bn\""
                if _outdegree(object) == 1
                    print(io, "[branch $inb]-->[internal $node]-->[branch $b]")
                elseif get(io, :compact, false)
                    if i == 1
                        print(io, "[branch $inb]-->[internal $node]-->[branches $b")
                    elseif i < _outdegree(object)
                        print(io, ", $b")
                    else
                        print(io, " and $b]")
                    end
                else # multiline view
                    if i == 1
                        print(io, "[branch $inb]-->[internal $node]-->[branch $b]\n")
                    elseif i < _outdegree(object)
                        print(io, "$blank-->[branch $b]\n")
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
    print(io, "[node $source]-->[$(_getlength(object)) length branch]-->[node $destination]")   
end  

function show(io::IO, p::Pair{BT, B}) where {BT, B <: AbstractBranch}
    NT = typeof(_src(p[2]))
    source = NT <: Number ? "$(_src(p[2]))" : "\"$(_src(p[2]))\""
    destination = NT <: Number ? "$(_dst(p[2]))" : "\"$(_dst(p[2]))\""
    branch = BT <: Number ? "$(p[1])" : "\"$(p[1])\""
    print(io, "[node $source]-->[$(_getlength(p[2])) length branch $branch]-->[node $destination]")
end

function show(io::IO, object::AbstractTree)
    print(io, "Phylogenetic tree with $(length(_getnodes(object))) nodes and $(length(_getbranches(object))) branches")
function show(io::IO, object::TreeSet)
    @printf(io, "PhyloSet with %d trees\n", ntrees(object))
end

function showall(io::IO, object::TreeSet)
    show(io, object)
    println(io, "Trees:")
    foreach(t -> println(io, t), treenameiter(object))
end

function showall(io::IO, object::AbstractTree)
    print(io, object)
    print(io, "Nodes:\n") 
    print(io, _getnodes(object))
    print(io, "Branches:\n") 
    print(io, _getbranches(object))
end

function show(io::IO, object::TREE) where TREE <: AbstractBranchTree
    if get(io, :compact, false)
        print(io, "$TREE phylogenetic tree with $(length(_getnodes(object))) nodes ($(length(getleafnames(object))) leaves) and $(length(_getbranches(object))) branches")
    else
        print(io, "$TREE phylogenetic tree with $(length(_getnodes(object))) nodes and $(length(_getbranches(object))) branches\n") 
        print(io, "Leaf names:\n")
        print(io, getleafnames(object))
    end
end

function showall(io::IO, object::BinaryTree{LI, ND}) where {LI, ND}
    print(io, object)
    print(io, "\nNodes:\n") 
    print(io, _getnodes(object))
    print(io, "\nBranches:\n") 
    print(io, _getbranches(object))
    print(io, "\nNodeData:\n") 
    print(io, Dict(map(nodename -> nodename => getnoderecord(object, nodename),
                       nodenameiter(object))))
end

function showall(io::IO, object::BinaryTree{LI, Nothing}) where LI
    print(io, object)
    print(io, "\nNodes:\n") 
    print(io, _getnodes(object))
    print(io, "\nBranches:\n") 
    print(io, _getbranches(object))
end

function showall(io::IO, object::PolytomousTree{LI, ND}) where {LI, ND}
    print(io, object)
    print(io, "\nNodes:\n") 
    print(io, _getnodes(object))
    print(io, "\nBranches:\n") 
    print(io, _getbranches(object))
    print(io, "\nNodeData:\n") 
    print(io, Dict(map(nodename -> nodename => getnoderecord(object, nodename),
                       nodenameiter(object))))
end

function showall(io::IO, object::PolytomousTree{LI, Nothing}) where LI
    print(io, object)
    print(io, "\nNodes:\n") 
    print(io, _getnodes(object))
    print(io, "\nBranches:\n") 
    print(io, _getbranches(object))
end
