importall Phylo.API

"""
    NodeTree

Binary phylogenetic tree object with known leaves and per node data
"""
type NodeTree{LI <: AbstractInfo, ND} <: AbstractTree{String, Int}
    nodes::Dict{String, BinaryNode{Int}}
    branches::Dict{Int, Branch{String}}
    leafrecords::Dict{String, LI}
    noderecords::Dict{String, ND}
    rootheight::Nullable{Float64}
end

function NodeTree{LI, ND}(lt::NodeTree{LI, ND}; copyinfo=true, empty=true)
    validate(lt) || error("Tree to copy is not valid")
    leafnames = getleafnames(lt)
    # Leaf records may be conserved across trees, as could be invariant?
    leafrecords = copyinfo ? deepcopy(getleafrecords(lt)) : getleafrecords(lt)
    if empty # Empty out everything else
        nodes = Dict(map(leaf -> leaf => BinaryNode{Int}(), leafnames))
        branches = Dict{Int, Branch{String}}()
        noderecords = Dict(map(leaf -> leaf => ND(), leafnames))
    else # Make copies of everything
        nodes = deepcopy(nodes)
        noderecords = deepcopy(getnoderecords(lt))
        branches = deepcopy(getbranches(lt))
    end
    return NodeTree{LI, ND}(nodes, branches, leafrecords, noderecords,
                            lt.rootheight)
end

function NodeTree(leaves::AbstractVector{String};
                  rootheight::Nullable{Float64} = Nullable{Float64}(),
                  leaftype::Type = LeafInfo,
                  nodetype::Type = Void)
    leaftype <: AbstractInfo ||
        error("Leaf information structure is not an subtype of AbstractInfo")
    nodes = Dict(map(leaf -> leaf => BinaryNode{Int}(), leaves))
    leafrecords = Dict(map(leaf -> leaf => leaftype(), leaves))
    noderecords = Dict(map(leaf -> leaf => nodetype(), leaves))
    return NodeTree{leaftype, nodetype}(nodes, Dict{Int, Branch{String}}(),
                                        leafrecords, noderecords, rootheight)
end

function NodeTree(numleaves::Int;
                  rootheight::Nullable{Float64} = Nullable{Float64}(),
                  leaftype::Type = LeafInfo,
                  nodetype::Type = Void)
    leaftype <: AbstractInfo ||
        error("Leaf information structure is not an subtype of AbstractInfo")
    leaves = map(num -> "Leaf $num", 1:numleaves)
    nodes = Dict(map(leaf -> leaf => BinaryNode{Int}(), leaves))
    leafrecords = Dict(map(leaf -> leaf => leaftype(), leaves))
    noderecords = Dict(map(leaf -> leaf => nodetype(), leaves))
    return NodeTree{leaftype, nodetype}(nodes, Dict{Int, Branch{String}}(),
                                        leafrecords, noderecords, rootheight)
end

function _nodetype{LI, ND}(::Type{NodeTree{LI, ND}})
    return BinaryNode{Int}
end

function _branchtype{LI, ND}(::Type{NodeTree{LI, ND}})
    return Branch{String}
end

function _getnodes(nt::NodeTree)
    return nt.nodes
end

function _getbranches(nt::NodeTree)
    return nt.branches
end

function _getleafnames(nt::NodeTree)
    return keys(nt.leafrecords)
end

function _getleafrecords(nt::NodeTree)
    return nt.leafrecords
end

function _getnoderecords(nt::NodeTree)
    return nt.noderecords
end

function _addnode!{LI, NR}(tree::NodeTree{LI, NR}, nodename)
    _setnode!(tree, nodename, BinaryNode{Int}())
    setnoderecord!(tree, nodename, NR())
    return nodename
end

function _deletenode!(tree::NodeTree, nodename)
    node = getnode(tree, nodename)
    if _hasinbound(node)
        deletebranch!(tree, _getinbound(node))
    end
    for b in _getoutbounds(node)
        deletebranch!(tree, n)
    end
    delete!(_getnodes(tree), nodename)
    delete!(_getnoderecords(tree), nodename)    
    return nodename
end

function _validate(tree::NodeTree)
    if Set(NodeNameIterator(tree, isleaf)) != Set(getleafnames(tree))
        warn("Leaf records do not match actual leaves of tree")
        return false
    end
    
    if Set(keys(_getnoderecords(tree))) != Set(keys(getnodes(tree)))
        warn("Leaf records do not match node records of tree")
        return false
    end
    
    rootheight = hasrootheight(tree) ? getrootheight(tree) : NaN
    for leaf in getleafnames(tree)
        if hasheight(tree, leaf)
            if isnan(rootheight)
                rootheight = getheight(tree, leaf) - heighttoroot(tree, leaf)
            end
            if !(getheight(tree, leaf) - rootheight ≈
                 heighttoroot(tree, leaf))
                warn("Leaf height ($(getheight(tree, leaf))) for $leaf does not match branches")
                return false
            end
        end
    end
    return true
end

function _hasrootheight(tree::NodeTree)
    return !isnull(tree.rootheight)
end

function _getrootheight(tree::NodeTree)
    return get(tree.rootheight)
end

function _setrootheight!(tree::NodeTree, height::Float64)
    tree.rootheight = height
    return height
end

function _clearrootheight!(tree::NodeTree)
    tree.rootheight = Nullable{Float64}()
end

"""
    NamedTree

Binary phylogenetic tree object with known leaves
"""
const NamedTree = NodeTree{LeafInfo, Void}

NamedTree(leaves::AbstractVector{String}) = NodeTree(leaves)
NamedTree(numleaves::Int) = NodeTree(numleaves)






_getnodenames(tree::AbstractTree) = collect(keys(_getnodes(tree)))
_getbranchnames(tree::AbstractTree) = collect(keys(_getbranches(tree)))
#  - _hasnode()
_hasnode(tree::AbstractTree, label) = haskey(_getnodes(tree), label)
#  - _hasbranch()
_hasbranch(tree::AbstractTree, label) = haskey(_getbranches(tree), label)
#  - _addbranch!()
function _addbranch!(tree::AbstractTree, source, target,
                     length::Float64, label)
    # Add the new branch
    _setbranch!(tree, label, Branch(source, target, length))
    
    # Update the associated source and target nodes
    _addoutbound!(getnode(tree, source), label)
    _setinbound!(getnode(tree, target), label)
    
    # Return updated tree
    return label
end
#  - _deletebranch!()
function _deletebranch!(tree::AbstractTree, label)
    # Find the branch
    branch = _getbranch(tree, label)
    # Remove branch reference from its source node
    _deleteoutbound!(_getsource(branch), label)
    # Remove branch reference from its target node
    _deleteinbound!(_gettarget(branch), label)
    # Remove branch itself
    delete!(_getbranches(tree), branch)
    # Return the branch label
    return label
end

#  - _deletenode!()
function _deletenode!(tree::AbstractTree, label)
    node = _getnode(tree, label)
    if _hasinbound(node)
        _deletebranch!(tree, _getinbound(node))
    end
    for b in _getoutbounds(node)
        _deletebranch!(tree, b)
    end
    delete!(_getnodes(tree), label)
    return label
end


#  - _getleafrecords()
function _getleafrecords(tree::AbstractTree)
    return Dict(map(leaf -> leaf=>nothing,
                    findleaves(tree) ∪ findunattacheds(tree)))
end
#  - _getleafnames()
function _getleafnames(tree::AbstractTree)
    return keys(Dict(map(leaf -> leaf=>nothing,
                         findleaves(tree) ∪ findunattacheds(tree))))
end
#  - _getleafrecord()
function _getleafrecord(tree::AbstractTree, label)
    return _getleafrecords(tree)[label]
end
#  - _setleafrecord()
function _setleafrecord!(tree::AbstractTree, label, value)
    _getleafrecords(tree)[label] = value
    return value
end
#  - _getnoderecords()
function _getnoderecords(tree::AbstractTree)
    return Dict(map(node -> node=>nothing, keys(_getnodes(tree))))
end
#  - _getnoderecord()
function _getnoderecord(tree::AbstractTree, label)
    return _getnoderecords(tree)[label]
end
#  - _setnoderecord!()
function _setnoderecord!(tree::AbstractTree, label, value)
    _getnoderecords(tree)[label] = value
    return value
end


"""
    clearrootheight(::AbstractTree)

Clears the tree's root height record.
"""
function clearrootheight!(tree::NodeTree)
    _clearrootheight!(tree)
end

function _getnode(tree::NodeTree, label)
    return _getnodes(tree)[label]
end

function _getbranch(tree::NodeTree, label)
    return _getbranches(tree)[label]
end

function _setnode!(tree::NodeTree, label, node)
    return _getnodes(tree)[label] = node
end

function _setbranch!(tree::NodeTree, label, branch)
    _hasbranch(tree, label) &&
        error("Branch $label already exists")
    return _getbranches(tree)[label] = branch
end

# Interfaces to an info
# ---------------------

function _hasheight(tree::NodeTree, label)
    return _hasheight(getleafrecord(tree, label))
end

function _getheight(tree::NodeTree, label)
    return _getheight(getleafrecord(tree, label))
end

function _setheight!(tree::NodeTree, label, height::Float64)
    ai = getleafrecord(tree, label)
    _setheight!(ai, height)
    setleafrecord!(tree, label, ai)
    return height
end
