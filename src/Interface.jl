using Phylo.API
using Compat: mapreduce

getnodes(tree::AbstractTree) = _getnodes(tree)
getbranches(tree::AbstractTree) = _getbranches(tree)

function ntrees(tree::AbstractTree)
    return _ntrees(tree)
end

# AbstractTree methods
"""
    nodetype(tree::AbstractTree)

Returns type of nodes in a tree.
"""
nodetype(tree::AbstractTree) = _nodetype(tree)

"""
    branchtype(tree::AbstractTree)

Returns type of branches in a tree.
"""
branchtype(tree::AbstractTree) = _branchtype(tree)

"""
    nodenametype(::AbstractTree)

Returns type of node names.
"""
nodenametype(::T) where {NL, BL, T <: AbstractTree{NL, BL}} = NL

"""
    branchnametype(::AbstractTree)

Returns type of branch names.
"""
branchnametype(::T) where {NL, BL, T <: AbstractTree{NL, BL}} = BL

"""
    addbranch!(tree::AbstractTree, source, destination[, length::Float64];
               branchname = _newbranchlabel(tree))

Add a branch from `source` to `destination` on `tree`.
"""
function addbranch!(tree::AbstractTree, source, destination, length::Float64 = NaN;
                    branchname = _newbranchlabel(tree))
    _hasnode(tree, source) ||
    error("Tree does not have an available source node called $source")
    hasoutboundspace(tree, source) ||
    error("$source already has maximum number of outbound connections ($(outdegree(tree, source)))")
    _hasnode(tree, destination) ||
        error("Tree does not have a destination node called $destination")
    !hasinbound(tree, destination) ||
            error("Tree does not have an available destination node called $destination")
    destination != source || error("Branch must connect different nodes")
    _hasbranch(tree, branchname) &&
        error("Tree already has a branch called $branchname")

    return _addbranch!(tree, source, destination, length, branchname)
end

"""
    deletebranch!(tree::AbstractTree, branchname)

Delete the branch `branchname` from `tree`.
"""
function deletebranch!(tree::AbstractTree, branchname)
    _hasbranch(tree, branchname) ||
        error("Tree does not have a branch called $branchname")
    return _deletebranch!(tree, branchname)
end

"""
    branch!(tree::AbstractTree, source[, length])
    branch!(tree::AbstractTree, source[, length]; destination)
    branch!(tree::AbstractTree, source[, length]; destination, branchname)

Branch from a source node `source` and create a destination node `destination`.
"""
function branch!(tree::AbstractTree, source, length::Float64 = NaN;
                 destination = _newnodelabel(tree),
                 branchname = _newbranchlabel(tree))
    _hasnode(tree, source) ||
        error("Node $source not present in tree")
    !_hasnode(tree, destination) ||
        error("Node $destination already present in tree")
    _hasoutboundspace(_getnode(tree, source)) ||
        error("Node $source has no space to add branches")

    return _branch!(tree, source, length, destination, branchname)
end

"""
    addnode!(tree::AbstractTree)
    addnode!(tree::AbstractTree, nodename)


"""
function addnode!(tree::AbstractTree, nodename = _newnodelabel(tree))
    !_hasnode(tree, nodename) ||
        error("Node $nodename already present in tree")
    return _addnode!(tree, nodename)
end

"""
    addnodes!(tree::AbstractTree, nodenames::AbstractVector)
    addnodes!(tree::AbstractTree, count::Integer)


"""
function addnodes! end

function addnodes!(tree::AbstractTree, nodenames::AbstractVector)
    all(map(name -> !_hasnode(tree, name), nodenames)) ||
        error("Some of nodes $nodenames already present in tree")
    return _addnodes!(tree, nodenames)
end

function addnodes!(tree::AbstractTree, count::Integer)
    return _addnodes!(tree, count)
end

"""
    deletenode!(tree::AbstractTree, nodename)


"""
function deletenode!(tree::AbstractTree, nodename)
    return _deletenode!(tree, nodename)
end

"""
    getnodenames(tree::AbstractTree)


"""
function getnodenames(tree::AbstractTree)
    return _getnodenames(tree)
end

"""
    hasnode(tree::AbstractTree, nodename)


"""
function hasnode(tree::AbstractTree, nodename)
    return _hasnode(tree, nodename)
end

"""
    getnode(tree::AbstractTree, nodename)


"""
function getnode(tree::AbstractTree, nodename)
    _hasnode(tree, nodename) ||
        error("Node $nodename does not exist")
    return _getnode(tree, nodename)
end

"""
    getbranchnames(tree::AbstractTree)


"""
function getbranchnames(tree::AbstractTree)
    return _getbranchnames(tree)
end

"""
    hasbranch(tree::AbstractTree, branchname)


"""
function hasbranch(tree::AbstractTree, branchname)
    return _hasbranch(tree, branchname)
end

"""
    getbranch(tree::AbstractTree, branchname)


"""
function getbranch(tree::AbstractTree, branchname)
    _hasbranch(tree, branchname) ||
        error("Branch $branchname does not exist")
    return _getbranch(tree, branchname)
end

"""
    hasrootheight(tree::AbstractTree)


"""
function hasrootheight(tree::AbstractTree)
    return _hasrootheight(tree)
end

"""
    getrootheight(tree::AbstractTree)


"""
function getrootheight(tree::AbstractTree)
    return _getrootheight(tree)
end

"""
    setrootheight!(tree::AbstractTree, height)


"""
function setrootheight!(tree::AbstractTree, height)
    return _setrootheight!(tree, height)
end

"""
    validate(tree::AbstractTree)


"""
function validate(tree::T) where {NL, BL, T <: AbstractTree{NL, BL}}
    nodes = _getnodes(tree)
    branches = _getbranches(tree)
    if !isempty(nodes) || !isempty(branches)
        # We need to validate the connections
        if Set(mapreduce(_getinbound, push!, nodefilter(_hasinbound, tree);
                         init = BL[])) != Set(keys(branches))
            warn("Inbound branches must exactly match Branch labels")
            return false
        end

        if Set(mapreduce(_getoutbounds, append!, nodeiter(tree);
                         init = BL[])) != Set(keys(branches))
            warn("Node outbound branches must exactly match Branch labels")
            return false
        end

        if !(mapreduce(_src, push!, branchiter(tree); init = NL[]) ⊆
             Set(keys(nodes)))
            warn("Branch sources must be node labels")
            return false
        end

        if !(mapreduce(_dst, push!, branchiter(tree); init = NL[]) ⊆
             Set(keys(nodes)))
            warn("Branch destinations must be node labels")
            return false
        end
    end

    return _validate(tree)
end


# AbstractNode methods
"""
    isleaf(node::AbstractNode)
    isleaf(tree::AbstractTree, nodename)


"""
function isleaf end

function isleaf(node::AbstractNode)
    return _isleaf(node)
end

function isleaf(tree::AbstractTree, nodename)
    return _isleaf(_getnode(tree, nodename))
end

"""
    isroot(node::AbstractNode)
    isroot(tree::AbstractTree, nodename)


"""
function isroot end

function isroot(node::AbstractNode)
    return _isroot(node)
end

function isroot(tree::AbstractTree, nodename)
    return _isroot(_getnode(tree, nodename))
end

"""
    isinternal(node::AbstractNode)
    isinternal(tree::AbstractTree, nodename)


"""
function isinternal end

function isinternal(node::AbstractNode)
    return _isinternal(node)
end

function isinternal(tree::AbstractTree, nodename)
    return _isinternal(_getnode(tree, nodename))
end

"""
    isunattached(node::AbstractNode)
    isunattached(tree::AbstractTree, nodename)


"""
function isunattached end

function isunattached(node::AbstractNode)
    return _isunattached(node)
end

function isunattached(tree::AbstractTree, nodename)
    return _isunattached(_getnode(tree, nodename))
end

"""
    indegree(node::AbstractNode)
    indegree(tree::AbstractTree, nodename)


"""
function indegree end

function indegree(node::AbstractNode)
    return _indegree(node)
end

function indegree(tree::AbstractTree, nodename)
    return _indegree(_getnode(tree, nodename))
end

"""
    outdegree(node::AbstractNode)
    outdegree(tree::AbstractTree, nodename)


"""
function outdegree end

function outdegree(node::AbstractNode)
    return _outdegree(node)
end

function outdegree(tree::AbstractTree, nodename)
    return _outdegree(_getnode(tree, nodename))
end

"""
    hasoutboundspace(node::AbstractNode)
    hasoutboundspace(tree::AbstractTree, nodename)

Does the node have space for an[other] outbound connection?
"""
function hasoutboundspace end

function hasoutboundspace(node::AbstractNode)
    return _hasoutboundspace(node)
end

function hasoutboundspace(tree::AbstractTree, nodename)
    return _hasoutboundspace(_getnode(tree, nodename))
end

"""
    hasinbound(node::AbstractNode)
    hasinbound(tree::AbstractTree, nodename)

Does the node have an inbound connection?
"""
function hasinbound end

function hasinbound(node::AbstractNode)
    return _hasinbound(node)
end

function hasinbound(tree::AbstractTree, nodename)
    return _hasinbound(_getnode(tree, nodename))
end

"""
    hasinboundspace(node::AbstractNode)
    hasinboundspace(tree::AbstractTree, nodename)

Does the node have space for an inbound connection?
"""
function hasinboundspace end

function hasinboundspace(node::AbstractNode)
    return _hasinboundspace(node)
end

function hasinboundspace(tree::AbstractTree, nodename)
    return _hasinboundspace(_getnode(tree, nodename))
end

"""
    getinbound(node::AbstractNode)
    getinbound(tree::AbstractTree, nodename)

return the name of the inbound branch to this node.
"""
function getinbound end

function getinbound(node::AbstractNode)
    return _getinbound(node)
end

function getinbound(tree::AbstractTree, nodename)
    return _getinbound(_getnode(tree, nodename))
end

"""
    getparent(tree::AbstractTree, nodename)

Return the name of the parent node for this node.
"""
function getparent(tree::AbstractTree, nodename)
    return src(tree, getinbound(tree, nodename))
end

"""
    getancestors(tree::AbstractTree, nodename)

Return the name of all of the nodes that are ancestral to this node.
"""
function getancestors(tree::AbstractTree, nodename)
    return _treepast(tree, nodename)[2][2:end]
end

"""
    getoutbounds(node::AbstractNode)
    getoutbounds(tree::AbstractTree, nodename)

Return the names of the outbound branches from this node.
"""
function getoutbounds end

function getoutbounds(node::AbstractNode)
    return _getoutbounds(node)
end

function getoutbounds(tree::AbstractTree, nodename)
    return _getoutbounds(_getnode(tree, nodename))
end

"""
    getchildren(tree::AbstractTree, nodename)

Return the name(s) of the child node(s) for this node.
"""
function getchildren(tree::AbstractTree, nodename)
    return map(branch -> dst(tree, branch), getoutbounds(tree, nodename))
end

"""
    getdescendants(tree::AbstractTree, nodename)

Return the names of all of the nodes that descend from this node.
"""
function getdescendants(tree::AbstractTree, nodename)
    return _treefuture(tree, nodename)[2][2:end]
end

"""
    hasheight(tree::AbstractTree, nodename)

Does the node have a height defined?
"""
function hasheight end

function hasheight(tree::AbstractTree, nodename)
    return _hasheight(tree, nodename) ||
        (_hasrootheight(tree) &&
         mapreduce(b -> haslength(tree, b), &, branchhistory(tree, nodename);
         init = _hasrootheight(tree)))
end

"""
    getheight(tree::AbstractTree, nodename)

Return the height of the node.
"""
function getheight(tree::AbstractTree, nodename)
    return _hasheight(tree, nodename) ? _getheight(tree, nodename) :
        mapreduce(b -> getlength(tree, b), +, branchhistory(tree, nodename);
                  init = getrootheight(tree))
end

"""
    setheight!(tree::AbstractTree, nodename, height)

Set the height of the node.
"""
function setheight!(tree::AbstractTree, nodename, height)
    return _setheight!(tree, nodename, height)
end


# AbstractBranch methods
"""
    src(branch::AbstractBranch)
    src(tree::AbstractTree, branchname)

Return the source node for this branch.
"""
function src end

function src(branch::AbstractBranch)
    return _src(branch)
end

function src(tree::AbstractTree, branchname)
    return _src(_getbranch(tree, branchname))
end

"""
    dst(branch::AbstractBranch)
    dst(tree::AbstractTree, branchname)

Return the destination node for this branch.
"""
function dst end

function dst(branch::AbstractBranch)
    return _dst(branch)
end

function dst(tree::AbstractTree, branchname)
    return _dst(_getbranch(tree, branchname))
end

"""
    Pair(branch::AbstractBranch)
    Pair(tree::AbstractTree, branchname)

Return a Pair containing the source and destination for this branch.
"""
function Pair end

function Pair(branch::AbstractBranch)
    return Pair(src(branch), dst(branch))
end

function Pair(tree::AbstractTree, branchname)
    return Pair(_getbranch(tree, branchname))
end

"""
    Tuple(branch::AbstractBranch)
    Tuple(tree::AbstractTree, branchname)

Return a Tuple containing the source and destination for this branch.
"""
function Tuple end

function Tuple(branch::AbstractBranch)
    return (src(branch), dst(branch))
end

function Tuple(tree::AbstractTree, branchname)
    return Tuple(_getbranch(tree, branchname))
end

"""
    getlength(branch::AbstractBranch)
    getlength(tree::AbstractTree, branchname)

Return the length of this branch.
"""
function getlength end

function getlength(branch::AbstractBranch)
    return _getlength(branch)
end

function getlength(tree::AbstractTree, branchname)
    return _getlength(_getbranch(tree, branchname))
end

"""
    changesrc!(tree::AbstractTree, branchname, source)

Change the source node for this branch.
"""
function changesrc!(tree::AbstractTree, branchname, source)
    _hasbranch(tree, branchname) ||
        error("Branch $branchname does not exist")
    _hasnode(tree, source) ||
        error("Node $source does not exist")
    branch = _getbranch(tree, branchname)
    oldsource = _src(branch)
    _setsrc!(branch, source)
    _deleteoutbound!(_getnode(tree, oldsource), branchname)
    _addoutbound!(_getnode(tree, source), branchname)
    return branchname
end

"""
    changedst!(tree::AbstractTree, branchname, destination)

Change the destination node for this node.
"""
function changedst!(tree::AbstractTree, branchname, destination)
    _hasbranch(tree, branchname) ||
        error("Branch $branchname does not exist")
    _hasnode(tree, destination) ||
        error("Node $destination does not exist")
    branch = _getbranch(tree, branchname)
    olddestination = _dst(branch)
    _setdst!(branch, destination)
    _deleteinbound!(_getnode(tree, olddestination), branchname)
    _setinbound!(_getnode(tree, destination), branchname)
    return branchname
end

"""
    resetleaves!(::AbstractTree)

Reset the leaf records to the current leaves, deleting all leaf records.
"""
function resetleaves!(tree::AbstractTree)
    return _resetleaves!(tree)
end


"""
    getleafnames(::AbstractTree)

Retrieve the leaf names from the tree.
"""
function getleafnames(tree::AbstractTree)
    return collect(_getleafnames(tree))
end

"""
    getleafinfo(::AbstractTree, label)

retrieve the leaf info for a leaf of the tree.
"""
function getleafinfo(tree::AbstractTree, label)
    return _getleafinfo(tree, label)
end
function getleafinfo(tree::AbstractTree)
    return _getleafinfo(tree)
end
function leafinfotype(tree::AbstractTree)
    return _leafinfotype(tree)
end

"""
    setleafinfo!(::AbstractTree, table)

Set the leaf info for the leaves of the tree.
"""
function setleafinfo!(tree::AbstractTree, table)
    return _setleafinfo!(tree, table)
end

"""
    getnoderecord(::AbstractTree, label)

retrieve the node record for a leaf of the tree.
"""
function getnoderecord(tree::AbstractTree, label)
    return _getnoderecord(tree, label)
end

"""
    setnoderecord(::AbstractTree, label, value)

Set the node record for a node of the tree.
"""
function setnoderecord!(tree::AbstractTree, label, value)
    return _setnoderecord!(tree, label, value)
end

"""
    nleaves(::AbstractTree)

Count the number of leaves in the tree.
"""
nleaves(tree::AbstractTree) = _nleaves(tree)
