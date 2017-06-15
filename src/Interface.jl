using Compat
using Phylo.API

getnodes(tree::AbstractTree) = _getnodes(tree)
getbranches(tree::AbstractTree) = _getbranches(tree)

# AbstractTree methods
"""
    addbranch!(tree::AbstractTree, source, target[, length::Float64];
               branchname = _newbranchlabel(tree))

Add a branch from `source` to `target` on `tree`.
"""
function addbranch!(tree::AbstractTree, source, target, length::Float64 = NaN;
                    branchname = _newbranchlabel(tree))
    _hasnode(tree, source) && hasoutboundspace(tree, source) ||
        error("Tree does not have an available source node called $source")
    _hasnode(tree, target) && !hasinbound(tree, target) ||
        error("Tree does not have an available target node called $target")
    target != source || error("Branch must connect different nodes")
    _hasbranch(tree, branchname) &&
        error("Tree already has a branch called $branchname")
    
    return _addbranch!(tree, source, target, length, branchname)
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
    branch!(tree::AbstractTree, source[, length]; target)
    branch!(tree::AbstractTree, source[, length]; target, branchname)

Branch from a source node `source` and create a target node `target`.
"""
function branch!(tree::AbstractTree, source, length::Float64 = NaN;
                 target = _newnodelabel(tree),
                 branchname = _newbranchlabel(tree))
    _hasnode(tree, source) ||
        error("Node $source not present in tree")
    !_hasnode(tree, target) ||
        error("Node $target already present in tree")
    _hasoutboundspace(_getnode(tree, source)) ||
        error("Node $source has no space to add branches")
    
    return _branch!(tree, source, length, target, branchname)
end

"""
    addnode!(tree::AbstractTree)
    addnode!(tree::AbstractTree, nodename)


"""
function addnode!(tree::AbstractTree, nodename = _newnodelabel(tree))
    return _addnode!(tree, nodename)
end

"""
    addnodes!(tree::AbstractTree, nodenames::AbstractVector)
    addnodes!(tree::AbstractTree, count::Integer)


"""
function addnodes! end

function addnodes!(tree::AbstractTree, nodenames::AbstractVector)
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
function validate{NL, BL}(tree::AbstractTree{NL, BL})
    nodes = _getnodes(tree)
    branches = _getbranches(tree)
    if !isempty(nodes) || !isempty(branches)
        # We need to validate the connections
        if Set(mapreduce(_getinbound, push!, BL[],
                         NodeIterator(tree, _hasinbound))) !=
                             Set(keys(branches))
            warn("Inbound branches must exactly match Branch labels")
            return false
        end
        
        if Set(mapreduce(_getoutbounds, append!, BL[], NodeIterator(tree))) !=
            Set(keys(branches))
            warn("Node outbound branches must exactly match Branch labels")
            return false
        end
        
        if !(mapreduce(_getsource, push!, NL[], BranchIterator(tree)) ⊆
             Set(keys(nodes)))
            warn("Branch sources must be node labels")
            return false
        end

        if !(mapreduce(_gettarget, push!, NL[], BranchIterator(tree)) ⊆
             Set(keys(nodes)))
            warn("Branch targets must be node labels")
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


"""
function getinbound end

function getinbound(node::AbstractNode)
    return _getinbound(node)
end

function getinbound(tree::AbstractTree, nodename)
    return _getinbound(_getnode(tree, nodename))
end

"""
    getoutbounds(node::AbstractNode)
    getoutbounds(tree::AbstractTree, nodename)


"""
function getoutbounds end

function getoutbounds(node::AbstractNode)
    return _getoutbounds(node)
end

function getoutbounds(tree::AbstractTree, nodename)
    return _getoutbounds(_getnode(tree, nodename))
end

"""
    hasheight(tree::AbstractTree, nodename)


"""
function hasheight end

function hasheight(tree::AbstractTree, nodename)
    return _hasheight(tree, nodename)
end

"""
    getheight(tree::AbstractTree, nodename)


"""
function getheight end

function getheight(tree::AbstractTree, nodename)
    return _getheight(tree, nodename)
end

"""
    setheight!(tree::AbstractTree, nodename, height)


"""
function setheight!(tree::AbstractTree, nodename, height)
    return _setheight!(tree, nodename, height)
end


# AbstractBranch methods
"""
    getsource(branch::AbstractBranch)
    getsource(tree::AbstractTree, branchname)


"""
function getsource end

function getsource(branch::AbstractBranch)
    return _getsource(branch)
end

function getsource(tree::AbstractTree, branchname)
    return _getsource(_getbranch(tree, branchname))
end

"""
    gettarget(branch::AbstractBranch)
    gettarget(tree::AbstractTree, branchname)


"""
function gettarget end

function gettarget(branch::AbstractBranch)
    return _gettarget(branch)
end

function gettarget(tree::AbstractTree, branchname)
    return _gettarget(_getbranch(tree, branchname))
end

"""
    getlength(branch::AbstractBranch)
    getlength(tree::AbstractTree, branchname)


"""
function getlength end

function getlength(branch::AbstractBranch)
    return _getlength(branch)
end

function getlength(tree::AbstractTree, branchname)
    return _getlength(_getbranch(tree, branchname))
end

"""
    changesource!(tree::AbstractTree, branchname, source)


"""
function changesource!(tree::AbstractTree, branchname, source)
    _hasbranch(tree, branchname) ||
        error("Branch $branchname does not exist")
    _hasnode(tree, source) ||
        error("Node $source does not exist")
    branch = _getbranch(tree, branchname)
    oldsource = _getsource(branch)
    _setsource!(branch, source)
    _deleteoutbound!(tree, oldsource, branchname)
    _addoutbound!(tree, source, branchname)
    return branchname
end

"""
    changetarget!(tree::AbstractTree, branchname, target)


"""
function changetarget!(tree::AbstractTree, branchname, target)
    _hasbranch(tree, branchname) ||
        error("Branch $branchname does not exist")
    _hasnode(tree, target) ||
        error("Node $target does not exist")
    branch = _getbranch(tree, branchname)
    oldtarget = _gettarget(branch)
    _settarget!(branch, target)
    _deleteinbound!(tree, oldtarget, branchname)
    _setinbound!(tree, target, branchname)
    return branchname
end


"""
    getleafnames(::AbstractTree)

retrieve the leaf names from the tree.
"""
function getleafnames(tree::AbstractTree)
    return collect(_getleafnames(tree))
end

"""
    getleafrecords(::AbstractTree)

retrieve the Dict containing the leaf records of the tree.
"""
function getleafrecords(tree::AbstractTree)
    return _getleafrecords(tree)
end

"""
    getleafrecord(::AbstractTree, label)

retrieve the leaf record for a leaf of the tree.
"""
function getleafrecord(tree::AbstractTree, label)
    return _getleafrecord(tree, label)
end

"""
    setleafrecord(::AbstractTree, label, value)

Set the leaf record for a leaf of the tree.
"""
function setleafrecord!(tree::AbstractTree, label, value)
    return _setleafrecord!(tree, label, value)
end

"""
    getnoderecords(::AbstractTree)

retrieve the Dict containing the node records of the tree.
"""
function getnoderecords(tree::AbstractTree)
    return _getnoderecords(tree)
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
