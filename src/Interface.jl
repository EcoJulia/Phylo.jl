using Compat
using AbstractPhylo.API

# AbstractTree methods
"""
    addbranch!(tree::AbstractTree)
    addbranch!(tree::AbstractTree, branchname)


"""
function addbranch!(tree::AbstractTree, branchname = _newbranchlabel(tree))
    return _addbranch!(tree, branchname)
end

"""
    deletebranch!(tree::AbstractTree, branchname)


"""
function deletebranch!(tree::AbstractTree, branchname)
    return _deletebranch!(tree, branchname)
end

"""
    branch!(tree::AbstractTree)
    branch!(tree::AbstractTree, nodename)
    branch!(tree::AbstractTree, nodename, branchname)


"""
function branch!(tree::AbstractTree,
                 nodename = _newnodelabel(tree),
                 branchname = _newbranchlabel(tree))
    return _branch!(tree, nodename, branchname)
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
    hasbranch(tree::AbstractTree, nodename)


"""
function hasbranch(tree::AbstractTree, nodename)
    return _hasbranch(tree, nodename)
end

"""
    getbranch(tree::AbstractTree, nodename)


"""
function getbranch(tree::AbstractTree, nodename)
    return _getbranch(tree, nodename)
end

"""
    hasrootheight


"""
function hasrootheight(tree::AbstractTree)

end

"""
    getrootheight


"""
function getrootheight(tree::AbstractTree)

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
function validate(tree::AbstractTree)
    nodes = getnodes(tree)
    branches = getbranches(tree)
    if !isempty(nodes) || !isempty(branches)
        # We need to validate the connections
        if Set(mapreduce(_getinbound, push!, BL[],
                         Compat.Iterators.filter(_hasinbound, values(nodes)))) !=
                             Set(keys(branches))
            warn("Inbound branches must exactly match Branch labels")
            return false
        end
        
        if Set(mapreduce(_getoutbounds, append!, BL[], values(nodes))) !=
            Set(keys(branches))
            warn("Node outbound branches must exactly match Branch labels")
            return false
        end
        
        if !(mapreduce(_getsource, push!, NL[], values(branches)) ⊆
             Set(keys(nodes)))
            warn("Branch sources must be node labels")
            return false
        end

        if !(mapreduce(_gettarget, push!, NL[], values(branches)) ⊆
             Set(keys(nodes)))
            warn("Branch targets must be node labels")
            return false
        end

        if length(findroots(tree) ∪ findunattacheds(tree)) == 0
            warn("This tree has no roots")
            return false
        end

        if length(findleaves(tree) ∪ findunattacheds(tree)) == 0
            warn("This tree has no leaves")
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
    return _setheight(tree, nodename, height)
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

