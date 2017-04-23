__precompile__()

"""
    AbstractPhylo package

The `AbstractPhylo` package provides a standard abstract interface to
phylogenetic trees, by defining `AbstractNode`, `AbstractBranch` and
`AbstractTree` supertypes, and methods to interface to them. It also
provides (through the `AbstractPhylo.API` submodule) methods to
(re)define to write your own phylogenetic type in a way that will
interact cleanly with other phylogenetic packages.
"""
module AbstractPhylo

"""
    AbstractPhylo.API submodule

The `AbstractPhylo.API` submodule should be `import[all]`ed if you want to
create a new phylogeny, node or branch subtype. Otherwise it can be
ignored.
"""
module API
include("API.jl")
export AbstractNode, AbstractBranch, AbstractTree
# AbstractTree methods
export _addbranch!, _deletebranch!, _branch!
export _addnode!, _addnodes!, _deletenode!
export _getnodenames, _hasnode, _getnode
export _getbranchnames, _hasbranch, _getbranch
export _hasrootheight, _getrootheight, _setrootheight!
export _nodetype, _branchtype, _extractnode, _extractbranch
export _validate

# AbstractNode methods
export _isleaf, _isroot, _isinternal, _isunattached
export _indegree, _outdegree, _getinbound, _getoutbounds
export _hasheight, _getheight, _setheight!

# AbstractBranch methods
export _getsource, _gettarget, _getlength
export _setsource!, _settarget!

# Label names
export _newnodelabel, _newbranchlabel

end

include("Interface.jl")
# AbstractTree methods
export addbranch!, deletebranch!, branch!
export addnode!, addnodes!, deletenode!
export getnodenames, hasnode, getnode
export getbranchnames, hasbranch, getbranch
export hasrootheight, getrootheight, setrootheight!
export getparent, getchildren, getancestors, getdescendants
export validate

# AbstractTree / AbstractNode methods
export isleaf, isroot, isinternal, isunattached
export indegree, outdegree, getinbound, getoutbounds
export hasheight, getheight, setheight!

# AbstractTree / AbstractBranch methods
export getsource, gettarget, getlength
export changesource!, changetarget!

include("Iterators.jl")
export NodeIterator, BranchIterator

end # module
