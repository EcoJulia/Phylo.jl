__precompile__()

"""
    Phylo package

The `Phylo` package provides a standard abstract interface to
phylogenetic trees, by defining `AbstractNode`, `AbstractBranch` and
`AbstractTree` supertypes, and methods to interface to them. It also
provides (through the `Phylo.API` submodule) methods to
(re)define to write your own phylogenetic type in a way that will
interact cleanly with other phylogenetic packages. Finally, it provides
a simple phylogenetics type.
"""
module Phylo

using Compat
@compat abstract type AbstractNode end
@compat abstract type AbstractBranch end
@compat abstract type AbstractTree{NodeLabel, BranchLabel} end
@compat abstract type AbstractInfo end
export AbstractNode, AbstractBranch, AbstractTree, AbstractInfo

"""
    Phylo.API submodule

The `Phylo.API` submodule should be `import[all]`ed if you want to
create a new phylogeny, node or branch subtype. Otherwise it can be
ignored.
"""
module API
include("API.jl")
# AbstractTree methods
export _addbranch!, _deletebranch!, _branch!
export _addnode!, _addnodes!, _deletenode!
export _getnodenames, _hasnode, _getnode, _getnodes
export _getbranchnames, _hasbranch, _getbranch, _getbranches
export _hasrootheight, _getrootheight, _setrootheight!
export _nodetype, _branchtype
export _extractnode, _extractbranch
export _extractnodename, _extractbranchname
export _getleafinfo, _setleafinfo!, _getnoderecord, _setnoderecord!
export _hasheight, _getheight, _setheight!
export _hasparent, _getparent, _getancestors
export _haschildren, _getchildren, _getdescendants
export _validate

# AbstractNode methods
export _isleaf, _isroot, _isinternal, _isunattached
export _indegree, _outdegree, _hasinbound, _getinbound, _getoutbounds
export _hasoutboundspace, _hasinboundspace

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
export getnodenames, hasnode, getnode, getnodes
export getbranchnames, hasbranch, getbranch, getbranches
export hasrootheight, getrootheight, setrootheight!
export hasparent, getparent, getancestors
export haschildren, getchildren, getdescendants
export validate

# AbstractTree / AbstractNode methods
export isleaf, isroot, isinternal, isunattached
export indegree, outdegree, hasinbound, getinbound, getoutbounds
export getleafnames
export getleafinfo, setleafinfo!, getnoderecord, setnoderecord!
export hasheight, getheight, setheight!

# AbstractTree / AbstractBranch methods
export getsource, gettarget, getlength
export changesource!, changetarget!

include("Info.jl")
export LeafInfo

include("Branch.jl")
export Branch

include("Node.jl")
export BinaryNode

include("Tree.jl")
export BinaryTree, NamedTree
export hasrootheight, getrootheight, setrootheight!, clearrootheight!
export getrootdistance

include("routes.jl")
export branchhistory, branchroute, nodehistory, noderoute
export distance, distances, heighttoroot, heightstoroot

# Iterator methods expanded
include("Iterators.jl")
export NodeIterator, NodeNameIterator, BranchIterator

# Random tree generator
include("rand.jl")
export Nonultrametric, Ultrametric

# Display methods expanded
include("show.jl")

end # module
