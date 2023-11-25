__precompile__()

"""
    Phylo package

The `Phylo` package provides some simple phylogenetics types (e.g. NamedTree)
to interface to the `Diversity` package for measuring phylogenetic diversity.
It also provides an interface to `R` for copying trees to and from that
language and can read newick and nexus tree files (including `TreeSet`s that
contain multiple trees).

Finally it also provides a standard abstract interface to
phylogenetic trees, by defining `AbstractNode`, `AbstractBranch` and
`AbstractTree` supertypes, and methods to interface to them. It also
provides (through the `Phylo.API` submodule) methods to
(re)define to write your own phylogenetic type in a way that will
interact cleanly with other phylogenetic packages.
"""
module Phylo

import Base: Pair, Tuple, show, eltype, length, getindex
import Graphs: src, dst, indegree, outdegree, degree
abstract type Rootedness end
struct Unrooted <: Rootedness end
abstract type Rooted <: Rootedness end
struct OneRoot <: Rooted end
struct ManyRoots <: Rooted end
export Unrooted, OneRoot, ManyRoots

abstract type TreeType end
struct OneTree <: TreeType end
struct ManyTrees <: TreeType end
export OneTree, ManyTrees

abstract type AbstractNode{RootType <: Rootedness, NodeLabel} end
abstract type AbstractBranch{RootType <: Rootedness, NodeLabel} end

using Distances
abstract type AbstractTree{TT <: TreeType, RT <: Rootedness, NL,
                           N <: AbstractNode{RT, NL},
                           B <: AbstractBranch{RT, NL}} <: Distances.UnionMetric
end

export AbstractTree

@enum TraversalOrder anyorder preorder inorder postorder breadthfirst
export TraversalOrder, anyorder, preorder, inorder, postorder, breadthfirst

"""
    Phylo.API submodule

The `Phylo.API` submodule should be `import`ed if you want to
create a new phylogeny, node or branch subtype. Otherwise it can be
ignored.
"""
module API
include("API.jl")
# AbstractTree methods
export _ntrees, _gettrees, _nroots, _getroots, _getroot
export _treenametype, _gettreenames, _gettree, _gettreename
export _createbranch!, _deletebranch!, _createnode!, _deletenode!
export _getnodenames, _getnodename, _hasnode, _getnode, _getnodes
export _getbranchnames, _getbranchname, _hasbranch, _getbranch, _getbranches
export _hasrootheight, _getrootheight, _setrootheight!, _clearrootheight!
export _getleafinfo, _setleafinfo!, _leafinfotype, _gettreeinfo
export _getnodedata, _setnodedata!, _nodedatatype
export _getbranchdata, _setbranchdata!, _branchdatatype
export _hasheight, _getheight, _setheight!
export _hasparent, _getparent, _getancestors
export _haschildren, _getchildren, _getdescendants
export _validate!, _traversal, _branchdims
export _getleafnames, _getleaves, _resetleaves!, _nleaves, _nnodes, _nbranches
export HoldsNodeData, MatchTreeNameType

# AbstractNode methods
export _isleaf, _isroot, _isinternal, _isunattached
export _indegree, _hasinboundspace, _outdegree, _hasoutboundspace, _hasspace, _degree
export _hasinbound, _getinbound, _addinbound!, _removeinbound!
export _getoutbounds, _addoutbound!, _removeoutbound!
export _getconnections, _addconnection!, _removeconnection!
export MatchNodeType, MatchNodeTypes, PreferNodeObjects, _prefernodeobjects

# AbstractBranch methods
export _src, _dst, _getlength, _haslength, _conn, _conns
export MatchBranchType, PreferBranchObjects, _preferbranchobjects
export MatchBranchNodeType

# Label names
export _newnodelabel, _newbranchlabel

end

include("Interface.jl")
# AbstractTree methods
export ntrees, gettrees, nroots, getroots, getroot, gettree
export treenametype, gettreenames, gettreename #, getonetree #unimplemented
export treetype, roottype, nodetype, nodedatatype, nodenametype
export branchtype, branchdatatype, branchnametype
export createbranch!, deletebranch!
export createnode!, createnodes!, deletenode!
export getnodenames, getnodename, hasnode, getnode, getnodes, nnodes
export getleafnames, getleaves, nleaves, getinternalnodes, ninternal
export getbranchnames, getbranchname, hasbranch, getbranch, getbranches, nbranches
export hasrootheight, getrootheight, setrootheight!
export validate!, traversal, branchdims

@deprecate addnode! createnode!
@deprecate addnodes! createnodes!
@deprecate addbranch! createbranch!
@deprecate(branch!(tree, source, length = missing;
                   destination = _newnodelabel(tree),
                   branchname = _newbranchlabel(tree)),
           createbranch!(tree, source, createnode!(tree, destination),
                         length; name = branchname))

# AbstractTree / AbstractNode methods
export isleaf, isroot, isinternal, isunattached
export degree, indegree, outdegree, hasinbound, getconnections, getinbound, getoutbounds
export hasoutboundspace, hasinboundspace
export getleafinfo, setleafinfo!, leafinfotype
export getnodedata, setnodedata!
export getparent, getancestors #, hasparent # unimplemented
export getchildren, getdescendants #, haschildren # unimplemented
export getsiblings
export hasheight, getheight, setheight!

# AbstractTree / AbstractBranch methods
export src, dst, getlength, haslength, conn, conns
export hasrootheight, getrootheight, setrootheight! #, clearrootheight! #unimplemented
export getbranchdata, setbranchdata!
# export getrootdistance # unimplemented

@deprecate getnoderecord getnodedata
@deprecate setnoderecord! setnodedata!
@deprecate getbranchrecord getbranchdata
@deprecate setbranchrecord! setbranchdata!
@deprecate getbranchinfo getbranchdata
@deprecate setbranchinfo! setbranchdata!

include("Branch.jl")
export Branch

include("Node.jl")
export BinaryNode, Node

include("Tree.jl")
export BinaryTree, NamedBinaryTree, NamedTree
export PolytomousTree, NamedPolytomousTree

include("LinkTree.jl")
export LinkBranch, LinkNode, LinkTree
export RootedTree, ManyRootTree, UnrootedTree

include("routes.jl")
export branchhistory, branchfuture, branchroute
export nodehistory, nodefuture, noderoute

# Iterator methods expanded
include("Iterators.jl")
export nodeiter, nodefilter, nodenameiter, nodenamefilter,
    branchiter, branchfilter, branchnameiter, branchnamefilter
@deprecate treeiter gettrees
@deprecate treenameiter gettreenames

# A set of multiple trees
include("TreeSet.jl")
export TreeSet, gettreeinfo

# Random tree generator
include("rand.jl")
export Nonultrametric, Ultrametric
export BrownianTrait, DiscreteTrait, SymmetricDiscreteTrait

# Read Newick Tree
include("newick.jl")
export parsenewick, parsenexus

# Display methods expanded
include("show.jl")

# Method for trimming trees
include("trim.jl")
export droptips!, keeptips!

# Plot recipes
include("plot.jl")
export map_depthfirst

# Metrics from the tree
include("metrics.jl")
export mrca, nodeheights
export distance, distances, heighttoroot, heightstoroot

# Path into package
path(path...; dir::String = "test") = joinpath(@__DIR__, "..", dir, path...)

# This symbol is only defined on Julia versions that support extensions
if !isdefined(Base, :get_extension)
    using Requires
end

@static if !isdefined(Base, :get_extension)
function __init__()
    @require RCall="6f49c342-dc21-5d91-9882-a32aef131414" include("../ext/PhyloRCallExt.jl")
end
end

end # module
