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
import Compat: IteratorSize, IteratorEltype

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

abstract type AbstractTree{TT <: TreeType, RT <: Rootedness, NL,
                           N <: AbstractNode{RT, NL},
                           B <: AbstractBranch{RT, NL}} end
export AbstractTree

@enum TraversalOrder preorder inorder postorder breadthfirst
export TraversalOrder

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
export _getnodenames, _hasnode, _getnode, _getnodes
export _getbranchnames, _getbranchname, _hasbranch, _getbranch, _getbranches
export _hasrootheight, _getrootheight, _setrootheight!, _clearrootheight!
export _getleafinfo, _setleafinfo!, _leafinfotype
export _getnoderecord, _setnoderecord!
export _hasheight, _getheight, _setheight!
export _hasparent, _getparent, _getancestors
export _haschildren, _getchildren, _getdescendants
export _validate, _traversal
export _getleafnames, _getleaves, _resetleaves!, _nleaves, _nnodes, _nbranches

# AbstractNode methods
export _isleaf, _isroot, _isinternal, _isunattached
export _indegree, _hasinboundspace, _outdegree, _hasoutboundspace, _degree
export _hasinbound, _getinbound, _addinbound!, _removeinbound!
export _getoutbounds, _addoutbound!, _removeoutbound!
export _getconnections, _addconnection!, _removeconnection!

# AbstractBranch methods
export _src, _dst, _getlength

# Label names
export _newnodelabel, _newbranchlabel

end

include("Interface.jl")
# AbstractTree methods
export ntrees, gettrees, nroots, getroots, getroot
export treenametype, gettreenames, getonetree, gettreename
export nodetype, branchtype, nodenametype, branchnametype
export createbranch!, deletebranch!, branch!
export createnode!, createnodes!, deletenode!
export getnodenames, getnodename, hasnode, getnode, getnodes
export getbranchnames, getbranchname, hasbranch, getbranch, getbranches
export hasrootheight, getrootheight, setrootheight!
export hasparent, getparent, getancestors
export haschildren, getchildren, getdescendants
export validate, traversal

# AbstractTree / AbstractNode methods
export isleaf, isroot, isinternal, isunattached
export indegree, outdegree, hasinbound, getinbound, getoutbounds
export hasoutboundspace, hasinboundspace
export getleafnames, getleaves, resetleaves, nleaves, nnodes
export getleafinfo, setleafinfo!, leafinfotype
export getnoderecord, setnoderecord!
export hasheight, getheight, setheight!

# AbstractTree / AbstractBranch methods
export src, dst, getlength
export hasrootheight, getrootheight, setrootheight!, clearrootheight!
export getrootdistance

#include("Info.jl")
#export LeafInfo

include("Branch.jl")
export Branch

include("Node.jl")
export BinaryNode, Node

include("Tree.jl")
export BinaryTree, NamedBinaryTree, NamedTree
export PolytomousTree, NamedPolytomousTree

#include("LinkTree.jl")
#export LinkBranch, LinkNode, LinkTree

#include("MetaTree.jl")
#export UnrootedMetaTree, RootedMetaTree, SimpleNode, SimpleBranch
#export RootedTree, ManyRootTree, UnrootedTree

include("routes.jl")
export branchhistory, branchfuture, branchroute
export nodehistory, nodefuture, noderoute
export distance, distances, heighttoroot, heightstoroot

# Iterator methods expanded
include("Iterators.jl")
export nodeiter, nodefilter, nodenameiter, nodenamefilter,
    branchiter, branchfilter, branchnameiter, branchnamefilter

# A set of multiple trees
include("TreeSet.jl")
export TreeSet, gettreeinfo

# Random tree generator
include("rand.jl")
export Nonultrametric, Ultrametric

# Read Newick Tree
include("newick.jl")
export parsenewick, parsenexus

# Display methods expanded
include("show.jl")

# Method for trimming trees
include("trim.jl")
export droptips!, keeptips!

# Plot recipes, only works on Julia v0.7 and up
@static if VERSION > v"0.7.0-"
    include("plot.jl")
end

# Path into package
path(path...; dir::String = "test") = joinpath(@__DIR__, "..", dir, path...)

using Requires
@static if VERSION < v"0.7.0-"
    @require RCall begin
        println("Creating Phylo RCall interface...")
        include("rcall.jl")
    end
else
    function __init__()
        @require RCall="6f49c342-dc21-5d91-9882-a32aef131414" begin
            println("Creating Phylo RCall interface...")
            include("rcall.jl")
        end
    end
end

end # module
