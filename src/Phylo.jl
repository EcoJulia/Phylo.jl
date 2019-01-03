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

abstract type AbstractNode end
abstract type AbstractBranch end
abstract type AbstractTree{NodeLabel, BranchLabel} end
abstract type AbstractBranchTree{NL, BL} <: AbstractTree{NL, BL} end
export AbstractNode, AbstractBranch, AbstractTree

"""
    Phylo.API submodule

The `Phylo.API` submodule should be `import`ed if you want to
create a new phylogeny, node or branch subtype. Otherwise it can be
ignored.
"""
module API
include("API.jl")
# AbstractTree methods
export _ntrees, _addbranch!, _deletebranch!, _branch!, _setbranch!
export _addnode!, _addnodes!, _deletenode!, _setnode!
export _getnodenames, _hasnode, _getnode, _getnodes
export _getbranchnames, _hasbranch, _getbranch, _getbranches
export _hasrootheight, _getrootheight, _setrootheight!, _clearrootheight!
export _nodetype, _branchtype
export _extractnode, _extractbranch
export _extractnodename, _extractbranchname
export _getleafinfo, _setleafinfo!, _leafinfotype
export _getnoderecord, _setnoderecord!
export _hasheight, _getheight, _setheight!
export _hasparent, _getparent, _getancestors
export _haschildren, _getchildren, _getdescendants
export _validate
export _getleafnames, _resetleaves!, _nleaves

# AbstractNode methods
export _isleaf, _isroot, _isinternal, _isunattached
export _indegree, _hasinbound, _getinbound, _setinbound!, _deleteoutbound!
export _outdegree, _getoutbounds, _addoutbound!, _deleteoutbound!
export _hasoutboundspace, _hasinboundspace

# AbstractBranch methods
export _src, _dst, _getlength
export _setsrc!, _setdst!

# Label names
export _newnodelabel, _newbranchlabel

end

include("Interface.jl")
# AbstractTree methods
export ntrees, nodetype, branchtype, nodenametype, branchnametype
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
export hasoutboundspace, hasinboundspace
export getleafnames, resetleaves, nleaves
export getleafinfo, setleafinfo!, leafinfotype
export getnoderecord, setnoderecord!
export hasheight, getheight, setheight!

# AbstractTree / AbstractBranch methods
export src, dst, getlength
export changesrc!, changedst!

include("Info.jl")
export LeafInfo

include("Branch.jl")
export Branch

include("Node.jl")
export BinaryNode, Node

include("Tree.jl")
export BinaryTree, NamedBinaryTree, NamedTree
export PolytomousTree, NamedPolytomousTree
export hasrootheight, getrootheight, setrootheight!, clearrootheight!
export getrootdistance

include("routes.jl")
export branchhistory, branchroute, nodehistory, noderoute
export distance, distances, heighttoroot, heightstoroot

# Iterator methods expanded
include("Iterators.jl")
export nodeiter, nodefilter, nodenameiter, nodenamefilter,
    branchiter, branchfilter, branchnameiter, branchnamefilter

# A set of multiple trees
include("TreeSet.jl")
export TreeSet, treeiter, treenameiter, treeinfoiter

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

# Plot recipes
include("plot.jl")

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
